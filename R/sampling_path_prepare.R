sampling_path_prepare <- function(
  tag,
  likelihood,
  movement,
  thr_likelihood,
  thr_gs,
  quiet = FALSE,
  profile = FALSE
) {
  GeoPressureR::tag_assert(tag, "setmap")

  likelihood <- GeoPressureR::tag2likelihood(tag, likelihood)

  # Build the light likelihood maps and the structural mask once.
  map <- GeoPressureR::tag_prepare_likelihood(
    tag,
    likelihood = likelihood,
    thr_likelihood = thr_likelihood,
    thr_gs = Inf,
    quiet = quiet
  )
  g <- GeoPressureR::map_expand(map$extent, map$scale)
  ncell <- prod(g$dim)
  struct_mask <- if ("mask_water" %in% names(map)) {
    !as.vector(map$mask_water)
  } else {
    rep(TRUE, ncell)
  }

  # Use the light HDR as a fixed state space. Movement affects probabilities
  # inside this support but does not remove or recover likelihood cells.
  lk <- Map(
    function(lk_prob, lk_keep, stap_id) {
      finite_keep <- lk_keep & is.finite(lk_prob) & lk_prob > 0 & struct_mask
      idx <- which(finite_keep)

      list(
        stap_id = stap_id,
        idx = idx,
        log_prob = log(lk_prob[idx])
      )
    },
    map$lk_norm,
    map$lk_mask,
    map$stap$stap_id[map$stap$include]
  )

  # A likelihood concentrated on the outer grid edge can indicate that the map
  # extent truncated the light surface before sampler support was constructed.
  map_boundary <- matrix(FALSE, nrow = g$dim[1], ncol = g$dim[2])
  map_boundary[c(1L, g$dim[1]), ] <- TRUE
  map_boundary[, c(1L, g$dim[2])] <- TRUE
  included_stap <- which(map$stap$include)
  unknown_stap <- is.na(map$stap$known_lat[included_stap]) |
    is.na(map$stap$known_lon[included_stap])
  boundary_mass <- vapply(
    map$lk_norm,
    function(x) sum(x[map_boundary], na.rm = TRUE),
    numeric(1)
  )
  boundary_mode <- vapply(
    map$lk_norm,
    function(x) map_boundary[which.max(x)],
    logical(1)
  )
  n_boundary_mode <- sum(boundary_mode[unknown_stap])
  n_boundary_mass <- sum(boundary_mass[unknown_stap] > 0.01)

  if (!quiet && (n_boundary_mode > 0L || n_boundary_mass > 0L)) {
    cli::cli_alert_warning(
      "The configured map extent may truncate light likelihoods: {n_boundary_mode} unknown-period mode{?s} lie on the grid edge and {n_boundary_mass} period{?s} place more than 1% of likelihood mass there."
    )
  }

  support_idx <- sort(unique(unlist(lapply(lk, `[[`, "idx"))))

  if (!quiet) {
    cli::cli_progress_done()
    cli::cli_progress_step(
      "Precompute movement kernel",
      msg_done = "Movement kernel ready"
    )
  }

  # Precompute the local movement footprints used by the Gibbs transitions.
  # Great-circle distances keep the hard movement support symmetric across rows.
  nrow <- g$dim[1]
  ncol <- g$dim[2]
  sampling_path_movement_validate(movement)
  move_probability <- movement$move_stay(0:(length(lk) - 1L))
  if (
    !all(
      is.finite(move_probability) &
        move_probability >= 0 &
        move_probability <= 1
    )
  ) {
    cli::cli_abort(
      "The movement model has `move_stay()` values outside [0, 1] for this track."
    )
  }
  ddeg_rad <- (1 / g$scale) * pi / 180
  earth_radius <- 6378137
  step_hours <- 24
  max_dist_m <- thr_gs * step_hours * 1000
  dy_m <- earth_radius * ddeg_rad
  max_row <- min(as.integer(ceiling(max_dist_m / dy_m)), nrow - 1L)
  drow_grid <- seq.int(-max_row, max_row)
  dlat_rad <- stats::median(diff(g$lat)) * pi / 180

  cell_row <- ((seq_len(ncell) - 1L) %% nrow) + 1L
  cell_col <- ((seq_len(ncell) - 1L) %/% nrow) + 1L
  movement_eps <- movement$approx_eps %||% 1e-6
  speed_model <- movement
  speed_model$move_stay <- NULL
  speed_model$move_stay_parameters <- NULL
  speed_model$approx_eps <- NULL
  kernels <- vector("list", nrow)
  drow <- vector("list", nrow)
  dcol <- vector("list", nrow)

  footprint_seconds <- system.time(
    for (i in seq_len(nrow)) {
      candidate_row <- i + drow_grid
      candidate_row <- candidate_row[
        candidate_row >= 1L & candidate_row <= nrow
      ]
      dx_m <- earth_radius *
        min(cos(g$lat[c(i, candidate_row)] * pi / 180)) *
        ddeg_rad
      max_col <- min(
        if (dx_m > 0) as.integer(ceiling(max_dist_m / dx_m)) else 0L,
        ncol - 1L
      )
      dcol_grid <- seq.int(-max_col, max_col)
      lat_from <- g$lat[i] * pi / 180
      lat_to <- lat_from + drow_grid * dlat_rad
      sin_dlat2 <- sin((lat_to - lat_from) / 2)^2
      sin_dlon2 <- sin(dcol_grid * ddeg_rad / 2)^2
      a <- matrix(
        sin_dlat2,
        nrow = length(drow_grid),
        ncol = length(dcol_grid)
      ) +
        outer(cos(lat_from) * cos(lat_to), sin_dlon2, `*`)
      dist_m <- 2 * earth_radius * atan2(sqrt(pmin(a, 1)), sqrt(pmax(1 - a, 0)))
      idx <- which(dist_m <= max_dist_m, arr.ind = TRUE)
      dr <- drow_grid[idx[, 1]]
      dc <- dcol_grid[idx[, 2]]
      dist_keep <- dist_m[idx]
      speed_prob <- GeoPressureR::speed2prob(
        dist_keep / 1000 / step_hours,
        speed_model
      )
      keep <- dist_keep == 0 | speed_prob >= movement_eps

      # Convert the fitted one-dimensional speed density to an isotropic
      # probability per grid cell. Dividing by distance removes the increasing
      # number of cells in wider rings, while cos(latitude) accounts for cell
      # area on the regular longitude-latitude grid. Zero distance is handled
      # separately by the stay component of the transition model.
      cell_area <- pmax(cos(lat_to[idx[, 1]]), 0)
      spatial_prob <- numeric(length(dist_keep))
      moving <- dist_keep > 0
      spatial_prob[moving] <- speed_prob[moving] *
        cell_area[moving] /
        (dist_keep[moving] / 1000)

      drow[[i]] <- dr[keep]
      dcol[[i]] <- dc[keep]
      kernels[[i]] <- log(spatial_prob[keep])
    }
  )[["elapsed"]]

  # Pack the exact row-specific great-circle weights into one dense lookup.
  # Gibbs updates can then index by source row and grid displacement without
  # recomputing distances or evaluating the speed density.
  max_dcol <- max(abs(unlist(dcol)))
  lookup_ndrow <- 2L * max_row + 1L
  lookup_ndcol <- 2L * max_dcol + 1L
  move_log_lookup <- rep(-Inf, nrow * lookup_ndrow * lookup_ndcol)

  lookup_seconds <- system.time(
    for (i in seq_len(nrow)) {
      lookup_idx <- i +
        (drow[[i]] + max_row) * nrow +
        (dcol[[i]] + max_dcol) * nrow * lookup_ndrow
      move_log_lookup[lookup_idx] <- kernels[[i]]
    }
  )[["elapsed"]]

  # The directed movement kernel is normalized at its source cell. Precompute
  # these constants for every cell that can occur in the fixed light support
  # so the future-neighbour Gibbs term can be evaluated in the correct
  # direction without repeated full-footprint scans.
  normalization_seconds <- system.time({
    move_log_norm <- sampling_path_movement_log_norm(
      support_idx = support_idx,
      struct_mask = struct_mask,
      drow = drow,
      dcol = dcol,
      kernels = kernels,
      cell_row = cell_row,
      cell_col = cell_col,
      nrow = nrow,
      ncol = ncol
    )
  })[["elapsed"]]

  movement_kernel <- list(
    move_log_lookup = move_log_lookup,
    move_log_norm = move_log_norm,
    max_drow = max_row,
    max_dcol = max_dcol,
    lookup_ndrow = lookup_ndrow
  )

  if (!quiet) {
    cli::cli_progress_done()
  }

  prep <- list(
    lk = lk,
    kt = list(
      move_probability = move_probability,
      movement_kernel = movement_kernel,
      nrow = nrow,
      ncol = ncol,
      struct_mask = struct_mask,
      cell_row = cell_row,
      cell_col = cell_col,
      lat = g$lat,
      lon = g$lon,
      lat_rad = g$lat * pi / 180,
      lon_rad = g$lon * pi / 180,
      cos_lat = cos(g$lat * pi / 180)
    )
  )
  if (profile) {
    prep$profile <- list(
      grid_cells = ncell,
      support_cells = length(support_idx),
      kernel_timing = tibble::tibble(
        support_cells = length(support_idx),
        footprint_seconds = footprint_seconds,
        lookup_seconds = lookup_seconds,
        normalization_seconds = normalization_seconds,
        lookup_cells = length(move_log_lookup)
      )
    )
  }

  structure(prep, class = "sampling_path_prep")
}

#' Validate the Movement Model
#'
#' @return `TRUE` invisibly.
#' @noRd
sampling_path_movement_validate <- function(movement) {
  if ("seasons" %in% names(movement)) {
    cli::cli_abort(
      "`movement$seasons` is not supported; provide one movement model for the complete track."
    )
  }
  if (is.null(movement$move_stay) || !is.function(movement$move_stay)) {
    cli::cli_abort("The movement model must define a `move_stay` function.")
  }

  invisible(TRUE)
}

#' Normalize the Movement Kernel Over the Structural Domain
#'
#' @description
#' Compute source-cell normalization constants by convolving each destination
#' mask row with the row-specific movement weights. This keeps the outer loop at
#' the number of latitude rows instead of the number of supported grid cells.
#'
#' @return Numeric log-normalization vector with values for `support_idx`.
#' @noRd
sampling_path_movement_log_norm <- function(
  support_idx,
  struct_mask,
  drow,
  dcol,
  kernels,
  cell_row,
  cell_col,
  nrow,
  ncol
) {
  struct_matrix <- matrix(struct_mask, nrow = nrow, ncol = ncol)
  move_log_norm <- rep(NA_real_, nrow * ncol)
  support_by_row <- split(support_idx, cell_row[support_idx])

  for (source_row_name in names(support_by_row)) {
    source_row <- as.integer(source_row_name)
    norm_by_col <- numeric(ncol)

    for (dr in unique(drow[[source_row]])) {
      destination_row <- source_row + dr
      if (destination_row < 1L || destination_row > nrow) {
        next
      }

      keep <- drow[[source_row]] == dr
      dc <- dcol[[source_row]][keep]
      weight <- exp(kernels[[source_row]][keep])
      moving <- dr != 0L | dc != 0L
      dc <- dc[moving]
      weight <- weight[moving]
      if (!length(dc)) {
        next
      }

      max_dc <- max(abs(dc))
      filter_weight <- numeric(2L * max_dc + 1L)
      filter_weight[max_dc - dc + 1L] <- weight
      padded_mask <- c(
        numeric(max_dc),
        struct_matrix[destination_row, ],
        numeric(max_dc)
      )
      filtered <- stats::filter(padded_mask, filter_weight, sides = 2)
      norm_by_col <- norm_by_col +
        as.numeric(
          filtered[max_dc + seq_len(ncol)]
        )
    }

    source_idx <- support_by_row[[source_row_name]]
    norm <- norm_by_col[cell_col[source_idx]]
    move_log_norm[source_idx] <- ifelse(norm > 0, log(norm), -Inf)
  }

  move_log_norm
}
