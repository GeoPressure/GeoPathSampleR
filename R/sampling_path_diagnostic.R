# Check whether saved paths approach the sampler's computational support.
sampling_path_diagnostic_support <- function(paths, tag) {
  params <- attr(paths, "sampling_parameters")
  boundary_fraction <- 0.9
  map_boundary_mass_threshold <- 0.01
  likelihood <- GeoPressureR::tag2likelihood(
    tag,
    params$likelihood %||% "map_light"
  )
  map <- GeoPressureR::tag_prepare_likelihood(
    tag,
    likelihood = likelihood,
    thr_likelihood = params$thr_likelihood,
    thr_gs = Inf,
    quiet = TRUE
  )
  g <- GeoPressureR::map_expand(map$extent, map$scale)
  stap_id <- map$stap$stap_id[map$stap$include]
  stap_row <- match(stap_id, tag$stap$stap_id)
  unknown_stap <- is.na(tag$stap$known_lat[stap_row]) |
    is.na(tag$stap$known_lon[stap_row])

  map_boundary <- matrix(FALSE, nrow = g$dim[1], ncol = g$dim[2])
  map_boundary[c(1L, g$dim[1]), ] <- TRUE
  map_boundary[, c(1L, g$dim[2])] <- TRUE
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

  n_outside_likelihood <- sum(unlist(Map(
    function(id, mask) sum(!paths$ind[paths$stap_id == id] %in% which(mask)),
    stap_id,
    map$lk_mask
  )))

  grid_lat <- g$lat[(paths$ind - 1L) %% g$dim[1] + 1L]
  grid_lon <- g$lon[(paths$ind - 1L) %/% g$dim[1] + 1L]
  sampled_row <- (paths$ind - 1L) %% g$dim[1] + 1L
  sampled_col <- (paths$ind - 1L) %/% g$dim[1] + 1L
  sampled_unknown <- paths$stap_id %in% stap_id[unknown_stap]
  sampled_at_map_boundary <- sampled_unknown &
    map_boundary[cbind(sampled_row, sampled_col)]
  n <- nrow(paths)
  adjacent <- paths$chain[-n] == paths$chain[-1L] & paths$j[-n] == paths$j[-1L]
  distance_km <- if (any(adjacent)) {
    GeoPressureR::haversine_distance(
      cbind(grid_lon[-n][adjacent], grid_lat[-n][adjacent]),
      cbind(grid_lon[-1L][adjacent], grid_lat[-1L][adjacent])
    )
  } else {
    numeric()
  }
  move_distance_km <- distance_km[distance_km > 0]
  max_radius_km <- params$thr_gs * 24
  near_radius <- move_distance_km >= boundary_fraction * max_radius_km

  list(
    map_extent = map$extent,
    thr_likelihood = params$thr_likelihood,
    max_radius_km = max_radius_km,
    boundary_fraction = boundary_fraction,
    n_locations = n,
    n_outside_likelihood = n_outside_likelihood,
    share_outside_likelihood = n_outside_likelihood / n,
    map_boundary_mass_threshold = map_boundary_mass_threshold,
    n_staps_boundary_mode = sum(boundary_mode[unknown_stap]),
    n_staps_boundary_mass = sum(
      boundary_mass[unknown_stap] > map_boundary_mass_threshold
    ),
    max_boundary_mass = max(boundary_mass[unknown_stap]),
    n_locations_at_map_boundary = sum(sampled_at_map_boundary),
    share_locations_at_map_boundary = mean(sampled_at_map_boundary),
    n_moves = length(move_distance_km),
    n_near_radius = sum(near_radius),
    share_near_radius = if (length(move_distance_km)) mean(near_radius) else 0,
    max_move_distance_km = if (length(move_distance_km)) {
      max(move_distance_km)
    } else {
      0
    }
  )
}

#' Diagnostics for sampled paths
#'
#' @description
#' Compute standard MCMC, spatial-agreement, iteration, and computational-
#' support diagnostics for paths produced by [sampling_path()]. Known-location
#' staps are excluded, while daily and long staps are reported separately. Tags
#' without a `stap0` column are reported as daily only.
#' Values of R-hat above 1.05 need attention, and the effective-sample-size
#' target is approximately 100 draws per chain. Inspect computational-support
#' warnings before deciding that additional iterations are the appropriate
#' response. See the [interactive example diagnostic report](https://geopressure.org/GeoPathSampleR/sampling-path-diagnostic.html).
#'
#' @param paths data.frame returned by [sampling_path()].
#' @param tag the tag used by [sampling_path()].
#' @param report logical; show an interactive diagnostic report in the browser.
#'   The default is `TRUE` in an interactive R session and `FALSE` in scripts.
#' @param quiet logical to hide printed diagnostics.
#'
#' @return Invisibly, a list containing all computed diagnostics.
#' @family sampling_path
#' @export
sampling_path_diagnostic <- function(
  paths,
  tag,
  report = interactive(),
  quiet = FALSE
) {
  time_start <- proc.time()[["elapsed"]]
  if (!quiet) {
    cli::cli_alert_info("Computing sampling diagnostics")
  }
  stap <- tag$stap
  sampler_params <- attr(paths, "sampling_parameters") %||% list()
  out <- sampling_path_diagnostic_compute(paths, stap)

  saved_j <- sort(unique(paths$j))
  checkpoint_share <- c(0.25, 0.5, 0.75, 1)
  checkpoint_n <- unique(pmax(1L, floor(length(saved_j) * checkpoint_share)))
  checkpoint_label <- paste0(round(checkpoint_n / length(saved_j) * 100), "%")
  checkpoint_label[checkpoint_n == length(saved_j)] <- "Full"
  out$iteration <- do.call(
    rbind,
    lapply(seq_along(checkpoint_n), function(i) {
      checkpoint <- sampling_path_diagnostic_compute(
        paths[paths$j <= saved_j[checkpoint_n[i]], , drop = FALSE],
        stap,
        spatial = FALSE,
        sample_summary = out$sample_summary[
          out$sample_summary$j <= saved_j[checkpoint_n[i]],
          ,
          drop = FALSE
        ]
      )
      sampling_path_diagnostic_checkpoint(
        checkpoint,
        checkpoint_label[i],
        checkpoint_n[i]
      )
    })
  )

  out$support <- if (
    all(c("thr_likelihood", "thr_gs") %in% names(sampler_params))
  ) {
    sampling_path_diagnostic_support(paths, tag)
  } else {
    NULL
  }
  out$n_chain <- length(unique(paths$chain))
  out$n_iter <- length(saved_j)
  out$n_known_stap <- nrow(stap) - nrow(out$marginal)
  out$sampler_params <- sampler_params
  out$tag_id <- tag$param$id

  if (!quiet) {
    sampling_path_diagnostic_print(out)
  }
  if (report) {
    if (!quiet) {
      cli::cli_alert_info("Rendering diagnostic report")
    }
    sampling_path_diagnostic_render(out, paths, tag)
    if (!quiet) cli::cli_alert_success("Diagnostic report ready")
  } else if (!quiet) {
    cli::cli_alert_success(
      "Sampling diagnostics computed in {round(proc.time()[['elapsed']] - time_start, 1)} s"
    )
  }

  invisible(out)
}

sampling_path_diagnostic_compute <- function(
  paths,
  stap,
  path_summary = TRUE,
  spatial = TRUE,
  sample_summary = NULL
) {
  assertthat::assert_that(is.data.frame(paths))
  assertthat::assert_that(all(
    c("j", "chain", "stap_id", "ind", "lat", "lon") %in% names(paths)
  ))
  assertthat::assert_that(all(
    c(
      "stap_id",
      "start",
      "end",
      "known_lat",
      "known_lon"
    ) %in%
      names(stap)
  ))

  stap0 <- stap$stap0 %||% rep(FALSE, nrow(stap))
  sampled_stap <- stap[
    !(is.finite(stap$known_lat) & is.finite(stap$known_lon)),
    "stap_id",
    drop = FALSE
  ]
  sampled_stap$stap0 <- stap0[match(sampled_stap$stap_id, stap$stap_id)]
  sampled_stap$stap_type <- ifelse(sampled_stap$stap0, "long", "daily")
  sampled_paths <- paths[
    paths$stap_id %in% sampled_stap$stap_id,
    ,
    drop = FALSE
  ]

  marginal <- do.call(
    rbind,
    lapply(split(sampled_paths, sampled_paths$stap_id), function(df) {
      out <- lapply(c("lat", "lon"), function(v) {
        mat <- sampling_path_iter_chain_matrix(df, v)
        sampling_path_posterior_diagnostic(mat)
      })
      data.frame(
        stap_id = df$stap_id[1],
        stap_type = sampled_stap$stap_type[match(
          df$stap_id[1],
          sampled_stap$stap_id
        )],
        lat_rhat = out[[1]][["rhat"]],
        lat_ess_bulk = out[[1]][["ess_bulk"]],
        lat_ess_tail = out[[1]][["ess_tail"]],
        lat_mcse_mean = out[[1]][["mcse_mean"]],
        lon_rhat = out[[2]][["rhat"]],
        lon_ess_bulk = out[[2]][["ess_bulk"]],
        lon_ess_tail = out[[2]][["ess_tail"]],
        lon_mcse_mean = out[[2]][["mcse_mean"]],
        row.names = NULL,
        check.names = FALSE
      )
    })
  )

  if (path_summary) {
    if (is.null(sample_summary)) {
      sample_summary <- sampling_path_diagnostic_sample_summary(paths)
    }

    summary_diagnostics <- do.call(
      rbind,
      lapply(c("n_stays", "n_unique_cells", "path_length_km"), function(var) {
        mat <- sampling_path_iter_chain_matrix(sample_summary, var)
        diagnostic <- sampling_path_posterior_diagnostic(mat)
        data.frame(
          metric = var,
          rhat = diagnostic[["rhat"]],
          ess_bulk = diagnostic[["ess_bulk"]],
          ess_tail = diagnostic[["ess_tail"]],
          mcse_mean = diagnostic[["mcse_mean"]],
          row.names = NULL,
          check.names = FALSE
        )
      })
    )
  } else {
    sample_summary <- NULL
    summary_diagnostics <- NULL
  }

  per_stap <- if (!spatial) {
    NULL
  } else {
    do.call(
      rbind,
      lapply(split(sampled_paths, sampled_paths$stap_id), function(df) {
        chain_tab <- split(df$ind, df$chain)
        all_ind <- sort(unique(df$ind))
        prob_mat <- do.call(
          cbind,
          lapply(chain_tab, function(ind) {
            as.numeric(table(factor(ind, levels = all_ind))) / length(ind)
          })
        )
        if (is.null(dim(prob_mat))) {
          prob_mat <- matrix(prob_mat, ncol = 1L)
        }
        pair_overlap <- if (ncol(prob_mat) < 2L) {
          NA_real_
        } else {
          pairs <- utils::combn(seq_len(ncol(prob_mat)), 2L)
          mean(apply(pairs, 2L, function(idx) {
            sum(pmin(prob_mat[, idx[1]], prob_mat[, idx[2]]))
          }))
        }
        chain_centroid <- stats::aggregate(
          cbind(lat, lon) ~ chain,
          data = df,
          FUN = mean
        )
        chain_separation_km <- if (nrow(chain_centroid) < 2L) {
          NA_real_
        } else {
          chain_pairs <- utils::combn(seq_len(nrow(chain_centroid)), 2L)
          GeoPressureR::haversine_distance(
            cbind(
              chain_centroid$lon[chain_pairs[1, ]],
              chain_centroid$lat[chain_pairs[1, ]]
            ),
            cbind(
              chain_centroid$lon[chain_pairs[2, ]],
              chain_centroid$lat[chain_pairs[2, ]]
            )
          )
        }
        overall_prob <- rowMeans(prob_mat)
        data.frame(
          stap_id = df$stap_id[1],
          stap_type = sampled_stap$stap_type[match(
            df$stap_id[1],
            sampled_stap$stap_id
          )],
          n_unique_cells = length(all_ind),
          modal_cell_share = max(overall_prob),
          occupancy_entropy = -sum(
            overall_prob[overall_prob > 0] * log(overall_prob[overall_prob > 0])
          ),
          mean_pairwise_overlap = pair_overlap,
          n_distinct_chain_modes = length(unique(apply(
            prob_mat,
            2L,
            function(p) all_ind[which.max(p)]
          ))),
          median_chain_centroid_separation_km = stats::median(
            chain_separation_km
          ),
          maximum_chain_centroid_separation_km = max(chain_separation_km),
          row.names = NULL,
          check.names = FALSE
        )
      })
    )
  }

  list(
    marginal = marginal,
    sample_summary = sample_summary,
    summary_diagnostics = summary_diagnostics,
    per_stap = per_stap,
    path_mode_share = NA_real_
  )
}

sampling_path_diagnostic_sample_summary <- function(paths) {
  paths <- paths[order(paths$chain, paths$j, paths$stap_id), , drop = FALSE]
  sample_start <- c(
    TRUE,
    paths$chain[-1L] != paths$chain[-nrow(paths)] |
      paths$j[-1L] != paths$j[-nrow(paths)]
  )
  sample_i <- cumsum(sample_start)
  n_sample <- max(sample_i)
  stay_start <- sample_start
  stay_start[-1L] <- sample_start[-1L] |
    paths$ind[-1L] != paths$ind[-nrow(paths)]

  n_stays <- tabulate(sample_i[stay_start], nbins = n_sample)
  unique_cell <- !duplicated(cbind(sample_i, paths$ind))
  n_unique_cells <- tabulate(sample_i[unique_cell], nbins = n_sample)

  transition_i <- which(stay_start & !sample_start)
  path_length_km <- numeric(n_sample)
  if (length(transition_i)) {
    transition_distance <- GeoPressureR::haversine_distance(
      cbind(paths$lon[transition_i - 1L], paths$lat[transition_i - 1L]),
      cbind(paths$lon[transition_i], paths$lat[transition_i])
    )
    transition_total <- rowsum(
      transition_distance,
      sample_i[transition_i],
      reorder = FALSE
    )
    path_length_km[as.integer(rownames(transition_total))] <-
      transition_total[, 1]
  }

  sample_first <- which(sample_start)
  data.frame(
    j = paths$j[sample_first],
    chain = paths$chain[sample_first],
    n_stays = n_stays,
    n_unique_cells = n_unique_cells,
    path_length_km = path_length_km,
    row.names = NULL,
    check.names = FALSE
  )
}

sampling_path_diagnostic_checkpoint <- function(
  x,
  checkpoint,
  saved_per_chain
) {
  finite_max <- function(value) {
    if (any(is.finite(value))) max(value, na.rm = TRUE) else NA_real_
  }
  daily <- x$marginal[x$marginal$stap_type == "daily", , drop = FALSE]
  long <- x$marginal[x$marginal$stap_type == "long", , drop = FALSE]
  daily_rhat <- c(daily$lat_rhat, daily$lon_rhat)
  long_rhat <- c(long$lat_rhat, long$lon_rhat)

  data.frame(
    checkpoint = checkpoint,
    saved_per_chain = saved_per_chain,
    daily_max_rhat = finite_max(daily_rhat),
    daily_rhat_over_1_05 = sum(daily_rhat > 1.05, na.rm = TRUE),
    daily_median_bulk_ess = stats::median(
      c(daily$lat_ess_bulk, daily$lon_ess_bulk),
      na.rm = TRUE
    ),
    long_max_rhat = finite_max(long_rhat),
    path_max_rhat = finite_max(x$summary_diagnostics$rhat),
    row.names = NULL,
    check.names = FALSE
  )
}

sampling_path_diagnostic_print <- function(x) {
  rhat <- c(
    x$marginal$lat_rhat,
    x$marginal$lon_rhat,
    x$summary_diagnostics$rhat
  )
  bulk_ess <- c(x$marginal$lat_ess_bulk, x$marginal$lon_ess_bulk)
  tail_ess <- c(x$marginal$lat_ess_tail, x$marginal$lon_ess_tail)
  ess_target <- 100 * x$n_chain
  support <- if (is.null(x$support)) {
    "not assessed (sampler settings unavailable)"
  } else if (
    x$support$n_outside_likelihood > 0L ||
      x$support$n_near_radius > 0L ||
      x$support$n_staps_boundary_mode > 0L ||
      x$support$n_staps_boundary_mass > 0L
  ) {
    "needs attention"
  } else {
    "no likelihood, map-extent, or movement-radius warning"
  }
  cli::cli_h2("Sampling diagnostics")
  cli::cli_bullets(c(
    "Chains: {x$n_chain}; saved iterations per chain: {x$n_iter}",
    "Sampled staps: {nrow(x$marginal)} ({sum(x$marginal$stap_type == 'daily')} daily, {sum(x$marginal$stap_type == 'long')} long)",
    "Maximum R-hat: {signif(max(rhat, na.rm = TRUE), 3)}; values above 1.05 need attention",
    "Minimum bulk ESS: {round(min(bulk_ess, na.rm = TRUE))}; minimum tail ESS: {round(min(tail_ess, na.rm = TRUE))}; target: {ess_target}",
    "Computational support: {support}"
  ))
}

sampling_path_iter_chain_matrix <- function(df, value_col) {
  df <- df[!is.na(df[[value_col]]), c("j", "chain", value_col)]
  df <- df[order(df$j, df$chain), , drop = FALSE]
  out <- matrix(
    NA_real_,
    nrow = length(unique(df$j)),
    ncol = length(unique(df$chain))
  )
  out[cbind(
    match(df$j, sort(unique(df$j))),
    match(df$chain, sort(unique(df$chain)))
  )] <- df[[value_col]]
  out
}

sampling_path_posterior_diagnostic <- function(x) {
  c(
    rhat = posterior::rhat(x),
    ess_bulk = posterior::ess_bulk(x),
    ess_tail = posterior::ess_tail(x),
    mcse_mean = posterior::mcse_mean(x)
  )
}

sampling_path_diagnostic_render <- function(
  diagnostic,
  paths,
  tag,
  browse = TRUE,
  output_file = NULL
) {
  if (is.null(output_file)) {
    report_dir <- tempfile("sampling_path_diagnostic_")
    dir.create(report_dir)
    report_file <- file.path(report_dir, "sampling_path_diagnostic.html")
    on.exit(unlink(report_file), add = TRUE)
  } else {
    report_file <- normalizePath(output_file, mustWork = FALSE)
    dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)
  }
  diagnostic_file <- tempfile(fileext = ".rds")
  params_file <- tempfile(fileext = ".json")
  saveRDS(
    list(
      diagnostic = diagnostic,
      paths = paths,
      stap = tag$stap,
      map_extent = tag$param$tag_set_map$extent %||% NULL
    ),
    diagnostic_file
  )
  on.exit(unlink(diagnostic_file), add = TRUE)
  on.exit(unlink(params_file), add = TRUE)
  jsonlite::write_json(
    list(diagnostic_file = diagnostic_file),
    params_file,
    auto_unbox = TRUE
  )

  template <- system.file(
    "report",
    "sampling_path_diagnostic.qmd",
    package = "GeoPathSampleR"
  )
  if (!nzchar(template)) {
    cli::cli_abort(
      "The GeoPathSampleR diagnostic-report template is unavailable."
    )
  }
  render_name <- paste0("sampling_path_diagnostic_", Sys.getpid(), ".html")
  rendered_file <- file.path(dirname(template), render_name)
  resource_dir <- file.path(
    dirname(template),
    paste0(tools::file_path_sans_ext(render_name), "_files")
  )
  on.exit(unlink(rendered_file), add = TRUE)
  on.exit(unlink(resource_dir, recursive = TRUE), add = TRUE)
  render_output <- withr::with_dir(
    dirname(template),
    system2(
      "quarto",
      c(
        "render",
        shQuote(basename(template)),
        "--to",
        "html",
        "--output",
        shQuote(render_name),
        "--execute-params",
        shQuote(params_file)
      ),
      stdout = TRUE,
      stderr = TRUE
    )
  )
  status <- attr(render_output, "status") %||% 0L
  if (status != 0L) {
    cli::cli_abort(c(
      "Failed to render the sampling-path diagnostic report.",
      "x" = paste(utils::tail(render_output, 8), collapse = "\n")
    ))
  }
  if (!file.copy(rendered_file, report_file, overwrite = TRUE)) {
    cli::cli_abort(
      "Failed to open the rendered sampling-path diagnostic report."
    )
  }
  if (browse) {
    utils::browseURL(report_file)
  }

  invisible(report_file)
}
