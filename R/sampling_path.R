#' Gibbs sampling for stationary paths
#'
#' @description
#' Prepare the likelihood and movement inputs, run the Gibbs sampler for one or
#' more chains, and return the sampled path locations in long format.
#'
#' @param tag a GeoPressureR `tag` object containing the likelihood maps and
#'   `tag_set_map` parameters.
#' @param movement one movement model applied across the complete track. The
#'   parameter list must contain the fields passed to `speed2prob()` and a
#'   `move_stay` function returning the probability of moving after a stay
#'   duration. Users can optionally provide `approx_eps` to drop extremely small
#'   movement weights from the precomputed neighbourhood (default `1e-6`).
#' @param iter number of Gibbs iterations.
#' @param chains number of chains to run.
#' @param warmup number of iterations to discard at the beginning.
#' @param thin thinning interval for saved samples.
#' @param refresh progress message update frequency, in iterations.
#' @param block_interval interval between block iterations. All remaining
#'   iterations are site iterations that update each eligible period once.
#'   Use `0` or `Inf` to disable block iterations. Values greater than one are
#'   required so site iterations remain part of the chain.
#' @param long_period_light_only logical. `TRUE` (the default) draws unknown
#'   `stap0` locations directly from their retained light likelihood and samples
#'   the daily path conditionally without allowing movement to update them.
#'   `FALSE` updates those locations with both light and movement information.
#'   Tags without a `stap0` column have no long periods.
#' @param component_weights named non-negative weights for the light likelihood,
#'   daily movement model, and interval route model. The default gives all three
#'   components their calibrated weight. A zero removes that component's soft
#'   contribution while retaining structural support constraints.
#' @param route_detour non-negative multiplier of the empirically predicted
#'   excess route distance. `0` targets direct routes, `1` retains the empirical
#'   prediction, and values above one favour more detoured routes.
#' @param route_model calibrated route-model parameter list. The default uses
#'   the maintained duration-and-direct-distance calibration. Supply cross-fitted
#'   parameters only for validation workflows.
#' @param likelihood field of the `tag` list containing the likelihood map.
#' @param thr_likelihood threshold of percentile to keep likely locations in
#'   each stationary period. Lower values keep fewer nodes and speed up the
#'   sampler, but exclude more of the light likelihood surface entirely. The
#'   default `0.99` is a conservative state-space reduction; lower values such
#'   as `0.95` are faster but more aggressive.
#' @param thr_gs maximum ground speed, in km/h, used as the hard radius of the
#'   movement kernel. Larger values allow broader transitions and can improve
#'   mixing when movements are uncertain, but they also increase computation.
#'   Smaller values are faster but enforce a tighter movement radius.
#' @param workers number of parallel workers used for independent chains.
#'   Use 1 for sequential execution.
#' @param seed optional RNG seed for reproducible chains.
#' @param quiet logical to hide messages.
#'
#' @return A data.frame with columns `j`, `chain`, `stap_id`, `ind`, `lat`, and
#'   `lon`. The result carries `type = "sampling"` and the scalar sampler
#'   settings in a `sampling_parameters` attribute for reproducible diagnostics.
#'
#' @details
#' `sampling_path()` separates structural constraints from computational
#' approximations. Impossible geography, such as water after masking, is removed
#' from the state space. The arguments `thr_likelihood`, `thr_gs`,
#' and `movement$approx_eps` control how aggressively the sampler reduces the
#' fixed light support and movement kernel for speed.
#'
#' A practical tuning strategy is:
#' - start with the defaults for final inference;
#' - use a smaller `thr_likelihood` only for explicitly approximate exploratory
#'   runs;
#' - if no movement-feasible path exists, increase `thr_likelihood` or
#'   `thr_gs`.
#'
#' @examples
#' example_dir <- system.file("extdata", package = "GeoPathSampleR")
#' tag <- GeoPressureR::tag_create(
#'   "14OI",
#'   crop_start = "2015-07-17",
#'   crop_end = "2016-07-11",
#'   directory = file.path(example_dir, "data/raw-tag/14OI"),
#'   assert_pressure = FALSE,
#'   quiet = TRUE
#' )
#' tag <- GeoPressureR::twilight_create(tag)
#' tag <- GeoPressureR::twilight_label_read(
#'   tag,
#'   file = file.path(example_dir, "data/twilight-label/14OI-labeled.csv")
#' )
#' tag <- GeoPressureR::tag_stap_daily(
#'   tag,
#'   stap_long = file.path(example_dir, "data/stap-label/14OI.csv"),
#'   movement_period = "day",
#'   quiet = TRUE
#' )
#' tag <- GeoPressureR::tag_set_map(
#'   tag,
#'   extent = c(-20, 40, -40, 60),
#'   scale = 1
#' )
#' tag <- GeoPressureR::geolight_map(
#'   tag,
#'   twl_calib_adjust = 1,
#'   fitted_location_duration = 30,
#'   twl_llp = \(n) 1.75 * log(n) / n,
#'   quiet = TRUE
#' )
#'
#' paths <- sampling_path(
#'   tag,
#'   iter = 3,
#'   warmup = 1,
#'   seed = 1,
#'   quiet = TRUE
#' )
#' @family sampling_path
#' @noRd
sampling_path_with_route_model <- function(
  tag,
  iter,
  likelihood = "map_light",
  movement = list(
    method = "gamma",
    shape = 1.412246,
    scale = 8.909248,
    low_speed_fix = 0.001,
    zero_speed_ratio = 0,
    move_stay_parameters = list(
      p_0 = 0.647893,
      p_inf = 0.1318134,
      tau = 0.8624705
    ),
    move_stay = function(t) {
      0.1318134 + (0.647893 - 0.1318134) * exp(-t / 0.8624705)
    }
  ),
  component_weights = c(light = 1, movement = 1, route = 1),
  route_detour = 1,
  long_period_light_only = TRUE,
  chains = 1,
  warmup = floor(iter / 4),
  thin = 1,
  block_interval = 2,
  thr_likelihood = 0.99,
  thr_gs = 2000 / 24,
  refresh = 10,
  workers = 1,
  seed = NULL,
  quiet = FALSE,
  route_model = sampling_path_route_model()
) {
  time_start <- proc.time()[["elapsed"]]

  assertthat::assert_that(
    is.numeric(chains),
    length(chains) == 1,
    chains > 0,
    chains == floor(chains),
    is.numeric(workers),
    length(workers) == 1,
    workers > 0,
    workers == floor(workers),
    is.numeric(iter),
    length(iter) == 1,
    iter > 0,
    iter == floor(iter),
    is.numeric(warmup),
    length(warmup) == 1,
    warmup >= 0,
    warmup == floor(warmup),
    is.numeric(thin),
    length(thin) == 1,
    thin > 0,
    thin == floor(thin),
    is.numeric(refresh),
    length(refresh) == 1,
    refresh > 0,
    refresh == floor(refresh),
    is.numeric(block_interval),
    length(block_interval) == 1,
    block_interval >= 0,
    is.infinite(block_interval) ||
      block_interval == 0 ||
      block_interval >= 2,
    is.infinite(block_interval) ||
      block_interval == floor(block_interval),
    iter > warmup,
    is.null(seed) ||
      (is.numeric(seed) && length(seed) == 1 && seed == floor(seed)),
    is.logical(long_period_light_only),
    length(long_period_light_only) == 1,
    !is.na(long_period_light_only),
    is.numeric(component_weights),
    identical(sort(names(component_weights)), c("light", "movement", "route")),
    all(is.finite(component_weights)),
    all(component_weights >= 0),
    is.numeric(route_detour),
    length(route_detour) == 1,
    is.finite(route_detour),
    route_detour >= 0
  )
  component_weights <- component_weights[c("light", "movement", "route")]

  chains <- as.integer(chains)
  if (is.infinite(block_interval)) {
    block_interval <- 0L
  } else {
    block_interval <- as.integer(block_interval)
  }
  workers <- min(as.integer(workers), chains)

  prep <- sampling_path_prepare(
    tag = tag,
    likelihood = likelihood,
    movement = movement,
    thr_likelihood = thr_likelihood,
    thr_gs = thr_gs,
    quiet = quiet
  )

  lk <- prep$lk
  kt <- prep$kt
  stap_ids <- vapply(lk, `[[`, numeric(1), "stap_id")
  nstap <- length(lk)
  stap_row <- match(stap_ids, tag$stap$stap_id)
  stap0 <- tag$stap$stap0 %||% rep(FALSE, nrow(tag$stap))
  known_stap <- if (all(c("known_lat", "known_lon") %in% names(tag$stap))) {
    is.finite(tag$stap$known_lat[stap_row]) &
      is.finite(tag$stap$known_lon[stap_row])
  } else {
    rep(FALSE, nstap)
  }
  long_period_i <- if (long_period_light_only) {
    which(stap0[stap_row] & !known_stap)
  } else {
    integer()
  }
  prepared_route_prior <- sampling_path_prepare_route_prior(
    tag = tag,
    stap_ids = stap_ids,
    route_prior = if (component_weights[["route"]] > 0) route_model else NULL,
    weight = component_weights[["route"]],
    detour = route_detour
  )
  n_save <- floor((iter - warmup) / thin)

  if (n_save <= 0L) {
    cli::cli_abort(
      "No samples to save with the current `iter`, `warmup` and `thin`."
    )
  }

  chain_seeds <- if (is.null(seed)) {
    NULL
  } else {
    withr::with_seed(seed, sample.int(.Machine$integer.max, chains))
  }

  # Keep progress handling simple: one bar for sequential runs, one message for
  # parallel runs.
  use_parallel <- workers > 1L && chains > 1L && .Platform$OS.type != "windows"
  progress_bar <- NULL
  progress_closed <- FALSE
  progress_pos <- 0L
  progress_fn <- NULL

  if (!quiet) {
    cli::cli_alert_info(
      if (use_parallel) {
        "Starting sampling: {chains} chain{?s}, {iter} iteration{?s}, warmup {warmup}, thin {thin}, {workers} worker{?s}."
      } else {
        "Starting sampling: {chains} chain{?s}, {iter} iteration{?s}, warmup {warmup}, thin {thin}."
      }
    )
  }

  if (!quiet && !use_parallel) {
    progress_bar <- utils::txtProgressBar(
      min = 0,
      max = chains * iter,
      style = 3
    )
    on.exit(
      if (!progress_closed) {
        close(progress_bar)
      },
      add = TRUE
    )
    progress_fn <- function(amount) {
      progress_pos <<- progress_pos + amount
      utils::setTxtProgressBar(progress_bar, progress_pos)
    }
  }

  chain_ids <- seq_len(chains)

  # Run one chain at a time, reusing the same chain runner in sequential and
  # parallel execution.
  chain_samples <- if (use_parallel) {
    parallel::mclapply(
      chain_ids,
      sampling_path_run_chain,
      lk = lk,
      kt = kt,
      iter = iter,
      warmup = warmup,
      thin = thin,
      block_interval = block_interval,
      long_period_i = long_period_i,
      component_weights = component_weights,
      route_prior = prepared_route_prior,
      refresh = refresh,
      quiet = quiet,
      progress_fn = NULL,
      seed_i = chain_seeds,
      mc.cores = workers
    )
  } else {
    lapply(
      chain_ids,
      sampling_path_run_chain,
      lk = lk,
      kt = kt,
      iter = iter,
      warmup = warmup,
      thin = thin,
      block_interval = block_interval,
      long_period_i = long_period_i,
      component_weights = component_weights,
      route_prior = prepared_route_prior,
      refresh = refresh,
      quiet = quiet,
      progress_fn = progress_fn,
      seed_i = chain_seeds
    )
  }

  invalid_chain <- which(!vapply(chain_samples, is.matrix, logical(1)))

  if (length(invalid_chain)) {
    bad_i <- invalid_chain[1]
    bad_value <- chain_samples[[bad_i]]

    if (inherits(bad_value, "try-error")) {
      cli::cli_abort(c(
        "Parallel sampling failed in chain {.val {bad_i}}.",
        "x" = as.character(bad_value)
      ))
    }

    cli::cli_abort(
      "Sampling chain {.val {bad_i}} returned an invalid result of class {.val {class(bad_value)[1]}}."
    )
  }

  repair_diagnostics <- lapply(
    chain_samples,
    attr,
    which = "repair_diagnostics"
  )
  repair_diagnostics <- repair_diagnostics[
    !vapply(repair_diagnostics, is.null, logical(1))
  ]
  n_repair_iterations <- sum(vapply(
    repair_diagnostics,
    `[[`,
    integer(1),
    "n_repair_iterations"
  ))
  n_repaired_sites <- sum(vapply(
    repair_diagnostics,
    `[[`,
    integer(1),
    "n_repaired_sites"
  ))
  max_repaired_sites <- if (length(repair_diagnostics)) {
    max(vapply(
      repair_diagnostics,
      `[[`,
      integer(1),
      "max_repaired_sites"
    ))
  } else {
    0L
  }

  # Convert each sampled matrix back to the long path format used elsewhere.
  paths <- do.call(
    rbind,
    lapply(chain_ids, function(chain_id) {
      samples <- chain_samples[[chain_id]]
      ind <- as.integer(as.vector(t(samples)))

      data.frame(
        j = rep(seq_len(nrow(samples)), each = nstap),
        chain = rep.int(chain_id, nrow(samples) * nstap),
        stap_id = rep.int(stap_ids, times = nrow(samples)),
        ind = ind,
        lat = kt$lat[(ind - 1L) %% kt$nrow + 1L],
        lon = kt$lon[(ind - 1L) %/% kt$nrow + 1L],
        row.names = NULL,
        check.names = FALSE
      )
    })
  )

  # Rows are already ordered by chain, saved iteration, and period because
  # `chain_ids`, `samples`, and `stap_ids` are each constructed in that order.
  # Avoid sorting this potentially large table after every validation run.

  # Known locations remain fixed in the final output even if the sampler kept
  # them in the discrete path index representation. The sampled rows repeat
  # `stap_ids` in order, so direct replacement avoids a full-table join.
  if (all(c("known_lat", "known_lon") %in% names(tag$stap))) {
    known_lat <- tag$stap$known_lat[stap_row]
    known_lon <- tag$stap$known_lon[stap_row]
    known_idx <- !is.na(known_lat) & !is.na(known_lon)

    if (any(known_idx)) {
      known_lat <- rep(known_lat, length.out = nrow(paths))
      known_lon <- rep(known_lon, length.out = nrow(paths))
      known_idx <- !is.na(known_lat) & !is.na(known_lon)
      paths$lat[known_idx] <- known_lat[known_idx]
      paths$lon[known_idx] <- known_lon[known_idx]
    }
  }

  row.names(paths) <- NULL
  class(paths) <- c("sampling_path", "data.frame")
  attr(paths, "type") <- "sampling"
  attr(paths, "sampling_parameters") <- list(
    iter = iter,
    chains = chains,
    warmup = warmup,
    thin = thin,
    refresh = refresh,
    block_interval = block_interval,
    long_period_light_only = long_period_light_only,
    n_light_only_long_periods = length(long_period_i),
    n_repair_iterations = n_repair_iterations,
    n_repaired_sites = n_repaired_sites,
    max_repaired_sites = max_repaired_sites,
    component_weights = component_weights,
    route_detour = route_detour,
    route_model = route_model,
    likelihood = likelihood,
    thr_likelihood = thr_likelihood,
    thr_gs = thr_gs,
    movement_approx_eps = movement$approx_eps %||% 1e-6,
    movement_method = movement$method,
    movement_shape = movement$shape,
    movement_scale = movement$scale,
    workers = workers,
    seed = seed
  )

  if (!is.null(progress_bar)) {
    close(progress_bar)
    progress_closed <- TRUE
  }

  if (!quiet) {
    elapsed <- proc.time()[['elapsed']] - time_start
    hours <- floor(elapsed / 3600)
    minutes <- floor((elapsed %% 3600) / 60)
    seconds <- elapsed %% 60
    cli::cli_alert_success(
      "Sampling finished in {hours} h {minutes} min {round(seconds, 1)} s."
    )
  }

  paths
}

#' Gibbs sampling for stationary paths
#'
#' @description
#' Sample stationary trajectories from light likelihood maps and calibrated
#' movement components.
#'
#' @param tag GeoPressureR tag object containing likelihood maps.
#' @param iter number of Gibbs iterations.
#' @param likelihood tag field containing the light likelihood map.
#' @param movement movement model containing a speed kernel and stay/move
#'   probabilities.
#' @param component_weights named non-negative weights for the light,
#'   movement, and route components.
#' @param route_detour non-negative multiplier of the empirical excess route
#'   distance; `1` retains the calibrated default.
#' @param long_period_light_only whether unknown long periods are drawn from
#'   their light likelihood without movement updates.
#' @param chains number of independent chains.
#' @param warmup number of initial iterations to discard.
#' @param thin interval between saved samples.
#' @param block_interval interval between block updates; use `0` to disable
#'   them.
#' @param thr_likelihood retained light-likelihood percentile.
#' @param thr_gs maximum movement speed in km/h used for hard support.
#' @param refresh progress update interval in iterations.
#' @param workers number of parallel workers for independent chains.
#' @param seed optional random seed.
#' @param quiet whether to suppress progress messages.
#' @return A data.frame with columns `j`, `chain`, `stap_id`, `ind`, `lat`, and
#'   `lon`. The result carries `type = "sampling"` and sampler settings in a
#'   `sampling_parameters` attribute.
#' @examples
#' \dontrun{
#' paths <- sampling_path(tag, iter = 1000, chains = 4, seed = 1)
#' }
#' @family sampling_path
#' @export
sampling_path <- function(
  tag,
  iter,
  likelihood = "map_light",
  movement = list(
    method = "gamma",
    shape = 1.412246,
    scale = 8.909248,
    low_speed_fix = 0.001,
    zero_speed_ratio = 0,
    move_stay_parameters = list(
      p_0 = 0.647893,
      p_inf = 0.1318134,
      tau = 0.8624705
    ),
    move_stay = function(t) {
      0.1318134 + (0.647893 - 0.1318134) * exp(-t / 0.8624705)
    }
  ),
  component_weights = c(light = 1, movement = 1, route = 1),
  route_detour = 1,
  long_period_light_only = TRUE,
  chains = 1,
  warmup = floor(iter / 4),
  thin = 1,
  block_interval = 2,
  thr_likelihood = 0.99,
  thr_gs = 2000 / 24,
  refresh = 10,
  workers = 1,
  seed = NULL,
  quiet = FALSE
) {
  sampling_path_with_route_model(
    tag = tag,
    iter = iter,
    likelihood = likelihood,
    movement = movement,
    component_weights = component_weights,
    route_detour = route_detour,
    long_period_light_only = long_period_light_only,
    chains = chains,
    warmup = warmup,
    thin = thin,
    block_interval = block_interval,
    thr_likelihood = thr_likelihood,
    thr_gs = thr_gs,
    refresh = refresh,
    workers = workers,
    seed = seed,
    quiet = quiet,
    route_model = sampling_path_route_model()
  )
}

#' Return the Maintained Route-Length Prior
#'
#' @return A named list of fixed calibration parameters, or `NULL` when the
#'   route-length prior is disabled.
#' @noRd
sampling_path_route_model <- function() {
  list(
    intercept = 0.127076604313848,
    duration_coefficient = 0.120519227570285,
    distance_coefficient = -0.0836208733274686,
    duration_center = 3.36513833184645,
    distance_center = 8.13763879220316,
    positive_link_scale = 10,
    residual_sd = 0.115013822553951,
    min_direct_distance_km = 100
  )
}

#' Run One Sampling Chain
#'
#' @description
#' Internal worker used by [sampling_path()] to initialize one chain, run the
#' Gibbs iterations, and return the saved sampled indices.
#'
#' @param chain_id integer chain identifier.
#' @param lk list of prepared per-stap likelihood objects.
#' @param kt prepared transition lookup object.
#' @inheritParams sampling_path
#' @param progress_fn optional progress callback for sequential runs.
#' @param seed_i optional vector of chain seeds.
#'
#' @return Integer matrix of sampled path indices, one row per saved draw.
#' @noRd
sampling_path_run_chain <- function(
  chain_id,
  lk,
  kt,
  iter,
  warmup,
  thin,
  block_interval,
  long_period_i,
  component_weights,
  route_prior,
  refresh,
  quiet,
  progress_fn,
  seed_i = NULL
) {
  if (!is.null(seed_i)) {
    withr::local_seed(seed_i[chain_id])
  }

  nstap <- length(lk)
  n_save <- floor((iter - warmup) / thin)

  # Initialize the path so every adjacent stap is movement-feasible before the
  # stochastic iterations start.
  path_idx <- sampling_path_initialize(
    lk = lk,
    kt = kt,
    chain_id = chain_id,
    init_power = 4
  )
  route_distance <- sampling_path_route_distance_state(
    path_idx,
    kt,
    route_prior
  )
  route_edge_distance <- route_distance$edge
  route_interval_distance <- route_distance$interval

  samples <- matrix(NA_integer_, nrow = n_save, ncol = nstap)
  save_idx <- 0L
  site_iteration_i <- 0L
  repair_iteration_count <- 0L
  repaired_site_count <- 0L
  max_repaired_sites <- 0L
  progress_every <- max(1L, as.integer(refresh))
  last_progress <- 0L

  for (iter_i in seq_len(iter)) {
    if (length(long_period_i)) {
      refreshed <- sampling_path_refresh_long_periods(
        path_idx = path_idx,
        long_period_i = long_period_i,
        lk = lk,
        kt = kt,
        chain_id = chain_id,
        light_weight = component_weights[["light"]]
      )
      path_idx <- refreshed$path_idx
      repair_iteration_count <-
        repair_iteration_count + as.integer(refreshed$n_repaired_sites > 0L)
      repaired_site_count <-
        repaired_site_count + refreshed$n_repaired_sites
      max_repaired_sites <- max(
        max_repaired_sites,
        refreshed$n_repaired_sites
      )
      if (!is.null(route_edge_distance)) {
        route_distance <- sampling_path_route_distance_state(
          path_idx,
          kt,
          route_prior
        )
        route_edge_distance <- route_distance$edge
        route_interval_distance <- route_distance$interval
      }
    }

    block_iteration <-
      block_interval > 0L &&
      iter_i %% block_interval == 0L &&
      nstap > 1L

    if (block_iteration) {
      # Freeze the current stay structure and update each block in one go.
      block_length <- rle(path_idx)$lengths
      block_end <- cumsum(block_length)
      block_start <- block_end - block_length + 1L

      for (block_i in seq_along(block_start)) {
        block_range <- block_start[block_i]:block_end[block_i]
        if (any(block_range %in% long_period_i)) {
          next
        }

        new_idx <- sampling_path_update_block(
          path_idx = path_idx,
          block_start = block_start[block_i],
          block_end = block_end[block_i],
          lk = lk,
          kt = kt,
          chain_id = chain_id,
          iter_i = iter_i,
          component_weights = component_weights,
          route_prior = route_prior,
          route_edge_distance = route_edge_distance,
          route_interval_distance = route_interval_distance
        )

        path_idx[block_range] <- new_idx
        if (!is.null(route_edge_distance)) {
          changed_edge <- seq.int(
            max(1L, block_start[block_i] - 1L),
            min(nstap - 1L, block_end[block_i])
          )
          route_distance <- sampling_path_update_route_distance(
            path_idx,
            changed_edge,
            kt,
            route_prior,
            route_edge_distance,
            route_interval_distance
          )
          route_edge_distance <- route_distance$edge
          route_interval_distance <- route_distance$interval
        }
      }
    } else {
      site_iteration_i <- site_iteration_i + 1L
      site_order <- if (site_iteration_i %% 2L == 1L) {
        seq_len(nstap)
      } else {
        rev(seq_len(nstap))
      }

      for (stap_i in setdiff(site_order, long_period_i)) {
        path_idx[stap_i] <- sampling_path_update_site(
          path_idx = path_idx,
          stap_i = stap_i,
          lk = lk,
          kt = kt,
          chain_id = chain_id,
          component_weights = component_weights,
          route_prior = route_prior,
          route_edge_distance = route_edge_distance,
          route_interval_distance = route_interval_distance,
          iter_i = iter_i
        )
        if (!is.null(route_edge_distance)) {
          changed_edge <- intersect(
            c(stap_i - 1L, stap_i),
            seq_len(nstap - 1L)
          )
          route_distance <- sampling_path_update_route_distance(
            path_idx,
            changed_edge,
            kt,
            route_prior,
            route_edge_distance,
            route_interval_distance
          )
          route_edge_distance <- route_distance$edge
          route_interval_distance <- route_distance$interval
        }
      }
    }

    if (iter_i > warmup && ((iter_i - warmup) %% thin == 0L)) {
      save_idx <- save_idx + 1L
      samples[save_idx, ] <- path_idx
    }

    if (
      !quiet &&
        !is.null(progress_fn) &&
        (iter_i %% progress_every == 0L || iter_i == iter)
    ) {
      progress_fn(iter_i - last_progress)
      last_progress <- iter_i
    }
  }

  attr(samples, "repair_diagnostics") <- list(
    n_repair_iterations = repair_iteration_count,
    n_repaired_sites = repaired_site_count,
    max_repaired_sites = max_repaired_sites
  )
  samples
}

#' Update One Stay Block
#'
#' @description
#' Internal block Gibbs update used by `sampling_path_run_chain()`. It freezes
#' the current stay structure, aggregates the light likelihood over one stay
#' block, and samples one shared location for the whole block.
#'
#' @param path_idx integer vector of current path indices.
#' @param block_start integer start position of the stay block.
#' @param block_end integer end position of the stay block.
#' @param lk list of per-stap likelihood objects.
#' @param kt precomputed transition lookup object.
#' @param chain_id optional chain identifier used in diagnostics.
#' @param iter_i optional Gibbs iteration identifier used in diagnostics.
#' @param route_prior prepared route-prior parameters, or `NULL`.
#' @param route_edge_distance current distance for every path edge.
#' @param route_interval_distance current cumulative distance for every
#'   long-period interval.
#'
#' @return Integer index chosen for the updated stay block.
#' @noRd
sampling_path_update_block <- function(
  path_idx,
  block_start,
  block_end,
  lk,
  kt,
  chain_id = NULL,
  iter_i = NULL,
  component_weights = c(light = 1, movement = 1, route = 1),
  route_prior = NULL,
  route_edge_distance = NULL,
  route_interval_distance = NULL
) {
  prev_idx <- if (block_start > 1L) path_idx[block_start - 1L] else NULL
  next_idx <- if (block_end < length(path_idx)) {
    path_idx[block_end + 1L]
  } else {
    NULL
  }

  # Aggregate the block likelihood on the intersection of its fixed supports.
  block_states <- lk[block_start:block_end]
  block_idx <- Reduce(intersect, lapply(block_states, `[[`, "idx"))
  block_log_prob <- sampling_path_weight_log_prob(
    Reduce(
      `+`,
      lapply(block_states, function(state) {
        state$log_prob[match(block_idx, state$idx)]
      })
    ),
    component_weights[["light"]]
  )

  state <- build_transition_state(
    lk_state = list(
      stap_id = glue::glue(
        "{lk[[block_start]]$stap_id}:{lk[[block_end]]$stap_id}"
      ),
      idx = block_idx,
      log_prob = block_log_prob
    )
  )

  state_prev <- state
  if (!is.null(prev_idx)) {
    state_prev$log_prob <- state_prev$log_prob +
      sampling_path_weight_log_prob(
        sampling_path_transition_log_prob(
          prev_idx,
          state_prev$idx,
          kt
        ),
        component_weights[["movement"]]
      )
  }

  state_next <- state_prev
  if (!is.null(next_idx)) {
    state_next$log_prob <- state_next$log_prob +
      sampling_path_weight_log_prob(
        sampling_path_transition_log_prob(
          state_next$idx,
          next_idx,
          kt
        ),
        component_weights[["movement"]]
      )
  }

  state_next$log_prob <- state_next$log_prob +
    sampling_path_route_log_prior(
      candidate_idx = state_next$idx,
      path_idx = path_idx,
      replace_i = block_start:block_end,
      kt = kt,
      route_prior = route_prior,
      route_edge_distance = route_edge_distance,
      route_interval_distance = route_interval_distance
    )

  max_log <- max(state_next$log_prob)

  if (!is.finite(max_log)) {
    ctx <- if (!is.null(iter_i)) {
      glue::glue(" (chain {chain_id}, iter {iter_i})")
    } else {
      ""
    }
    cli::cli_abort(
      "Block transition probability is zero for stap_id {.val {state$stap_id}}{ctx}."
    )
  }

  prob <- exp(state_next$log_prob - max_log)
  sampled_idx <- state_next$idx[sampling_path_sample_prob(prob)]

  sampled_idx
}

#' Update One Site
#'
#' @description
#' Internal local update used by `sampling_path_run_chain()`. It reads the
#' fixed candidate support for one period, applies neighbouring transition
#' weights, and samples a new index.
#'
#' @param path_idx integer vector of current path indices.
#' @param stap_i integer stap position to update.
#' @param lk list of per-stap likelihood objects.
#' @param kt precomputed transition lookup object.
#' @param chain_id optional chain identifier used in diagnostics.
#' @param iter_i optional Gibbs iteration identifier used in diagnostics.
#' @return Integer index chosen for the updated stap.
#' @noRd
sampling_path_update_site <- function(
  path_idx,
  stap_i,
  lk,
  kt,
  chain_id = NULL,
  component_weights = c(light = 1, movement = 1, route = 1),
  route_prior = NULL,
  route_edge_distance = NULL,
  route_interval_distance = NULL,
  iter_i = NULL
) {
  # Read the fixed light-likelihood support for this stap.
  state <- build_transition_state(lk[[stap_i]])
  state$log_prob <- sampling_path_weight_log_prob(
    state$log_prob,
    component_weights[["light"]]
  )

  state$log_prob <- state$log_prob +
    sampling_path_weight_log_prob(
      sampling_path_local_log_prior(
        candidate_idx = state$idx,
        path_idx = path_idx,
        stap_i = stap_i,
        kt = kt
      ),
      component_weights[["movement"]]
    )
  state$log_prob <- state$log_prob +
    sampling_path_route_log_prior(
      candidate_idx = state$idx,
      path_idx = path_idx,
      replace_i = stap_i,
      kt = kt,
      route_prior = route_prior,
      route_edge_distance = route_edge_distance,
      route_interval_distance = route_interval_distance
    )

  max_log <- max(state$log_prob)

  if (!is.finite(max_log)) {
    ctx <- if (!is.null(iter_i)) {
      glue::glue(" (chain {chain_id}, iter {iter_i})")
    } else {
      ""
    }
    cli::cli_abort(
      "Transition probability is zero for stap_id {.val {state$stap_id}}{ctx}."
    )
  }

  prob <- exp(state$log_prob - max_log)
  sampled_idx <- state$idx[sampling_path_sample_prob(prob)]

  sampled_idx
}

#' Exact Local Residence Prior for One Stap Update
#'
#' @description
#' Recalculate the transition prior only across the neighbouring residence runs
#' whose lengths can change when `stap_i` is reassigned. Ordinary candidates
#' share one residence pattern and retain vectorized movement-kernel lookup;
#' only candidates equal to a neighbour require separate scalar evaluation.
#'
#' @return Numeric log prior, one value per candidate index.
#' @noRd
sampling_path_local_log_prior <- function(
  candidate_idx,
  path_idx,
  stap_i,
  kt
) {
  nstap <- length(path_idx)
  if (nstap == 1L) {
    return(rep(0, length(candidate_idx)))
  }

  prev_idx <- if (stap_i > 1L) path_idx[stap_i - 1L] else NULL
  next_idx <- if (stap_i < nstap) path_idx[stap_i + 1L] else NULL
  left_start <- if (stap_i > 1L) {
    i <- stap_i - 1L
    while (i > 1L && path_idx[i - 1L] == path_idx[i]) {
      i <- i - 1L
    }
    i
  } else {
    1L
  }
  right_end <- if (stap_i < nstap) {
    i <- stap_i + 1L
    while (i < nstap && path_idx[i + 1L] == path_idx[i]) {
      i <- i + 1L
    }
    i
  } else {
    nstap
  }
  edge_start <- left_start
  edge_end <- min(nstap - 1L, right_end)
  special_idx <- unique(c(prev_idx, next_idx))
  special <- candidate_idx %in% special_idx
  out <- rep(-Inf, length(candidate_idx))

  if (!all(special)) {
    placeholder <- candidate_idx[which(!special)[1L]]
    out[!special] <- sampling_path_local_path_log_prior(
      candidate_idx = placeholder,
      path_idx = path_idx,
      stap_i = stap_i,
      edge_start = edge_start,
      edge_end = edge_end,
      kt = kt,
      omit_kernel_edges = c(stap_i - 1L, stap_i)
    )

    if (!is.null(prev_idx)) {
      out[!special] <- out[!special] +
        sampling_path_transition_log_prob(
          prev_idx,
          candidate_idx[!special],
          kt
        )
    }
    if (!is.null(next_idx)) {
      out[!special] <- out[!special] +
        sampling_path_transition_log_prob(
          candidate_idx[!special],
          next_idx,
          kt
        )
    }
  }

  for (candidate_i in which(special)) {
    out[candidate_i] <- sampling_path_local_path_log_prior(
      candidate_idx = candidate_idx[candidate_i],
      path_idx = path_idx,
      stap_i = stap_i,
      edge_start = edge_start,
      edge_end = edge_end,
      kt = kt
    )
  }

  out
}

#' Transition Prior Across One Affected Residence Window
#'
#' @return Numeric log probability.
#' @noRd
sampling_path_local_path_log_prior <- function(
  candidate_idx,
  path_idx,
  stap_i,
  edge_start,
  edge_end,
  kt,
  omit_kernel_edges = integer()
) {
  local_path <- path_idx
  local_path[stap_i] <- candidate_idx
  residence <- 0L
  log_prob <- 0

  for (edge_i in edge_start:edge_end) {
    source_idx <- local_path[edge_i]
    destination_idx <- local_path[edge_i + 1L]
    move_prob <- kt$move_probability[residence + 1L]

    if (source_idx == destination_idx) {
      log_prob <- log_prob + log1p(-move_prob)
      residence <- residence + 1L
    } else {
      log_prob <- log_prob + log(move_prob)
      if (!edge_i %in% omit_kernel_edges) {
        log_prob <- log_prob +
          sampling_path_transition_log_prob(
            source_idx,
            destination_idx,
            kt
          )
      }
      residence <- 0L
    }
  }

  log_prob
}

#' Initialize a Sampling Path
#'
#' @description
#' Internal initialization routine used by `sampling_path_run_chain()`. It draws
#' one location per stap under a sharpened light likelihood, then repairs only
#' adjacent pairs that fall outside the movement kernel.
#'
#' @param lk list of per-stap likelihood objects.
#' @param kt precomputed transition lookup object.
#' @param chain_id optional chain identifier used in diagnostics.
#' @param init_power exponent applied to the light log-probabilities during the
#'   initial draw.
#'
#' @return Integer vector of initialized path indices.
#' @noRd
sampling_path_initialize <- function(
  lk,
  kt,
  chain_id = NULL,
  init_power = 4
) {
  lk_init <- lapply(lk, function(state) {
    state$log_prob <- init_power * state$log_prob
    state
  })
  path_idx <- vapply(lk_init, sampling_path_draw_light, integer(1))
  sampling_path_repair(
    path_idx = path_idx,
    fixed_i = integer(),
    lk = lk_init,
    kt = kt,
    chain_id = chain_id
  )$path_idx
}

#' Refresh Unknown Long-Periods
#'
#' @description
#' Draw each unknown long-period directly from its retained light likelihood,
#' retain the previous daily path where it remains movement-feasible, and
#' repair only daily locations needed to restore hard movement support.
#'
#' @return A list containing the repaired path and repair counts.
#' @noRd
sampling_path_refresh_long_periods <- function(
  path_idx,
  long_period_i,
  lk,
  kt,
  chain_id = NULL,
  light_weight = 1
) {
  path_idx[long_period_i] <- vapply(
    lk[long_period_i],
    sampling_path_draw_light,
    integer(1),
    weight = light_weight
  )
  drawn_path_idx <- path_idx
  repaired <- sampling_path_repair(
    path_idx = path_idx,
    fixed_i = long_period_i,
    lk = lk,
    kt = kt,
    chain_id = chain_id,
    changed_i = long_period_i
  )
  repaired_i <- setdiff(
    which(repaired$path_idx != drawn_path_idx),
    long_period_i
  )

  list(
    path_idx = repaired$path_idx,
    n_repaired_sites = length(repaired_i),
    repaired_i = repaired_i
  )
}

#' Repair Hard Movement-Support Violations
#'
#' @description
#' Retain valid locations and expand local repairs across adjacent daily
#' periods until the path reconnects. Fixed long-periods are never altered.
#' This establishes a feasible state; it is not a posterior update.
#'
#' @return A list containing the repaired path and changed period indices.
#' @noRd
sampling_path_repair <- function(
  path_idx,
  fixed_i,
  lk,
  kt,
  chain_id = NULL,
  max_passes = NULL,
  changed_i = NULL
) {
  nstap <- length(lk)
  original_path_idx <- path_idx
  fixed <- seq_len(nstap) %in% fixed_i

  if (nstap <= 1L) {
    return(list(path_idx = path_idx, changed_i = integer()))
  }

  if (is.null(max_passes)) {
    max_passes <- max(10L, 2L * nstap)
  }

  edge_id <- seq_len(nstap - 1L)
  edge_ok <- rep(TRUE, nstap - 1L)
  edge_check <- if (is.null(changed_i)) {
    edge_id
  } else {
    intersect(c(changed_i - 1L, changed_i), edge_id)
  }
  edge_ok[edge_check] <- sampling_path_edges_valid(path_idx, kt, edge_check)

  repair_position <- function(stap_i, required_side) {
    repaired_idx <- sampling_path_repair_stap(
      path_idx = path_idx,
      stap_i = stap_i,
      lk = lk,
      kt = kt,
      use_prev = stap_i > 1L,
      use_next = stap_i < nstap
    )

    if (is.na(repaired_idx)) {
      repaired_idx <- sampling_path_repair_stap(
        path_idx = path_idx,
        stap_i = stap_i,
        lk = lk,
        kt = kt,
        use_prev = required_side == "previous",
        use_next = required_side == "next"
      )
    }

    repaired_idx
  }

  for (pass_i in seq_len(max_passes)) {
    if (all(edge_ok)) {
      return(list(
        path_idx = path_idx,
        changed_i = which(path_idx != original_path_idx)
      ))
    }

    if (pass_i %% 2L == 1L) {
      for (edge_i in edge_id) {
        if (edge_ok[edge_i]) {
          next
        }

        stap_i <- if (!fixed[edge_i + 1L]) edge_i + 1L else edge_i
        if (fixed[stap_i]) {
          next
        }
        required_side <- if (stap_i == edge_i + 1L) "previous" else "next"
        repaired_idx <- repair_position(stap_i, required_side)

        if (!is.na(repaired_idx)) {
          path_idx[stap_i] <- repaired_idx
          changed_edge <- intersect(c(stap_i - 1L, stap_i), edge_id)
          edge_ok[changed_edge] <- sampling_path_edges_valid(
            path_idx,
            kt,
            changed_edge
          )
        }
      }
    } else {
      for (edge_i in rev(edge_id)) {
        if (edge_ok[edge_i]) {
          next
        }

        stap_i <- if (!fixed[edge_i]) edge_i else edge_i + 1L
        if (fixed[stap_i]) {
          next
        }
        required_side <- if (stap_i == edge_i) "next" else "previous"
        repaired_idx <- repair_position(stap_i, required_side)

        if (!is.na(repaired_idx)) {
          path_idx[stap_i] <- repaired_idx
          changed_edge <- intersect(c(stap_i - 1L, stap_i), edge_id)
          edge_ok[changed_edge] <- sampling_path_edges_valid(
            path_idx,
            kt,
            changed_edge
          )
        }
      }
    }
  }

  bad_edge <- which(!edge_ok)[1]
  cli::cli_abort(
    "Could not initialize a movement-valid path near stap_id {.val {lk[[bad_edge]]$stap_id}} -> {.val {lk[[bad_edge + 1L]]$stap_id}}."
  )
}

#' Draw One Cell from a Retained Light Likelihood
#'
#' @return Integer grid-cell index.
#' @noRd
sampling_path_draw_light <- function(state, weight = 1) {
  log_prob <- sampling_path_weight_log_prob(state$log_prob, weight)
  max_log <- max(log_prob)
  prob <- exp(log_prob - max_log)
  state$idx[sampling_path_sample_prob(prob)]
}

#' Weight a Log-Probability Without Relaxing Hard Support
#'
#' @return Numeric log-probabilities.
#' @noRd
sampling_path_weight_log_prob <- function(log_prob, weight) {
  out <- weight * log_prob
  out[!is.finite(log_prob)] <- -Inf
  out
}

#' Prepare Span-Level Route Prior
#'
#' @return A route-prior list with sampler row indices, or `NULL`.
#' @noRd
sampling_path_prepare_route_prior <- function(
  tag,
  stap_ids,
  route_prior,
  weight = 1,
  detour = 1
) {
  if (is.null(route_prior)) {
    return(NULL)
  }

  parameter_names <- c(
    "intercept",
    "duration_coefficient",
    "distance_coefficient",
    "duration_center",
    "distance_center",
    "positive_link_scale",
    "residual_sd"
  )
  assertthat::assert_that(
    all(parameter_names %in% names(route_prior)),
    all(is.finite(unlist(route_prior[parameter_names]))),
    route_prior$residual_sd > 0,
    is.finite(weight),
    weight > 0,
    is.finite(detour),
    detour >= 0
  )
  route_prior$weight <- weight
  route_prior$detour <- detour

  stap_row <- match(stap_ids, tag$stap$stap_id)
  stap0 <- tag$stap$stap0 %||% rep(FALSE, nrow(tag$stap))
  long_period_i <- which(stap0[stap_row])
  if (length(long_period_i) < 2L) {
    return(NULL)
  }

  route_prior$intervals <- tibble::tibble(
    start_i = long_period_i[-length(long_period_i)],
    end_i = long_period_i[-1L]
  )
  start_row <- stap_row[route_prior$intervals$start_i]
  end_row <- stap_row[route_prior$intervals$end_i]
  interval_start <- tag$stap$end[start_row]
  interval_end <- tag$stap$start[end_row]
  route_prior$intervals$duration_days <- as.numeric(
    difftime(interval_end, interval_start, units = "days")
  )
  route_prior$intervals$duration_centered <-
    log1p(route_prior$intervals$duration_days) -
    route_prior$duration_center
  route_prior$interval_edges <- Map(
    seq.int,
    route_prior$intervals$start_i,
    route_prior$intervals$end_i - 1L
  )
  route_prior$edge_interval <- rep(NA_integer_, length(stap_ids) - 1L)
  for (interval_i in seq_along(route_prior$interval_edges)) {
    route_prior$edge_interval[route_prior$interval_edges[[interval_i]]] <-
      interval_i
  }
  route_prior
}

#' Initialize Cached Route Distances
#'
#' @return A list containing edge and interval distances, or two `NULL` values.
#' @noRd
sampling_path_route_distance_state <- function(path_idx, kt, route_prior) {
  if (is.null(route_prior) || length(path_idx) <= 1L) {
    return(list(edge = NULL, interval = NULL))
  }

  edge_distance <- sampling_path_index_distance_km(
    path_idx[-length(path_idx)],
    path_idx[-1L],
    kt
  )
  interval_distance <- vapply(
    route_prior$interval_edges,
    function(edge_i) sum(edge_distance[edge_i]),
    numeric(1)
  )
  list(edge = edge_distance, interval = interval_distance)
}

#' Update Cached Route Distances After a Local Path Change
#'
#' @return Updated edge and interval distances.
#' @noRd
sampling_path_update_route_distance <- function(
  path_idx,
  changed_edge,
  kt,
  route_prior,
  edge_distance,
  interval_distance
) {
  edge_distance[changed_edge] <- sampling_path_index_distance_km(
    path_idx[changed_edge],
    path_idx[changed_edge + 1L],
    kt
  )
  affected_interval <- unique(route_prior$edge_interval[changed_edge])
  affected_interval <- affected_interval[!is.na(affected_interval)]
  interval_distance[affected_interval] <- vapply(
    route_prior$interval_edges[affected_interval],
    function(edge_i) sum(edge_distance[edge_i]),
    numeric(1)
  )
  list(edge = edge_distance, interval = interval_distance)
}

#' Evaluate Span-Level Route Prior for One Local Update
#'
#' @return Numeric log-prior vector, one value per candidate index.
#' @noRd
sampling_path_route_log_prior <- function(
  candidate_idx,
  path_idx,
  replace_i,
  kt,
  route_prior,
  route_edge_distance = NULL,
  route_interval_distance = NULL
) {
  if (is.null(route_prior)) {
    return(rep(0, length(candidate_idx)))
  }

  out <- rep(0, length(candidate_idx))
  affected_interval <- which(
    route_prior$intervals$start_i <= max(replace_i) &
      route_prior$intervals$end_i >= min(replace_i)
  )

  for (interval_i in affected_interval) {
    interval <- route_prior$intervals[interval_i, ]
    edge_i <- route_prior$interval_edges[[interval_i]]
    current_edge_distance <- if (is.null(route_edge_distance)) {
      sampling_path_index_distance_km(
        path_idx[edge_i],
        path_idx[edge_i + 1L],
        kt
      )
    } else {
      route_edge_distance[edge_i]
    }
    changed_edge <- edge_i[
      edge_i %in% replace_i | (edge_i + 1L) %in% replace_i
    ]
    current_interval_distance <- if (is.null(route_interval_distance)) {
      sum(current_edge_distance)
    } else {
      route_interval_distance[interval_i]
    }
    route_distance <- rep(
      current_interval_distance -
        sum(current_edge_distance[match(changed_edge, edge_i)]),
      length(candidate_idx)
    )

    for (edge in changed_edge) {
      from_idx <- if (edge %in% replace_i) {
        candidate_idx
      } else {
        path_idx[edge]
      }
      to_idx <- if ((edge + 1L) %in% replace_i) {
        candidate_idx
      } else {
        path_idx[edge + 1L]
      }
      route_distance <- route_distance +
        sampling_path_index_distance_km(from_idx, to_idx, kt)
    }

    start_idx <- if (interval$start_i %in% replace_i) {
      candidate_idx
    } else {
      path_idx[interval$start_i]
    }
    end_idx <- if (interval$end_i %in% replace_i) {
      candidate_idx
    } else {
      path_idx[interval$end_i]
    }
    direct_distance <- sampling_path_index_distance_km(
      start_idx,
      end_idx,
      kt
    )
    apply_prior <- direct_distance >= route_prior$min_direct_distance_km
    log_ratio <- log(pmax(route_distance / direct_distance, 1))
    linear_predictor <-
      route_prior$intercept +
      route_prior$duration_coefficient * interval$duration_centered +
      route_prior$distance_coefficient *
        (log(direct_distance) - route_prior$distance_center)
    mean_log_ratio <- sampling_path_positive_linear_link(
      linear_predictor,
      route_prior$positive_link_scale
    )
    mean_log_ratio <- log1p(
      (route_prior$detour %||% 1) * expm1(mean_log_ratio)
    )
    out[apply_prior] <- out[apply_prior] +
      stats::dnorm(
        log_ratio[apply_prior],
        mean = mean_log_ratio[apply_prior],
        sd = route_prior$residual_sd,
        log = TRUE
      )
  }

  out * (route_prior$weight %||% 1)
}

#' Positive-Linear Link
#'
#' @return Positive numeric vector that is approximately linear above zero.
#' @noRd
sampling_path_positive_linear_link <- function(x, scale) {
  out <- x
  smooth <- scale * x <= 30
  out[smooth] <- log1p(exp(scale * x[smooth])) / scale
  out
}

#' Great-Circle Distance Between Grid Indices
#'
#' @return Numeric distance in kilometres.
#' @noRd
sampling_path_index_distance_km <- function(from_idx, to_idx, kt) {
  n <- max(length(from_idx), length(to_idx))
  from_idx <- rep_len(from_idx, n)
  to_idx <- rep_len(to_idx, n)
  from_row <- (from_idx - 1L) %% kt$nrow + 1L
  to_row <- (to_idx - 1L) %% kt$nrow + 1L
  from_lat <- kt$lat_rad[from_row]
  to_lat <- kt$lat_rad[to_row]
  from_lon <- kt$lon_rad[(from_idx - 1L) %/% kt$nrow + 1L]
  to_lon <- kt$lon_rad[(to_idx - 1L) %/% kt$nrow + 1L]
  dlat <- to_lat - from_lat
  dlon <- to_lon - from_lon
  a <- sin(dlat / 2)^2 +
    kt$cos_lat[from_row] * kt$cos_lat[to_row] * sin(dlon / 2)^2
  2 * 6378.137 * atan2(sqrt(pmin(a, 1)), sqrt(pmax(1 - a, 0)))
}

#' Repair One Initialized Stap
#'
#' @description
#' Draw one replacement location for an initialized path, using the sharpened
#' likelihood while requiring the requested neighbouring edges to be valid.
#'
#' @return Integer grid index, or `NA_integer_` when no local repair is possible.
#' @noRd
sampling_path_repair_stap <- function(
  path_idx,
  stap_i,
  lk,
  kt,
  use_prev,
  use_next
) {
  prev_idx <- if (use_prev && stap_i > 1L) path_idx[stap_i - 1L] else NULL
  next_idx <- if (use_next && stap_i < length(path_idx)) {
    path_idx[stap_i + 1L]
  } else {
    NULL
  }
  state <- build_transition_state(lk[[stap_i]])

  keep <- is.finite(state$log_prob)

  if (!is.null(prev_idx)) {
    keep <- keep &
      is.finite(sampling_path_transition_log_prob(
        prev_idx,
        state$idx,
        kt
      ))
  }

  if (!is.null(next_idx)) {
    keep <- keep &
      is.finite(sampling_path_transition_log_prob(
        state$idx,
        next_idx,
        kt
      ))
  }

  if (!any(keep)) {
    return(NA_integer_)
  }

  log_prob <- state$log_prob[keep]
  max_log <- max(log_prob)
  prob <- exp(log_prob - max_log)
  state$idx[keep][sampling_path_sample_prob(prob)]
}

#' Draw One Index from Unnormalized Probabilities
#'
#' @return Integer scalar.
#' @noRd
sampling_path_sample_prob <- function(prob) {
  cumulative_prob <- cumsum(prob)
  findInterval(
    stats::runif(1L, max = cumulative_prob[length(cumulative_prob)]),
    cumulative_prob
  ) +
    1L
}

#' Check Movement Support Across Selected Path Edges
#'
#' @return Logical vector aligned with `edge_i`.
#' @noRd
sampling_path_edges_valid <- function(path_idx, kt, edge_i) {
  from_idx <- path_idx[edge_i]
  to_idx <- path_idx[edge_i + 1L]
  valid <- !is.na(from_idx) & !is.na(to_idx) & kt$struct_mask[to_idx]
  stay <- valid & from_idx == to_idx
  move <- valid & !stay
  valid[move] <- is.finite(sampling_path_transition_log_prob(
    from_idx[move],
    to_idx[move],
    kt
  ))
  valid
}

#' Build One Local Transition State
#'
#' @description
#' Internal helper that constructs the fixed light-support state for one stap.
#'
#' @param lk_state per-stap likelihood object.
#' @return A list with candidate indices, log-probabilities, and `stap_id`.
#' @noRd
build_transition_state <- function(lk_state) {
  list(
    idx = lk_state$idx,
    log_prob = lk_state$log_prob,
    stap_id = lk_state$stap_id
  )
}

#' Evaluate Directed Movement Probabilities
#'
#' @description
#' Evaluate the normalized movement component from `source_idx` to
#' `destination_idx`. Either argument can be scalar and is recycled to the
#' length of the other.
#'
#' @return Numeric vector of log probabilities.
#' @noRd
sampling_path_transition_log_prob <- function(
  source_idx,
  destination_idx,
  kt
) {
  kernel <- kt$movement_kernel
  n <- max(length(source_idx), length(destination_idx))
  source_idx <- rep_len(source_idx, n)
  destination_idx <- rep_len(destination_idx, n)
  source_row <- kt$cell_row[source_idx]
  drow <- kt$cell_row[destination_idx] - source_row
  dcol <- kt$cell_col[destination_idx] - kt$cell_col[source_idx]
  out <- rep(-Inf, n)
  move_pos <- source_idx != destination_idx &
    abs(drow) <= kernel$max_drow &
    abs(dcol) <= kernel$max_dcol &
    kt$struct_mask[destination_idx]

  move_i <- which(move_pos)
  lookup_idx <- source_row[move_i] +
    (drow[move_i] + kernel$max_drow) * kt$nrow +
    (dcol[move_i] + kernel$max_dcol) * kt$nrow * kernel$lookup_ndrow
  move_log_prob <- kernel$move_log_lookup[lookup_idx]
  finite_move <- is.finite(move_log_prob)
  move_i <- move_i[finite_move]
  out[move_i] <- move_log_prob[finite_move] -
    kernel$move_log_norm[source_idx[move_i]]

  out
}
