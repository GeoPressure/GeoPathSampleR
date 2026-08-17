sampling_test_lookup <- function(ncol = 4L) {
  max_dcol <- ncol - 1L
  move_log_lookup <- rep(0, 2L * max_dcol + 1L)
  move_log_lookup[max_dcol + 1L] <- -Inf

  list(
    move_probability = rep(0.3, 10),
    movement_kernel = list(
      move_log_lookup = move_log_lookup,
      move_log_norm = rep(log(ncol - 1L), ncol),
      max_drow = 0L,
      max_dcol = max_dcol,
      lookup_ndrow = 1L
    ),
    nrow = 1L,
    ncol = ncol,
    struct_mask = rep(TRUE, ncol),
    cell_row = rep(1L, ncol),
    cell_col = seq_len(ncol),
    lat = 0,
    lon = seq_len(ncol),
    lat_rad = 0,
    lon_rad = seq_len(ncol) * pi / 180,
    cos_lat = 1
  )
}

test_that("batched movement normalization matches direct summation", {
  nrow <- 3L
  ncol <- 5L
  ncell <- nrow * ncol
  struct_mask <- c(
    TRUE,
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    TRUE,
    TRUE,
    TRUE
  )
  cell_row <- ((seq_len(ncell) - 1L) %% nrow) + 1L
  cell_col <- ((seq_len(ncell) - 1L) %/% nrow) + 1L
  drow <- list(
    c(0L, 0L, 1L),
    c(-1L, 0L, 1L),
    c(-1L, 0L, 0L)
  )
  dcol <- list(
    c(-1L, 1L, 0L),
    c(0L, -1L, 1L),
    c(0L, -1L, 1L)
  )
  kernels <- lapply(drow, function(x) log(seq_along(x) + 1))
  support_idx <- which(struct_mask)

  actual <- sampling_path_movement_log_norm(
    support_idx,
    struct_mask,
    drow,
    dcol,
    kernels,
    cell_row,
    cell_col,
    nrow,
    ncol
  )
  expected <- rep(NA_real_, ncell)

  for (source_idx in support_idx) {
    source_row <- cell_row[source_idx]
    source_col <- cell_col[source_idx]
    destination_row <- source_row + drow[[source_row]]
    destination_col <- source_col + dcol[[source_row]]
    keep <- destination_row >= 1L &
      destination_row <= nrow &
      destination_col >= 1L &
      destination_col <= ncol &
      (drow[[source_row]] != 0L | dcol[[source_row]] != 0L)
    destination_idx <- destination_row[keep] +
      (destination_col[keep] - 1L) * nrow
    weight <- exp(kernels[[source_row]][keep])
    total <- sum(weight[struct_mask[destination_idx]])
    expected[source_idx] <- if (total > 0) log(total) else -Inf
  }

  expect_equal(actual, expected)
})

test_that("local movement prior matches full-path probability ratios", {
  kt <- sampling_test_lookup()
  path_idx <- c(1L, 1L, 2L, 3L)
  candidate_idx <- c(1L, 2L, 4L)

  brute_log_prior <- function(candidate) {
    local_path <- path_idx
    local_path[2L] <- candidate
    residence <- 0L
    out <- 0

    for (edge_i in seq_len(length(local_path) - 1L)) {
      if (local_path[edge_i] == local_path[edge_i + 1L]) {
        out <- out + log1p(-kt$move_probability[residence + 1L])
        residence <- residence + 1L
      } else {
        out <- out +
          log(kt$move_probability[residence + 1L]) +
          sampling_path_transition_log_prob(
            local_path[edge_i],
            local_path[edge_i + 1L],
            kt
          )
        residence <- 0L
      }
    }
    out
  }

  local <- sampling_path_local_log_prior(candidate_idx, path_idx, 2L, kt)
  brute <- vapply(candidate_idx, brute_log_prior, numeric(1))

  expect_equal(local - local[1L], brute - brute[1L])
})

test_that("block updates use the intersection of sparse likelihood supports", {
  kt <- sampling_test_lookup()
  lk <- list(
    list(stap_id = 1L, idx = 1:3, log_prob = log(c(0.2, 0.5, 0.3))),
    list(stap_id = 2L, idx = 2:4, log_prob = log(c(0.4, 0.4, 0.2))),
    list(stap_id = 3L, idx = 1:4, log_prob = rep(log(0.25), 4))
  )
  path_idx <- c(2L, 2L, 4L)

  sampled <- vapply(
    1:20,
    function(seed) {
      withr::with_seed(
        seed,
        sampling_path_update_block(
          path_idx = path_idx,
          block_start = 1L,
          block_end = 2L,
          lk = lk,
          kt = kt
        )
      )
    },
    integer(1)
  )

  expect_true(all(sampled %in% c(2L, 3L)))

  block_idx <- c(2L, 3L)
  log_probability <- lk[[1L]]$log_prob[match(block_idx, lk[[1L]]$idx)] +
    lk[[2L]]$log_prob[match(block_idx, lk[[2L]]$idx)] +
    log(kt$move_probability[1L]) +
    sampling_path_transition_log_prob(block_idx, path_idx[3L], kt)
  probability <- exp(log_probability - max(log_probability))
  expected <- withr::with_seed(
    4,
    block_idx[sampling_path_sample_prob(probability)]
  )
  actual <- withr::with_seed(
    4,
    sampling_path_update_block(
      path_idx = path_idx,
      block_start = 1L,
      block_end = 2L,
      lk = lk,
      kt = kt
    )
  )
  expect_equal(actual, expected)
})

test_that("component weights scale light and movement contributions", {
  kt <- sampling_test_lookup()
  lk <- list(
    list(stap_id = 1L, idx = 1:3, log_prob = log(c(0.2, 0.5, 0.3)))
  )
  path_idx <- 2L

  expected_uniform <- withr::with_seed(
    8,
    lk[[1]]$idx[sampling_path_sample_prob(rep(1, 3))]
  )
  actual_uniform <- withr::with_seed(
    8,
    sampling_path_update_site(
      path_idx = path_idx,
      stap_i = 1L,
      lk = lk,
      kt = kt,
      component_weights = c(light = 0, movement = 0, route = 0)
    )
  )
  expect_equal(actual_uniform, expected_uniform)

  light_weight <- 2
  expected_light <- withr::with_seed(
    8,
    lk[[1]]$idx[sampling_path_sample_prob(exp(light_weight * lk[[1]]$log_prob))]
  )
  actual_light <- withr::with_seed(
    8,
    sampling_path_update_site(
      path_idx = path_idx,
      stap_i = 1L,
      lk = lk,
      kt = kt,
      component_weights = c(light = light_weight, movement = 0, route = 0)
    )
  )
  expect_equal(actual_light, expected_light)
})

test_that("zero weights retain impossible states as hard exclusions", {
  expect_equal(
    sampling_path_weight_log_prob(c(-Inf, -2, 0), 0),
    c(-Inf, 0, 0)
  )
})

test_that("incremental edge checks agree with a complete path check", {
  kt <- sampling_test_lookup()
  kt$struct_mask[4L] <- FALSE
  path_idx <- c(1L, 1L, 3L, 4L)
  complete <- sampling_path_edges_valid(path_idx, kt, 1:3)
  changed_edges <- 2:3

  expect_equal(complete, c(TRUE, TRUE, FALSE))
  expect_equal(
    sampling_path_edges_valid(path_idx, kt, changed_edges),
    complete[changed_edges]
  )
})

test_that("cached route distances preserve the local route prior", {
  kt <- sampling_test_lookup()
  stap <- data.frame(
    stap_id = 1:4,
    start = as.POSIXct("2026-01-01", tz = "UTC") + 0:3 * 86400,
    end = as.POSIXct("2026-01-02", tz = "UTC") + 0:3 * 86400,
    stap0 = c(TRUE, FALSE, FALSE, TRUE)
  )
  route_prior <- sampling_path_prepare_route_prior(
    list(stap = stap),
    stap_ids = stap$stap_id,
    route_prior = sampling_path_route_model()
  )
  path_idx <- c(1L, 2L, 3L, 4L)
  route_distance <- sampling_path_route_distance_state(
    path_idx,
    kt,
    route_prior
  )
  candidate_idx <- 1:4

  uncached <- sampling_path_route_log_prior(
    candidate_idx,
    path_idx,
    2L,
    kt,
    route_prior
  )
  cached <- sampling_path_route_log_prior(
    candidate_idx,
    path_idx,
    2L,
    kt,
    route_prior,
    route_distance$edge,
    route_distance$interval
  )

  expect_equal(cached, uncached)

  path_updated <- path_idx
  path_updated[2] <- 4L
  updated_distance <- sampling_path_update_route_distance(
    path_updated,
    changed_edge = 1:2,
    kt,
    route_prior,
    route_distance$edge,
    route_distance$interval
  )
  fresh_distance <- sampling_path_route_distance_state(
    path_updated,
    kt,
    route_prior
  )

  expect_equal(updated_distance$edge, fresh_distance$edge)
  expect_equal(updated_distance$interval, fresh_distance$interval)
})

test_that("route weight and detour adjust the calibrated route term", {
  kt <- sampling_test_lookup()
  stap <- data.frame(
    stap_id = 1:4,
    start = as.POSIXct("2026-01-01", tz = "UTC") + 0:3 * 86400,
    end = as.POSIXct("2026-01-02", tz = "UTC") + 0:3 * 86400,
    stap0 = c(TRUE, FALSE, FALSE, TRUE)
  )
  path_idx <- c(1L, 2L, 3L, 4L)
  candidate_idx <- 2L

  route_default <- sampling_path_prepare_route_prior(
    list(stap = stap),
    stap$stap_id,
    sampling_path_route_model()
  )
  route_weighted <- sampling_path_prepare_route_prior(
    list(stap = stap),
    stap$stap_id,
    sampling_path_route_model(),
    weight = 2
  )
  route_direct <- sampling_path_prepare_route_prior(
    list(stap = stap),
    stap$stap_id,
    sampling_path_route_model(),
    detour = 0
  )

  default_value <- sampling_path_route_log_prior(
    candidate_idx,
    path_idx,
    2L,
    kt,
    route_default
  )
  weighted_value <- sampling_path_route_log_prior(
    candidate_idx,
    path_idx,
    2L,
    kt,
    route_weighted
  )
  direct_value <- sampling_path_route_log_prior(
    candidate_idx,
    path_idx,
    2L,
    kt,
    route_direct
  )

  expect_equal(weighted_value, 2 * default_value)
  expect_gt(direct_value, default_value)
})

test_that("route support clamps unsupported covariates to its boundary", {
  support <- list(
    duration_log = log1p(c(1, 10)),
    log_distance_lower = log(c(100, 100)),
    log_distance_upper = log(c(1000, 1000))
  )

  projected <- sampling_path_route_support_projection(
    5,
    c(50, 500, 1500),
    support
  )

  expect_equal(projected$duration_log, rep(log1p(5), 3))
  expect_equal(projected$distance_log, log(c(100, 500, 1000)))
  expect_equal(
    sampling_path_route_support_projection(20, 500, support)$duration_log,
    log1p(10)
  )
})

test_that("default route model includes the direct-distance calibration support", {
  route_model <- sampling_path_route_model()

  expect_equal(route_model$intercept, 0.144588068514403)
  expect_equal(route_model$distance_coefficient, -0.0419027867961805)
  expect_equal(route_model$min_direct_distance_km, 300)
  expect_length(route_model$support$duration_log, 80)
  expect_length(route_model$support$log_distance_lower, 80)
  expect_length(route_model$support$log_distance_upper, 80)
})
