test_that("path summaries retain stationary-period structure", {
  stap <- data.frame(
    stap_id = 1:3,
    start = as.POSIXct("2026-01-01", tz = "UTC") + 0:2 * 86400,
    end = as.POSIXct("2026-01-02", tz = "UTC") + 0:2 * 86400,
    stap0 = c(TRUE, FALSE, TRUE),
    known_lat = NA_real_,
    known_lon = NA_real_
  )
  paths <- expand.grid(j = 1:4, chain = 1:2, stap_id = stap$stap_id)
  paths <- paths[order(paths$chain, paths$j, paths$stap_id), ]
  paths$ind <- ifelse(paths$stap_id == 2, 2L, 1L)
  paths$lat <- ifelse(paths$ind == 1L, 45, 46)
  paths$lon <- ifelse(paths$ind == 1L, 5, 6)

  collapsed <- path_collapse(paths, stap)
  consensus <- path_summary(paths, stap, by = "consensus_stay")
  diagnostic_summary <- sampling_path_diagnostic_sample_summary(paths)

  expect_equal(nrow(collapsed), 24)
  expect_equal(consensus$stap_id_start, c(1L, 2L, 3L))
  expect_equal(consensus$lat, c(45, 46, 45))
  expect_equal(diagnostic_summary$n_stays, rep(3L, 8L))
  expect_equal(diagnostic_summary$n_unique_cells, rep(2L, 8L))
})

test_that("diagnostics return marginal summaries without report rendering", {
  stap <- data.frame(
    stap_id = 1:2,
    start = as.POSIXct("2026-01-01", tz = "UTC") + 0:1 * 86400,
    end = as.POSIXct("2026-01-02", tz = "UTC") + 0:1 * 86400,
    stap0 = c(TRUE, FALSE),
    known_lat = NA_real_,
    known_lon = NA_real_
  )
  paths <- expand.grid(j = 1:8, chain = 1:4, stap_id = stap$stap_id)
  paths <- paths[order(paths$chain, paths$j, paths$stap_id), ]
  paths$ind <- ifelse(paths$stap_id == 1L, 1L, 2L)
  paths$lat <- ifelse(paths$stap_id == 1L, 45, 46)
  paths$lon <- ifelse(paths$stap_id == 1L, 5, 6)

  diagnostic <- sampling_path_diagnostic(
    paths,
    tag = list(stap = stap, param = list(id = "test")),
    report = FALSE,
    quiet = TRUE
  )

  expect_equal(diagnostic$tag_id, "test")
  expect_equal(nrow(diagnostic$marginal), 2L)
  expect_equal(diagnostic$iteration$checkpoint, c("25%", "50%", "75%", "Full"))
  expect_equal(diagnostic$iteration$saved_per_chain, c(2L, 4L, 6L, 8L))
  expect_null(diagnostic$support)
})

test_that("plot_movement returns stacked speed and move/stay panels", {
  movement <- list(
    method = "gamma",
    shape = 2,
    scale = 5,
    low_speed_fix = 0.001,
    zero_speed_ratio = 0,
    move_stay_parameters = list(p_0 = 0.7, p_inf = 0.2, tau = 1),
    move_stay = function(t) 0.2 + 0.5 * exp(-t)
  )

  plot <- plot_movement(movement, speed = 0:10, stay_duration = 0:5)

  expect_s3_class(plot, "patchwork")
  expect_equal(length(plot$patches$plots), 1L)
})

test_that("14OI extdata completes the light-only sampling workflow", {
  example_dir <- system.file("extdata", package = "GeoPathSampleR")
  tag <- GeoPressureR::tag_create(
    "14OI",
    crop_start = "2015-07-17",
    crop_end = "2016-07-11",
    directory = file.path(example_dir, "data/raw-tag/14OI"),
    assert_pressure = FALSE,
    quiet = TRUE
  )
  tag <- GeoPressureR::twilight_create(tag)
  tag <- GeoPressureR::twilight_label_read(
    tag,
    file = file.path(example_dir, "data/twilight-label/14OI-labeled.csv")
  )
  tag <- GeoPressureR::tag_stap_daily(
    tag,
    stap_long = file.path(example_dir, "data/stap-label/14OI.csv"),
    movement_period = "day",
    quiet = TRUE
  )
  tag <- GeoPressureR::tag_set_map(
    tag,
    extent = c(-20, 40, -40, 60),
    scale = 1
  )
  tag <- GeoPressureR::geolight_map(
    tag,
    twl_calib_adjust = 1,
    fitted_location_duration = 30,
    twl_llp = \(n) 1.5 * log(n) / n,
    quiet = TRUE
  )

  movement <- eval(formals(sampling_path)$movement)
  paths <- sampling_path(
    tag,
    iter = 3,
    warmup = 1,
    chains = 1,
    seed = 1,
    quiet = TRUE
  )

  expect_equal(nrow(paths), 2L * nrow(tag$stap))
  expect_equal(sort(unique(paths$stap_id)), tag$stap$stap_id)
  expect_true(all(is.finite(paths$lat) & is.finite(paths$lon)))
})
