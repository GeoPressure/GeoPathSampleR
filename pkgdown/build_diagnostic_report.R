# Build the frozen pkgdown diagnostic example from the getting-started setup.

repo_dir <- normalizePath(".", mustWork = TRUE)
devtools::load_all(repo_dir, quiet = TRUE)

example_dir <- system.file("extdata", package = "GeoPathSampleR")
withr::with_dir(example_dir, {
  tag <- GeoPressureR::tag_create(
    "14OI",
    crop_start = "2015-07-17",
    crop_end = "2016-07-11",
    directory = file.path(example_dir, "data/raw-tag/14OI"),
    assert_pressure = FALSE,
    quiet = TRUE
  )
  tag <- GeoPressureR::twilight_create(tag)
  tag <- GeoPressureR::twilight_label_read(tag)
  tag <- GeoPressureR::tag_stap_daily(
    tag,
    stap_long = tag$param$id,
    movement_period = "day",
    quiet = TRUE
  )
  tag <- GeoPressureR::tag_set_map(
    tag,
    extent = c(-20, 40, -11, 60),
    scale = 1
  )
  tag <- GeoPressureR::geolight_map(
    tag,
    twl_calib_adjust = 1,
    fitted_location_duration = 30,
    twl_llp = \(n) 1.5 * log(n) / n,
    quiet = TRUE
  )
  tag$map_light <- GeoPressureR::map_add_mask_water(tag$map_light)

  paths <- GeoPathSampleR::sampling_path(
    tag,
    iter = 500,
    chains = 4,
    warmup = 100,
    seed = 1,
    quiet = TRUE
  )
  diagnostic <- GeoPathSampleR::sampling_path_diagnostic(
    paths,
    tag,
    report = FALSE,
    quiet = TRUE
  )

  GeoPathSampleR::sampling_path_diagnostic_render(
    diagnostic,
    paths,
    tag,
    browse = FALSE,
    output_file = file.path(
      repo_dir,
      "pkgdown",
      "assets",
      "sampling-path-diagnostic.html"
    )
  )
})
