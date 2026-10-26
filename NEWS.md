# GeoPathSampleR v0.2.0

## Main

- [Add boundary sweeps to resample movement timing while retaining residence-segment locations](https://github.com/GeoPressure/GeoPathSampleR/commit/f15e91e), enabled every iteration by default with `boundary_interval = 1`. Use `0` or `Inf` to disable them.

- [Update the default ground-speed and residence-dependent departure parameters](https://github.com/GeoPressure/GeoPathSampleR/commit/d40999d) from the pooled multi-sensor calibration. Results sampled with the previous defaults should be recomputed before comparison.
- [Replace the route-detour model with a Gamma prior on route excess](https://github.com/GeoPressure/GeoPathSampleR/commit/7d4fb3f) whose scale depends on interval duration and direct endpoint distance ([3798978](https://github.com/GeoPressure/GeoPathSampleR/commit/3798978)). The prior applies to intervals with at least two movement days and endpoints at least 300 km apart, with predictors limited to the empirical calibration support.
- Export `sampling_path_default_movement()` and document the movement model, `route_detour` and the modular cut of `long_period_light_only` in `sampling_path()`.

## Minor

- Exclude the redundant repository GPL text from the R source package to resolve its package-check note.

- [Fix known locations on water cells, stabilise the prior for direct routes, and reject excluded stationary periods](https://github.com/GeoPressure/GeoPathSampleR/commit/49e7735).
- [Speed up route-prior evaluation, path summaries and diagnostics](https://github.com/GeoPressure/GeoPathSampleR/commit/6da2d15), and [avoid replicating scalar movement endpoints](https://github.com/GeoPressure/GeoPathSampleR/commit/0b989d1).
- Document boundary sweeps, the complete sampler update schedule, route-prior support and computational controls in the function reference and getting-started vignette.

- [Fix the diagnostic report rendering and browser opening](https://github.com/GeoPressure/GeoPathSampleR/commit/fb0d61f), and [declare its optional report dependencies](https://github.com/GeoPressure/GeoPathSampleR/commit/3ff7a89).
- [Align the README and vignette with the movement model and sampler controls](https://github.com/GeoPressure/GeoPathSampleR/commit/ae3f0de), and add a provisional citation for the submitted methods manuscript alongside the software citation.
- [Format the elapsed time compactly in the sampling message](https://github.com/GeoPressure/GeoPathSampleR/commit/38fe161).

**Full Changelog**: <https://github.com/GeoPressure/GeoPathSampleR/compare/v0.1.0...v0.2.0>

# GeoPathSampleR v0.1.0

* Initial public release of the stationary-path Gibbs sampler and its core
  summaries and diagnostics.
* Simplified path sampling around one movement model, sparse period
  likelihoods, batched kernel normalization, and incremental path repair.
