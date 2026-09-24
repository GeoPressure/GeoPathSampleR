# GeoPathSampleR 0.2.0 (development)

* Updated the default ground-speed and residence-dependent departure parameters
  from the pooled multi-sensor calibration. Results sampled with the previous
  defaults should be recomputed before comparison.
* Replaced the route-detour model with a Gamma prior on route excess. Its scale
  depends on interval duration and direct endpoint distance. The prior applies
  to intervals with at least two movement days and endpoints at least 300 km
  apart, with predictors limited to the empirical calibration support.
* Fixed the diagnostic report rendering and browser opening, and declared its
  optional report dependencies.
* Aligned the README and vignette with the movement model and sampler controls.
  Added a provisional citation for the submitted methods manuscript alongside
  the software citation.

# GeoPathSampleR 0.1.0

* Initial public release of the stationary-path Gibbs sampler and its core
  summaries and diagnostics.
* Simplified path sampling around one movement model, sparse period
  likelihoods, batched kernel normalization, and incremental path repair.
