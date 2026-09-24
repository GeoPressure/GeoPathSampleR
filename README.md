
<!-- README.md is generated from README.Rmd. Please edit that file and use devtools::build_readme() -->

# GeoPathSampleR <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/GeoPressure/GeoPathSampleR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/GeoPressure/GeoPathSampleR/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/GeoPressure/GeoPathSampleR/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/GeoPressure/GeoPathSampleR/actions/workflows/pkgdown.yaml)
[![Codecov](https://codecov.io/gh/GeoPressure/GeoPathSampleR/graph/badge.svg)](https://app.codecov.io/gh/GeoPressure/GeoPathSampleR)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21944194.svg)](https://doi.org/10.5281/zenodo.21944194)
[![lint](https://github.com/GeoPressure/GeoPathSampleR/actions/workflows/jarl-check.yml/badge.svg)](https://github.com/GeoPressure/GeoPathSampleR/actions/workflows/jarl-check.yml)
[![format](https://github.com/GeoPressure/GeoPathSampleR/actions/workflows/format-check.yml/badge.svg)](https://github.com/GeoPressure/GeoPathSampleR/actions/workflows/format-check.yml)
<!-- badges: end -->

GeoPathSampleR reconstructs bird migration trajectories from
stationary-period likelihood maps built with GeoPressureR. Its Gibbs
sampler combines residence-dependent departure, a ground-speed
distribution conditional on movement, and a prior on cumulative route
distance between long periods. Unknown long-period locations can be
sampled from light evidence alone. The package also provides posterior
stay summaries and sampling diagnostics.

## 📦 Installation

To install the latest version from GitHub:

``` r
# install.packages("pak")
pak::pkg_install("GeoPressure/GeoPathSampleR")
```

## 📘 Vignettes

See the [getting-started
vignette](https://geopressure.org/GeoPathSampleR/articles/getting-started.html)
for inputs, movement settings, and interpretation.

## 📚 Citation

For the software, cite:

> Nussbaumer, R. (2026). *GeoPathSampleR: Bayesian reconstruction of
> animal trajectories from geolocation likelihood maps.* Zenodo.
> <https://doi.org/10.5281/zenodo.21944194>. Available at:
> <https://github.com/GeoPressure/GeoPathSampleR>

The methods and validation are described in the following submitted
manuscript (citation details will be updated after publication):

> Nussbaumer, R., Gravey, M., Wong, J. B. & Benoit, L. (2026).
> Reconstructing bird migration trajectories from light: inspectable
> likelihoods and biologically informed movement priors. Manuscript
> submitted to *Methods in Ecology and Evolution*.
