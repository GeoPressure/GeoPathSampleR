#' Collapse consecutive identical locations in paths
#'
#' Summarise each realised path by merging consecutive stationary periods that
#' were assigned to the same grid cell. Paths are identified by `j` and, when
#' present, `chain`.
#'
#' @param paths A path data frame containing `stap_id`, `ind`, `lat`, and
#'   `lon`, and optionally `j`, `chain`, `start`, and `end`.
#' @param stap A stationary-period data frame containing `stap_id`, `start`,
#'   and `end`. If `NULL`, these columns are taken from `paths`.
#' @param label_stap_id logical; if `TRUE`, replace `stap_id` in the output
#'   with a character label of the form `"{stap_id_start}->{stap_id_end}`.
#'   The numeric `stap_id_start` and `stap_id_end` columns are always kept.
#'
#' @return A data frame with one row per collapsed stay and columns
#'   `stap_id`, `stap_id_start`, `stap_id_end`, `ind`, `lat`, `lon`, `start`,
#'   `end`, `known`, `j`, and `chain`.
#' @family path
#' @export
path_collapse <- function(paths, stap = NULL, label_stap_id = FALSE) {
  stopifnot(all(c("stap_id", "ind", "lat", "lon") %in% names(paths)))

  if (is.null(stap)) {
    stopifnot(all(c("start", "end") %in% names(paths)))
    stap <- unique(paths[c("stap_id", "start", "end")])
  }
  stopifnot(all(c("stap_id", "start", "end") %in% names(stap)))

  if (!"j" %in% names(paths)) {
    paths$j <- 1L
  }
  if (!"chain" %in% names(paths)) {
    paths$chain <- 1L
  }

  paths <- paths[order(paths$chain, paths$j, paths$stap_id), , drop = FALSE]
  stap_row <- match(paths$stap_id, stap$stap_id)
  known <- if (all(c("known_lat", "known_lon") %in% names(stap))) {
    is.finite(stap$known_lat[stap_row]) & is.finite(stap$known_lon[stap_row])
  } else if ("known" %in% names(paths)) {
    as.logical(paths$known)
  } else {
    rep(FALSE, nrow(paths))
  }

  path_start <- c(
    TRUE,
    paths$chain[-1L] != paths$chain[-nrow(paths)] |
      paths$j[-1L] != paths$j[-nrow(paths)]
  )
  stay_start <- path_start
  stay_start[-1L] <- path_start[-1L] |
    paths$ind[-1L] != paths$ind[-nrow(paths)]
  stay_i <- cumsum(stay_start)
  first_i <- which(stay_start)
  last_i <- c(first_i[-1L] - 1L, nrow(paths))

  out <- data.frame(
    stap_id_start = paths$stap_id[first_i],
    stap_id_end = paths$stap_id[last_i],
    ind = paths$ind[first_i],
    lat = paths$lat[first_i],
    lon = paths$lon[first_i],
    start = stap$start[stap_row[first_i]],
    end = stap$end[stap_row[last_i]],
    known = as.logical(rowsum(as.integer(known), stay_i, reorder = FALSE) > 0),
    j = paths$j[first_i],
    chain = paths$chain[first_i]
  )
  out$stap_id <- if (isTRUE(label_stap_id)) {
    paste0(out$stap_id_start, "->", out$stap_id_end)
  } else {
    out$stap_id_start
  }
  out <- out[, c(
    "stap_id",
    "stap_id_start",
    "stap_id_end",
    "ind",
    "lat",
    "lon",
    "start",
    "end",
    "known",
    "j",
    "chain"
  )]
  row.names(out) <- NULL
  out
}
