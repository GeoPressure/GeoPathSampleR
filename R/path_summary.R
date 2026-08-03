#' Summarise an ensemble of paths
#'
#' Compute either one posterior location per stationary period or one path
#' using the consensus stay/move structure across the ensemble. The function
#' accepts both [sampling_path()] output and simulation paths with `j` but no
#' `chain` column.
#'
#' @param paths A path data frame containing `stap_id`, `ind`, `lat`, and
#'   `lon`, and optionally `j`, `chain`, `start`, and `end`.
#' @param stap A stationary-period data frame containing `stap_id`, `start`,
#'   and `end`. If `NULL`, these columns are taken from `paths`.
#' @param by Summary resolution: `"stap"` returns one row per stationary
#'   period; `"consensus_stay"` returns one row per consensus stay.
#' @param position Position summary: `"mean"`, `"median"`, or `"mode"`.
#'   The mode selects the most frequent grid cell and therefore retains a
#'   valid `ind`; the mean and median return `NA` for `ind`.
#' @param move_threshold Probability above which a boundary is classified as a
#'   move when `by = "consensus_stay"`.
#'
#' @return A data frame with one row per requested summary unit and columns
#'   `stap_id`, `stap_id_start`, `stap_id_end`, `ind`, `lat`, `lon`, `start`,
#'   `end`, and `n_path`. Consensus summaries also contain movement
#'   probabilities at the start and end boundaries.
#' @family path
#' @export
path_summary <- function(
  paths,
  stap = NULL,
  by = c("stap", "consensus_stay"),
  position = c("mean", "median", "mode"),
  move_threshold = 0.5
) {
  stopifnot(all(c("stap_id", "ind", "lat", "lon") %in% names(paths)))
  by <- match.arg(by)
  position <- match.arg(position)

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
  path_id <- interaction(paths$chain, paths$j, drop = TRUE)
  stap <- stap[order(stap$stap_id), , drop = FALSE]
  stap_id <- stap$stap_id

  position_summary <- function(df) {
    df <- df[is.finite(df$lat) & is.finite(df$lon), , drop = FALSE]
    if (position == "mode") {
      df <- df[!is.na(df$ind), , drop = FALSE]
      if (!nrow(df)) {
        return(c(ind = NA_real_, lat = NA_real_, lon = NA_real_))
      }
      cell <- names(sort(table(df$ind), decreasing = TRUE))[1L]
      row <- df[which(df$ind == as.numeric(cell))[1L], , drop = FALSE]
      return(c(
        ind = as.numeric(cell),
        lat = row$lat,
        lon = row$lon
      ))
    }
    statistic <- if (position == "mean") mean else stats::median
    c(
      ind = NA_real_,
      lat = statistic(df$lat, na.rm = TRUE),
      lon = statistic(df$lon, na.rm = TRUE)
    )
  }

  if (by == "stap") {
    out <- do.call(
      rbind,
      lapply(stap_id, function(id) {
        df <- paths[paths$stap_id == id, , drop = FALSE]
        pos <- position_summary(df)
        data.frame(
          stap_id = id,
          stap_id_start = id,
          stap_id_end = id,
          ind = as.integer(pos[["ind"]]),
          lat = pos[["lat"]],
          lon = pos[["lon"]],
          start = stap$start[match(id, stap$stap_id)],
          end = stap$end[match(id, stap$stap_id)],
          n_path = length(unique(path_id[paths$stap_id == id]))
        )
      })
    )
    row.names(out) <- NULL
    return(out)
  }

  paths <- paths[order(path_id, paths$stap_id), , drop = FALSE]
  path_id <- interaction(paths$chain, paths$j, drop = TRUE)
  n_stap <- length(stap_id)
  ind <- matrix(paths$ind, nrow = n_stap)
  boundary_probability <- c(
    NA_real_,
    rowMeans(ind[-1L, , drop = FALSE] != ind[-n_stap, , drop = FALSE])
  )
  group <- cumsum(c(1L, boundary_probability[-1L] > move_threshold))
  group_start <- match(seq_len(max(group)), group)
  group_end <- c(group_start[-1L] - 1L, n_stap)

  out <- do.call(
    rbind,
    lapply(seq_along(group_start), function(i) {
      ids <- stap_id[group_start[i]:group_end[i]]
      df <- paths[paths$stap_id %in% ids, , drop = FALSE]
      pos <- position_summary(df)
      data.frame(
        stap_id = ids[1L],
        stap_id_start = ids[1L],
        stap_id_end = ids[length(ids)],
        ind = as.integer(pos[["ind"]]),
        lat = pos[["lat"]],
        lon = pos[["lon"]],
        start = stap$start[match(ids[1L], stap$stap_id)],
        end = stap$end[match(ids[length(ids)], stap$stap_id)],
        n_path = length(unique(path_id)),
        move_probability_start = boundary_probability[group_start[i]],
        move_probability_end = if (group_end[i] < n_stap) {
          boundary_probability[group_end[i] + 1L]
        } else {
          NA_real_
        }
      )
    })
  )
  row.names(out) <- NULL
  out
}
