#' Plot a path-sampling movement model
#'
#' Plot the speed kernel and residence-dependent probability of movement for a
#' `sampling_path()` movement specification.
#'
#' @param movement movement list passed to `sampling_path()`.
#' @param speed speed values in km/h used to evaluate the movement kernel.
#' @param stay_duration residence durations in days used to evaluate
#'   `movement$move_stay`.
#'
#' @return A patchwork plot with speed and move/stay panels. When
#'   `movement$move_stay_parameters` contains `p_0`, `p_inf`, and `tau`, the
#'   move/stay panel marks those parameters for the exponential-decay-to-plateau
#'   model.
#' @examples
#' movement <- list(
#'   method = "gamma",
#'   shape = 1.412246,
#'   scale = 8.909248,
#'   low_speed_fix = 0.001,
#'   zero_speed_ratio = 0,
#'   move_stay_parameters = list(
#'     p_0 = 0.647893,
#'     p_inf = 0.1318134,
#'     tau = 0.8624705
#'   ),
#'   move_stay = function(t) 0.1318134 +
#'     (0.647893 - 0.1318134) * exp(-t / 0.8624705)
#' )
#' plot_movement(movement)
#' @import patchwork
#' @importFrom rlang .data
#' @family sampling_path
#' @export
plot_movement <- function(
  movement,
  speed = seq(0, 120),
  stay_duration = seq(0, 20)
) {
  sampling_path_movement_validate(movement)

  speed_data <- data.frame(
    speed = speed,
    probability = GeoPressureR::speed2prob(speed, movement)
  )
  stay_data <- data.frame(
    duration = stay_duration,
    probability = movement$move_stay(stay_duration)
  )
  move_stay_parameters <- movement$move_stay_parameters
  p_0 <- move_stay_parameters$p_0 %||% movement$move_stay(0)
  p_inf <- move_stay_parameters$p_inf %||% movement$move_stay(Inf)
  tau <- move_stay_parameters$tau
  speed_breaks <- pretty(range(speed), n = 6)
  speed_y_max <- max(speed_data$probability, na.rm = TRUE) * 1.1

  speed_plot <- ggplot2::ggplot(
    speed_data,
    ggplot2::aes(x = .data$speed, y = .data$probability)
  ) +
    ggplot2::geom_line(color = "#1A1A1A", linewidth = 1.1) +
    ggplot2::scale_y_continuous(
      name = "Probability",
      limits = c(0, speed_y_max),
      expand = ggplot2::expansion(mult = c(0, 0.03))
    ) +
    ggplot2::labs(title = "Speed kernel") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(
      name = "Speed (km/h)",
      breaks = speed_breaks,
      sec.axis = ggplot2::sec_axis(
        ~ . * 24,
        name = "Distance over 24 h (km)",
        breaks = speed_breaks * 24
      )
    )

  stay_plot <- ggplot2::ggplot(
    stay_data,
    ggplot2::aes(x = .data$duration, y = .data$probability)
  ) +
    ggplot2::geom_line(color = "#1A1A1A", linewidth = 1.1) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_point(
      data = data.frame(duration = 0, probability = p_0),
      size = 2.6
    ) +
    ggplot2::geom_hline(
      yintercept = p_inf,
      linetype = "dotted",
      linewidth = 0.35
    ) +
    ggplot2::scale_y_continuous(
      name = "Probability",
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x = "Days since arrival at current stopover",
      y = "Probability of moving next day",
      title = "Move probability"
    ) +
    ggplot2::theme_bw() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = p_0,
      label = paste0("p[0] == ", format(p_0, digits = 3)),
      parse = TRUE,
      hjust = -0.3,
      vjust = -0.8
    ) +
    ggplot2::annotate(
      "text",
      x = max(stay_duration),
      y = p_inf,
      label = paste0("p[infinity] == ", format(p_inf, digits = 3)),
      parse = TRUE,
      hjust = 1.1,
      vjust = -0.8
    )

  if (!is.null(tau)) {
    tangent_data <- data.frame(
      duration = c(0, tau),
      probability = c(p_0, p_inf)
    )
    stay_plot <- stay_plot +
      ggplot2::geom_line(
        data = tangent_data,
        ggplot2::aes(x = .data$duration, y = .data$probability),
        linetype = "dotted",
        linewidth = 0.35
      ) +
      ggplot2::annotate(
        "text",
        x = tau,
        y = p_inf,
        label = paste0("tau == ", format(tau, digits = 3), "~\"days\""),
        parse = TRUE,
        hjust = -0.1,
        vjust = -0.8
      )
  }

  speed_plot / stay_plot
}
