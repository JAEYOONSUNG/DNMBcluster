#' Flower plot — per-genome core / unique / accessory / absent counts
#'
#' DNMB-native rewrite of the BPGAconverter flower_plot. Reads directly
#' from the DNMB Parquet outputs via [compute_per_genome_stats()] — no
#' stats.xls file required.
#'
#' @param dnmb Output of [load_dnmb()].
#' @param output_file Optional output PDF path. When NULL no file is written.
#' @param label_outer If TRUE align organism labels on the outside of
#'   the petals (default TRUE).
#' @return A ggplot object.
#' @export
flower_plot <- function(dnmb, output_file = NULL, label_outer = TRUE) {
  stats <- compute_per_genome_stats(dnmb)

  if (nrow(stats) < 2) {
    stop("flower_plot needs at least 2 genomes; got ", nrow(stats))
  }

  # Build the legacy BPGA data frame shape for the plotting code below.
  data <- data.frame(
    `Organism name`               = stats$organism %>% dplyr::coalesce(stats$genome_key),
    `No. of core genes`           = stats$core,
    `No. of unique genes`         = stats$unique,
    `No. of accessory genes`      = stats$accessory,
    `No. of exclusively absent genes` = stats$absent,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  categories      <- data$`Organism name`
  core_genes      <- data$`No. of core genes`
  unique_genes    <- data$`No. of unique genes`
  accessory_genes <- data$`No. of accessory genes`
  absent_genes    <- data$`No. of exclusively absent genes`

  num_petals   <- length(categories)
  petal_width  <- num_petals / 4
  petal_height <- petal_width * 5

  flower_data <- data.frame(
    category  = categories,
    unique    = unique_genes,
    accessory = accessory_genes,
    absent    = absent_genes,
    angle     = seq(0, 360, length.out = num_petals + 1)[-1]
  )

  create_petal <- function(center_x, center_y, width, height, angle, label) {
    t <- seq(0, pi, length.out = 100)
    x <- width  * cos(t)
    y <- height * sin(t)
    rotated_x <- x * cos(angle) - y * sin(angle) + center_x
    rotated_y <- x * sin(angle) + y * cos(angle) + center_y
    data.frame(x = rotated_x, y = rotated_y, label = label)
  }

  petals <- do.call(rbind, lapply(seq_len(nrow(flower_data)), function(i) {
    angle_rad <- flower_data$angle[i] * pi / 180
    shift_x <- -petal_height / 8 * sin(angle_rad)
    shift_y <-  petal_height / 8 * cos(angle_rad)
    create_petal(
      center_x = shift_x,
      center_y = shift_y,
      width    = petal_width,
      height   = petal_height,
      angle    = angle_rad,
      label    = flower_data$category[i]
    )
  }))

  core_circle_radius      <- petal_width * 1.2
  accessory_circle_radius <- core_circle_radius * 2.0
  core_circle_data <- data.frame(
    x = core_circle_radius * cos(seq(0, 2 * pi, length.out = 100)),
    y = core_circle_radius * sin(seq(0, 2 * pi, length.out = 100))
  )
  accessory_circle_data <- data.frame(
    x = accessory_circle_radius * cos(seq(0, 2 * pi, length.out = 100)),
    y = accessory_circle_radius * sin(seq(0, 2 * pi, length.out = 100))
  )

  flower_organism <- flower_data %>%
    dplyr::mutate(
      text_x     = (petal_height * 1.15) * -sin(angle * pi / 180),
      text_y     = (petal_height * 1.15) *  cos(angle * pi / 180),
      text_angle = ifelse(angle <= 180, angle - 90, angle + 90),
      hjust      = if (label_outer) {
        ifelse(angle <= 180, 1, 0)
      } else {
        ifelse(angle <= 180, 0, 1)
      }
    )

  flower_unique <- flower_data %>%
    dplyr::mutate(
      text_x     = (petal_height * 0.9) * -sin(angle * pi / 180),
      text_y     = (petal_height * 0.9) *  cos(angle * pi / 180),
      text_angle = ifelse(angle <= 180, angle - 90, angle + 90)
    )

  flower_accessory <- flower_data %>%
    dplyr::mutate(
      text_x     = (petal_height * 0.36) * -sin(angle * pi / 180),
      text_y     = (petal_height * 0.36) *  cos(angle * pi / 180),
      text_angle = ifelse(angle <= 180, angle - 90, angle + 90)
    )

  plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = petals,
      ggplot2::aes(x = x, y = y, fill = label, group = label),
      color = "white", alpha = 0.25
    ) +
    ggplot2::geom_polygon(
      data = accessory_circle_data,
      ggplot2::aes(x = x, y = y),
      fill = "yellow", color = "white", alpha = 0.5
    ) +
    ggplot2::geom_polygon(
      data = core_circle_data,
      ggplot2::aes(x = x, y = y),
      fill = "white", color = "black"
    ) +
    ggplot2::geom_text(
      data = flower_organism,
      ggplot2::aes(x = text_x, y = text_y, label = category, angle = text_angle, hjust = hjust),
      size = 5
    ) +
    ggplot2::geom_text(
      data = flower_unique,
      ggplot2::aes(
        x = text_x, y = text_y,
        label = paste(
          unique,
          paste0("(", formatC(as.numeric(absent), format = "f", digits = 0, big.mark = ","), ")"),
          sep = " "
        ),
        angle = text_angle
      ),
      size = 5, fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = flower_accessory,
      ggplot2::aes(
        x = text_x, y = text_y,
        label = formatC(as.numeric(accessory), format = "f", digits = 0, big.mark = ","),
        angle = text_angle
      ),
      size = 5, fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = data.frame(x = 0, y = 0),
      ggplot2::aes(
        x = x, y = y,
        label = paste(
          "Core gene",
          formatC(sum(as.numeric(core_genes)) / num_petals, format = "f", digits = 0, big.mark = ","),
          sep = "\n"
        )
      ),
      size = 6, fontface = "bold"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, plot, width = 10, height = 10)
    message("flower_plot written to: ", output_file)
  }

  invisible(plot)
}
