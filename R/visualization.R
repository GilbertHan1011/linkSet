#' Add Genome Links to Coverage Plot.
#'
#' @param linkSet A linkSet object
#' @param score.col Column index that contains score information, used when \code{file.type} is bedpe. Default: NULL.
#' @param score.threshold The score threshold, used when \code{score.col} is not NULL. Default: NULL.
#' @param score.color The score color vector. Default: c("grey70", "#56B1F7", "#132B43").
#' @param scale.range Scale the height of links according to width, should be greater than or equal to 1 (not scale). Default: 10.
#' @param plot.space Top and bottom margin. Default: 0.1.
#' @param plot.height The relative height of link to coverage plot. Default: 0.2.
#' @param show.rect Logical value, whether to add rect border to the plot. Default: FALSE.
#' @param extend.base Extend the base pair range to show more information. Default: 1000000.
#' @param x.range The range of x-axis to show. Default: NULL.
#' @param log.scale Logical value, whether to log1p the score. Default: TRUE.
#' @param arrow.size The size of the arrow head. Default: 0.02.
#' @param remove_x_axis Logical value, whether to remove the x-axis. Default: FALSE.
#' @param link_plot_on_top Logical value, whether to plot the link plot on top of the coverage plot. Default: FALSE.
#'
#' @return Plot.
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame start end
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom utils read.table
#' @importFrom scales rescale
#' @importFrom ggplot2 ggplot_add ggplot aes_string scale_color_gradientn
#'   labs theme_classic theme element_blank element_rect
#'   element_text margin scale_y_continuous scale_x_continuous expansion
#'   coord_cartesian geom_curve
#' @importFrom patchwork wrap_plots
#' @references \url{https://stuartlab.org/signac/articles/cicero.html}
#' @export
#'
#'
setMethod("geom_linkset", "linkSet", function(linkSet,
                      score.col = "count",
                      score.threshold = NULL,
                      score.color = c("grey70", "#56B1F7" , "#132B43"),
                      scale.range = 10,
                      plot.space = 0.1,
                      plot.height = 0.2,
                      arrow.size = 0.05,
                      remove_x_axis = FALSE,
                      link_plot_on_top = FALSE,
                      extend.base = 1000000,
                      show.rect = FALSE,
                      x.range = NULL,
                      log.scale = TRUE) {
  structure(
    list(
      linkSet = linkSet,
      score.col = score.col,
      score.threshold = score.threshold,
      score.color = score.color,
      scale.range = scale.range,
      plot.space = plot.space,
      plot.height = plot.height,
      show.rect = show.rect,
      arrow.size = arrow.size,
      remove_x_axis = remove_x_axis,
      link_plot_on_top = link_plot_on_top,
      extend.base = extend.base,
      x.range = x.range,
      log.scale = log.scale
    ),
    class = "interSet"
  )
})



#' @export
ggplot_add.interSet <- function(object, plot, object_name) {
  # get plot data
  if ("patchwork" %in% class(plot)) {
    track.data <- plot[[1]]$layers[[1]]$data
  } else {
    track.data <- plot$layers[[1]]$data
  }

  # get parameters
  linkSet <- object$linkSet
  score.col <- object$score.col
  score.threshold <- object$score.threshold
  score.color <- object$score.color
  scale.range <- object$scale.range
  plot.curve <- object$plot.curve
  plot.space <- object$plot.space
  plot.height <- object$plot.height
  show.rect <- object$show.rect
  arrow.size <- object$arrow.size
  remove_x_axis <- object$remove_x_axis
  link_plot_on_top <- object$link_plot_on_top
  flip_arrow <- link_plot_on_top
  top_margin <- bottom_margin <- plot.space
  extend.base <- object$extend.base
  x.range <- object$x.range
  log.scale <- object$log.scale


  # prepare plot range
  plot.range.chr <- as.character(seqnames(regionsBait(object$linkSet))[1])
  plot.range.start <- min(start(regions(object$linkSet))) - extend.base
  plot.range.end <- max(end(regions(object$linkSet))) + extend.base
  if (is.null(x.range)) {
    x.range <- c(plot.range.start, plot.range.end)
  }
  # prepare dataframe
  link.point.df <- data.frame(
    chr = as.character(seqnames(regionsBait(linkSet))),
    start = start(oe(linkSet)),
    end = start(regionsBait(linkSet))
  )

  # add score
  if (score.col %in% colnames(mcols(linkSet))) {
    link.point.df$score <- mcols(linkSet)[[score.col]]
    if (!is.null(score.threshold)) {
      link.point.df <- link.point.df[link.point.df$score > score.threshold, ]
    }
    if (log.scale) {
      link.point.df$score <- log1p(link.point.df$score)
    }
  }
  # filter link gr
  link.point.df <- link.point.df[link.point.df$start >= x.range[1] &
                                link.point.df$end >= x.range[1] &
                                link.point.df$end <= x.range[2]&
                                link.point.df$start <= x.range[2], ]
  if(nrow(link.point.df) == 0){
    warning("There are no valid links in the given region!")
    return(NULL)
  }
  rownames(link.point.df) <- 1:nrow(link.point.df)
  # check dataframe
  if (nrow(link.point.df) < 1) {
    warning("There are no valid links in the given region!")
    # create empty plot
    link.basic.plot <- ggplot2::ggplot(data = link.point.df)
  } else {
    # prepare plot dataframe
    link.point.df$group <- seq_len(length.out = nrow(link.point.df))
    link.point.plot <- link.point.df
    link.point.plot$width <- link.point.df$end - link.point.df$start
    #browser()
    # scale width to range
    link.point.plot$rw <- scales::rescale(link.point.plot$width, to = c(1, scale.range))

    if ("score" %in% colnames(link.point.plot)) {
      group_color <- "score"
      scale_color <- ggplot2::scale_color_gradientn(
        colors = score.color,
        limits = c(0, max(link.point.plot$score))
      )
    } else {
      group_color <- NULL
      scale_color <- ggplot2::scale_color_manual()
    }

    y_limit <- ifelse(flip_arrow, 0, 1)
    link.point.plot.pos <- link.point.plot[link.point.plot$width > 0,]
    link.point.plot.neg <- link.point.plot[link.point.plot$width < 0,]
    link.basic.plot <-
      ggplot2::ggplot(data = link.point.plot) +
      ggplot2::geom_curve(
        data = link.point.plot.pos,
        ggplot2::aes_string(
          x = "start",
          xend = "end",
          y = y_limit,
          yend = y_limit,
          color = group_color,
          size = "score"
        ),
        curvature = ifelse(flip_arrow, -0.2, 0.2),
        angle = 90,
        ncp = 15,
        arrow = ggplot2::arrow(length = ggplot2::unit(arrow.size, "npc"))
      ) +
      ggplot2::geom_curve(
        data = link.point.plot.neg,
        ggplot2::aes_string(
          x = "start",
          xend = "end",
          y = y_limit,
          yend = y_limit,
          color = group_color,
          size = "score"
        ),
        curvature = ifelse(flip_arrow, 0.2, -0.2),
        angle = 90,
        ncp = 15,
        arrow = ggplot2::arrow(length = ggplot2::unit(arrow.size, "npc"))
      ) +
      scale_color +
      ggplot2::scale_y_continuous(limits = c(0,1)) +
      ggplot2::scale_size_continuous(range = c(0.5, 2))
  }

  # create plot
  link.plot <-
    link.basic.plot +
    ggplot2::labs(y = "Links") +
    theme_linkset(
      x.range = x.range,
      margin.len = plot.space,
      show.rect = show.rect
    ) +
    ggplot2::guides(size = "none")  # Remove legend for arrow size

  # Add chromosome name to the side
  link.plot <- link.plot +
    ggplot2::annotate("text", x = min(x.range), y = 1,
                      label = plot.range.chr, hjust = 0, vjust = 1)

  # assemble plot
  patchwork::wrap_plots(
    plot + ggplot2::theme(plot.margin = ggplot2::margin(t = plot.space, b = plot.space)),
    link.plot,
    ncol = 1,
    heights = c(1, plot.height)
  )
  # Create a function to adjust plot margins and remove x-axis elements
  adjust_plot <- function(p, top_margin, bottom_margin, remove_x_axis = FALSE) {
    p <- p + ggplot2::theme(plot.margin = ggplot2::margin(t = top_margin, b = bottom_margin))
    if (remove_x_axis) {
      p <- p + ggplot2::theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    }
    return(p)
  }

  # Adjust plots based on user preferences
  adjusted_link_plot <- adjust_plot(link.plot, top_margin = top_margin, bottom_margin = bottom_margin)
  adjusted_plot <- adjust_plot(plot, top_margin = top_margin, bottom_margin = bottom_margin, remove_x_axis = object$remove_x_axis)

  # Determine plot order and create list of plots
  plot_list <- if (object$link_plot_on_top) {
    list(adjusted_link_plot, adjusted_plot)
  } else {
    list(adjusted_plot, adjusted_link_plot)
  }

  # Calculate heights based on plot order and overlap
  total_height <- 1 + object$plot.height
  heights <- if (object$link_plot_on_top) {
    c(object$plot.height, total_height - object$plot.height)
  } else {
    c(total_height - object$plot.height, object$plot.height)
  }

  # Combine plots
  combined_plot <- patchwork::wrap_plots(
    plot_list,
    ncol = 1,
    heights = heights
  ) +
    patchwork::plot_layout(guides = "collect")
  return(combined_plot)
}

#' Plot genomic ranges
#'
#' `geom_range()` and `geom_half_range()` draw tiles that are designed to
#' represent range-based genomic features, such as exons. In combination with
#' `geom_intron()`, these geoms form the core components for visualizing
#' transcript structures.
#'
#' `geom_range()` and `geom_half_range()` require the following `aes()`;
#' `xstart`, `xend` and `y` (e.g. transcript name). `geom_half_range()` takes
#' advantage of the vertical symmetry of transcript annotation by plotting only
#' half of a range on the top or bottom of a transcript structure. This can be
#' useful for comparing between two transcripts or free up plotting space for
#' other transcript annotations (e.g. `geom_junction()`).
#'
#' @param bait_col Color for bait regions. Default is "red".
#' @param oe_col Color for other end regions. Default is "DeepSkyBlue3".
#' @param default_col Default color for regions. Default is "grey".
#' @param minimal_width Minimal width for the range. Default is 0.01.
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#' @inheritParams ggplot2::geom_tile
#' @inheritParams ggplot2::geom_segment
#' @inheritParams grid::rectGrob
#' 
#' @importFrom rlang %||%
#'
#' @return the return value of a `geom_*` function is not intended to be
#'   directly handled by users. Therefore, `geom_*` functions should never be
#'   executed in isolation, rather used in combination with a
#'   `ggplot2::ggplot()` call.
#'
#' @export

geom_range <- function(mapping = NULL, data = NULL,
                       stat = "identity", position = "identity",
                       bait_col = "red",
                       oe_col = "DeepSkyBlue3",
                       default_col = "grey",
                       ...,
                       vjust = NULL,
                       linejoin = "mitre",
                       na.rm = FALSE,
                       minimal_width = 0.01,
                       show.legend = NA,
                       inherit.aes = TRUE) {
    ggplot2::layer(
        data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomRange,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            minimal_width = minimal_width,
            vjust = vjust,
            linejoin = linejoin,
            na.rm = na.rm,
            bait_col = bait_col,
            oe_col = oe_col,
            default_col = default_col,
            ...
        )
    )
}

#' `GeomRange` is `ggplot2::GeomTile` with modified `aes` to match genetic
#' nomenclature (`xstart`/`xend`)
#' @keywords internal
#' @noRd
GeomRange <- ggplot2::ggproto("GeomRange", ggplot2::GeomTile,
    required_aes = c("xstart", "xend","region"),
    default_aes = ggplot2::aes(
        fill = "grey",
        colour = "black",
        linewidth = 0.25,
        linetype = 1,
        alpha = NA,
        height = NA
    ),
    setup_data = function(data, params) {
        # modified from ggplot2::GeomTile

        data$height <- data$height %||% params$height %||% 0.5
        data$fill <- with(data, dplyr::case_when(
            region == "bait" ~ params$bait_col,
            region == "oe" ~ params$oe_col,
            TRUE ~ params$default_col  # Default color if not bait or oe
        ))

        transform(
            data,
            xmin = xstart,
            xmax = xend,
            ymin =  - data$height / 200,
            ymax =  data$height / 200,
            height = NULL
        )
    },
    draw_panel = function(self,
                          data,
                          panel_params,
                          coord,
                          vjust = NULL,
                          minimal_width = 0.01,
                          bait_col = "red",
                          oe_col = "DeepSkyBlue3",
                          default_col = "grey",
                          lineend = "butt",
                          linejoin = "mitre") {
        if (!coord$is_linear()) {
            # prefer to match geom_curve and warn
            # rather than copy the implementation from GeomRect for simplicity
            # also don'think geom_range would be used for non-linear coords
            warn("geom_ is not implemented for non-linear coordinates")
        }

        coords <- coord$transform(data, panel_params)
        grid::rectGrob(
            coords$xmin, coords$ymax,
            width = max(coords$xmax - coords$xmin, minimal_width),
            height = coords$ymax - coords$ymin,
            default.units = "native",
            just = c("left", "top"),
            vjust = vjust,
            gp = grid::gpar(
                col = coords$colour,
                fill = ggplot2::alpha(coords$fill, coords$alpha),
                lwd = coords$linewidth * ggplot2::.pt,
                lty = coords$linetype,
                linejoin = linejoin,
                lineend = lineend
            )
        )
    }
)

#' plot genomic ranges and links
#' @export
#' @aliases plot_genomic_ranges,linkSet-method
#' @aliases plot_genomic_ranges
#' @param linkset A `linkSet` object.
#' @param showBait A character vector specifying the bait region to be shown. Default: NULL.
#' @param showOE A `GRanges` object specifying the oe region to be shown. Default: NULL.
#' @param x.range A numeric vector of length 2 specifying the x-axis range. Default: NULL.
#' @param score.col A character string specifying the column name of the score. Default: "count".
#' @param show.rect Logical value, whether to show the rectangle. Default: TRUE.
#' @param extend.base A numeric value specifying the extension base. Default: 1000000.
#' @param ... Additional arguments.
#' @param bait_col A character string specifying the color of the bait region. Default: "red".
#' @param oe_col A character string specifying the color of the oe region. Default: "DeepSkyBlue3".
#' @param default_col A character string specifying the color of the default region. Default: "grey".
#' @param vjust A numeric value specifying the vertical justification. Default: NULL.
#' @param linejoin A character string specifying the line join. Default: "mitre".
#' @param na.rm Logical value, whether to remove NA values. Default: FALSE.
#' @param minimal_width A numeric value specifying the minimal width. Default: 0.01.
#' @param show.legend Logical value, whether to show the legend. Default: NA.
#' @param inherit.aes Logical value, whether to inherit the aesthetics. Default: TRUE.
#' @param link_plot_on_top Logical value, whether to plot the link plot on top of the coverage plot. Default: FALSE.
#' @param arrow.size A numeric value specifying the size of the arrow head. Default: 0.05.
#' @param remove_x_axis Logical value, whether to remove the x-axis. Default: FALSE.
#' @param plot.height A numeric value specifying the height of the plot. Default: 0.4.
#' @param plot.space A numeric value specifying the space between the plot and the link plot. Default: 0.1.
#' @param log.scale Logical value, whether to log scale the score. Default: TRUE.
#' @return A `ggplot` object.
#' @examples
#' data(linkExample)
#' plot_genomic_ranges(linkExample, extend.base = 10)
#' @export
#' 
setMethod("plot_genomic_ranges", "linkSet", function(linkset, showBait = NULL,
                                showOE = NULL, x.range = NULL,
                                score.col = "count",
                                show.rect = TRUE, extend.base = 1000000,
                                ...,
                                bait_col = "red",
                                oe_col = "DeepSkyBlue3",
                                default_col = "grey",
                                vjust = NULL,
                                linejoin = "mitre",
                                na.rm = FALSE,
                                minimal_width = 0.01,
                                show.legend = NA,
                                inherit.aes = TRUE,
                                link_plot_on_top = FALSE,
                                arrow.size = 0.05,
                                remove_x_axis = FALSE,
                                plot.height = 0.4,
                                plot.space = 0.1,
                                log.scale = TRUE) {

    # Subset linkset if bait or oe is provided
    if (!is.null(showBait)) {
        if (!is.character(showBait)) {
            stop("bait should be a character vector")
        }
        linkset <- subsetBait(linkset, showBait)
    }

    if (!is.null(showOE)) {
        if (!is(showOE, "GRanges")) {
            stop("oe should be a GRanges object")
        }
        linkset <- subsetOE(linkset, showOE)
    }

    # Extract data from linkset object
    if (is.null(x.range)) {
        plot.range.start <- min(start(regions(linkset))) - extend.base
        plot.range.end <- max(end(regions(linkset))) + extend.base
        x.range <- c(plot.range.start, plot.range.end)
    }

    if (is.null(regionsBait(linkset))) {
        stop("No bait region found, please annotate the bait region first")
    }

    if (!score.col %in% colnames(mcols(linkset))) {
      warning("score.col not found, using count as default")
      score.col <- "count"
      linkset <- countInteractions(linkset)
    }

    data <- extract_data_from_linkset(linkset)

    # Create the base plot
    p <- ggplot2::ggplot(data = data) +
        geom_range(
            ggplot2::aes(xstart = xstart, xend = xend, region = region),
            minimal_width = minimal_width,
            bait_col = bait_col,
            oe_col = oe_col,
            default_col = default_col,
            vjust = vjust,
            linejoin = linejoin,
            na.rm = na.rm,
            show.legend = show.legend,
            inherit.aes = inherit.aes,
            ...
        )

    # Apply the theme_linkset
    p <- p + ggplot2::labs(y = "Ranges") + theme_range(x.range, show.rect)
    p  <- p + geom_linkset(linkset, score.col = score.col, x.range = x.range, link_plot_on_top = link_plot_on_top,
                            show.rect = show.rect, extend.base = extend.base,
                            arrow.size = arrow.size, remove_x_axis = remove_x_axis,
                            plot.height = plot.height, plot.space = plot.space, log.scale = log.scale)
    return(p)
})

#' Extract data from linkSet object
#' @keywords internal
#' @noRd
extract_data_from_linkset <- function(linkset) {
    # Extract data from linkset object

    region_bait <- as.data.frame(regionsBait(linkset))
    region_oe <- as.data.frame(oe(linkset))

    region_bait$region <- "bait"
    region_oe$region <- "oe"
    regionDf <- rbind(region_bait,region_oe)
    regionDf <- unique(regionDf)
    regionDf <- regionDf[,c("start","end","region")]
    colnames(regionDf) <- c("xstart","xend","region")

    return(regionDf)
}


#' linkSet-theme
#' @export
#' @rdname linkSet-theme
#' @param x.range A numeric vector specifying the x-axis range.
#' @param margin.len A numeric value specifying the margin length.
#' @param show.rect A logical value indicating whether to show the rectangle.
#' @examples
#' data(linkExample)
#' x.range <- c(0, 100)
#' margin.len <- 0.5
#' show.rect <- TRUE
#' p <- ggplot2::ggplot() + theme_linkset(x.range, margin.len, show.rect)
#' print(p)
theme_linkset <- function(x.range, margin.len, show.rect) {
  if (show.rect) {
    list(
      ggplot2::theme_classic(),
      ggplot2::theme(
        axis.line.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        plot.margin = margin(t = margin.len, b = margin.len)
      ),
      #ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(0.5, 0)), position = "right"),
      ggplot2::scale_x_continuous(expand = c(0, 0)),
      ggplot2::coord_cartesian(xlim = x.range)
    )
  } else {
    list(
      ggplot2::theme_classic(),
      ggplot2::theme(
        axis.line.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y.right = ggplot2::element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = margin.len, b = margin.len)
      ),
      ggplot2::scale_x_continuous(expand = c(0, 0)),
      ggplot2::coord_cartesian(xlim = x.range)
    )
  }
}

#' @export
#' @rdname linkSet-theme
#' @examples
#' data(linkExample)
#' x.range <- c(0, 100)
#' show.rect <- TRUE
#' p <- ggplot2::ggplot() + theme_range(x.range, show.rect)
#' print(p)
#' 
theme_range <- function(x.range, show.rect) {
  if (show.rect) {
    list(
      ggplot2::theme_classic(),
      ggplot2::theme(
        axis.line.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        #plot.margin = margin(t = margin.len, b = margin.len)
      ),
      #ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(0.5, 0)), position = "right"),
      ggplot2::scale_x_continuous(expand = c(0, 0)),
      ggplot2::coord_cartesian(xlim = x.range)
    )
  } else {
    list(
      ggplot2::theme_classic(),
      ggplot2::theme(
        axis.line.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y.right = ggplot2::element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        #plot.margin = ggplot2::margin(t = margin.len, b = margin.len)
      ),
      ggplot2::scale_x_continuous(expand = c(0, 0)),
      ggplot2::coord_cartesian(xlim = x.range)
    )
  }
}

#' Plot baits in a linkSet object
#' @title Plot Baits
#' @importFrom rlang .data
#' @importFrom GenomicRanges strand
#' @name plotBaits
#' @param linkset A linkSet object
#' @param ... Additional plotting parameters
#' @return A ggplot object
#' @export
plotBaits <- function(linkset, scoreCol = "score", countCol = "count", n = 4, baits = NULL, plotBaitNames = TRUE, 
                      plevel1 = 5, plevel2 = 3, outfile = NULL,
                      width = 20, height = 20, extend.base = 1e6, bgCol = "black", lev2Col = "blue", 
                      lev1Col = "red", ...) {

  if (!is(linkset, "linkSet")) {
    stop("Input must be a linkSet object")
  }
  if (is.null(baits)) {
    baits <- sample(bait(linkset), n)
  } else {
    n <- length(baits)
  }
  
  # Pre-compute color vector
  color_vector <- c(bgCol, lev2Col, lev1Col)
  
  plot_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    bait <- baits[i]
    this <- subsetBait(linkset, bait)
    
    if (!is.null(extend.base)) {
      baitGr <- regionsBait(this)
      baitGr<- unique(baitGr)
      new_start <- max(0, start(baitGr) - extend.base)
      new_end <- end(baitGr) + extend.base
      
      expandGr <- GRanges(
        seqnames = seqnames(baitGr),
        ranges = IRanges(start = new_start, end = new_end),
        strand = strand(baitGr)
      )
      this <- subsetOE(this, expandGr)
    }
    if (length(this) == 0) {
      warning(paste0("No interactions found for bait ", bait))
      next
    }

    plotDf <- as.data.frame(this)
    plotDf$oe_middle <- start(oe(this)) + (end(oe(this)) - start(oe(this))) / 2
    bait_middle <- start(baitGr) + (end(baitGr) - start(baitGr)) / 2
    
    # Compute color factor
    plotDf$color_factor <- cut(plotDf[[scoreCol]], 
                               breaks = c(-Inf, plevel2, plevel1, Inf), 
                               labels = c(1, 2, 3))
    
    title <- if (plotBaitNames) {
      baitName <- bait
      if (grepl(",", baitName)) {
        baitName <- sub("(\\S+,).+", "\\1...", baitName)
      }
      paste0(baitName, " (", as.character(bait), ")")
    } else {
      as.character(bait)
    }
    p <- ggplot2::ggplot(plotDf, 
                         ggplot2::aes(x = .data$oe_middle, y = .data[[countCol]], color = .data$color_factor)) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(xintercept = bait_middle, color = "grey", linetype = "dashed") +
      ggplot2::labs(title = title, x = "Distance from viewpoint", y = countCol) +
      ggplot2::theme_minimal() +
      ggplot2::scale_y_continuous(limits = c(0, NA)) +
      ggplot2::scale_color_manual(values = color_vector, guide = "none")
    
    plot_list[[i]] <- p
  }

  combined_plot <- patchwork::wrap_plots(plotlist = plot_list, ncol = n)

  if (!is.null(outfile)) {
    ggplot2::ggsave(outfile, combined_plot, width = width, height = height)
  } else {
    print(combined_plot)
  }
  invisible(baits)
}