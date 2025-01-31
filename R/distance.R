#' Annotate linkSet with inter/intra chromosome interactions
#'
#' @param x A linkSet object
#'
#' @return A linkSet object with an additional metadata column 'inter_type'
#' @export
#' @aliases annotateInter
#' 
#' @examples
#' data(linkExample)
#' linkExample <- annotateInter(linkExample)

setMethod("annotateInter", "linkSet", function(x) {

  if (is.null(regionsBait(x))) {
    stop("regionsBait is not annotated. Please run annotatePromoter first.")
  }

  bait_seqnames <- as.character(seqnames(regionsBait(x)))
  oe_seqnames <- as.character(seqnames(anchors(x, type="oe")))

  inter_type <- ifelse(bait_seqnames == oe_seqnames, "inter","intra")

  # Add the new column to mcols
  mcols(x)$inter_type <- inter_type

  return(x)
}
)

#' @importFrom S4Vectors mcols mcols<-
.exist_inter <- function(x){
  !is.null(mcols(x)$inter_type)
  }

#' @importFrom S4Vectors mcols mcols<-
.exist_distance <- function(x){
  !is.null(mcols(x)$distance)
  }

.get_dist_output <- function(regs, ai1, ai2, type, inter_type) {
  type <- match.arg(type, c("mid", "gap", "span"))
  chr <- as.character(seqnames(regs))

  # To get sensible distances
  swapped <- .enforce_order(ai1, ai2)
  larger <- swapped[[2]]
  smaller <- swapped[[1]]


  # Protection when all inter's.
  is.same <- inter_type=="inter"
  output <- rep(as.integer(NA), length(larger))
  if (!any(is.same)) { return(output) }

  st <- start(regs)
  en <- end(regs)
  larger <- larger[is.same]
  smaller <- smaller[is.same]
  all.as <- st[larger]
  all.ae <- en[larger]
  all.ts <- st[smaller]
  all.te <- en[smaller]

  if (type=="gap") {
    output[is.same] <- pmax(all.as, all.ts) - pmin(all.ae, all.te) - 1L
  } else if (type=="span") {
    output[is.same] <- pmax(all.ae, all.te) - pmin(all.as, all.ts) + 1L
  } else if (type=="mid") {
    output[is.same] <- as.integer(abs(all.as + all.ae - all.ts - all.te)/2L) # Need 'abs', in case later range has earlier midpoint (e.g., if nested).
  } else if (type=="diag") {
    output[is.same] <- larger - smaller
  }
  return(output)
}


#' Calculate the distance between bait and the other end
#' @aliases pairdist
#' @description
#' Outputs an integer vector specifying the distance between the interacting bins,
#' depending on the type of distance specified.
#' 
#' Example:
#' ```
#'    rangeA:  |---------|
#'    rangeB:                |---------|
#'    mid:           <----------->
#'    gap:               <-->
#'    span:    <----------------------->
#' ```
#'
#' * mid: Half the distance between the end of first range and start of second range
#' * gap: Distance between the end of first range and start of second range
#' * span: Total span from start of first range to end of second range
#'
#'
#' @param x A linkSet object
#' @param type The type of distance to calculate, either "mid", "gap", or "span"
#'
#'
#' @return A linkSet object with a new metadata column "distance"
#' @export
#'
#' @examples
#' data(linkExample)
#' linkExample <- pairdist(linkExample, type="mid")
#'
setMethod("pairdist", "linkSet", function(x, type="mid"){
    if (!is.null(regionsBait(x))) {
        ai1 <- anchor1(x)
        ai2 <- anchor2(x)
    } else {
        stop("No regionsBait found. Please annotate regionsBait first.")
    }
    if (!.exist_inter(x)){
      x <- annotateInter(x)
    }
    inter_type <- mcols(x)$inter_type
    mcols(x)$distance <- .get_dist_output(regions(x), ai1, ai2, type, inter_type)
    return(x)
})

#' Diagnose the linkSet object, return barplot of inter/intra interaction and distance distribution
#' @export
#' @aliases diagnoseLinkSet
#' @examples
#' data(linkExample)
#' diagnoseLinkSet(linkExample)
#' @param x A linkSet object
setMethod("diagnoseLinkSet", "linkSet", function(x){
  if (!.exist_inter(x)){
    x <- annotateInter(x)
  }
  if (!.exist_distance(x)){
    x <- pairdist(x)
  }
  inter_type <- mcols(x)$inter_type
  distance <- mcols(x)$distance
  inter_distance <- log10(distance[!is.na(distance)])
  df <- data.frame(distance = inter_distance)
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = distance)) +
    ggplot2::geom_histogram() +
    ggplot2::labs(title = "Distance Distribution", x = "log10 Distance(bp)", y = "Count")+ggplot2::theme_bw()
  
  # Calculate proportions manually
  type_counts <- table(inter_type)
  type_props <- prop.table(type_counts)
  df2 <- data.frame(
    inter_type = factor(names(type_props), levels = c("inter", "intra")),
    proportion = as.numeric(type_props)
  )
  
  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = inter_type, y = .data$proportion)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(title = "Inter/Intra Interaction", x = "Interaction Type", y = "Percentage") +
    ggplot2::scale_x_discrete(drop = FALSE) +  # This ensures both levels are shown
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::theme_bw()
  p <- p1|p2
  print(p)
  return(x)
})