#' @rdname countInteractions
#' @export
setMethod("countInteractions", "linkSet", function(x, baitRegions = TRUE) {
  # Get the regions and bait names
  reg <- regions(x)
  bait_names <- bait(x)
  bait_regions <- paste(regionsBait(x))
  oe_regions <- paste(oe(x))

  if (baitRegions) {
    interaction_ids <- paste(bait_names, bait_regions, oe_regions, sep = "___")
  } else {
    interaction_ids <- paste(bait_names, oe_regions, sep = "___")
  }

  # Count occurrences of each unique interaction
  interaction_counts <- table(interaction_ids)

  # Identify unique interactions
  unique_indices <- !duplicated(interaction_ids)

  # Create new linkSet with unique interactions
  unique_linkSet <- x[unique_indices]
  unique_linkSet <- clean_unused_regions(unique_linkSet)
  # Add interaction counts as metadata
  indiceName <- interaction_ids[unique_indices]
  mcols(unique_linkSet)$count <- as.vector(interaction_counts[indiceName])

  # Return the results
  return(unique_linkSet)
})


#' @rdname countInteractibility
#' @export
setMethod("countInteractibility", "linkSet", function(x, baitRegions = TRUE) {
  # Ensure inter_type and count are present
  if (!"inter_type" %in% colnames(mcols(x))) {
    x <- annotateInter(x)
  }
  if (!"count" %in% colnames(mcols(x))) {
    x <- countInteractions(x)
  }

  # Get bait and oe identifiers
  bait_names <- bait(x)
  bait_regions <- paste(regionsBait(x))
  oe_regions <- paste(oe(x))
  trans <- x[x$inter_type == "intra"]

  if(length(trans) == 0){
    warning("No intra-chromosomal interactions found. Please run this function before you filterLinks.")
    return(x)
  }

  bait_names_trans <- bait(trans)
  bait_regions_trans <- paste(regionsBait(trans))
  oe_regions_trans <- paste(oe(trans))

  if (baitRegions) {
    bait_ids <- paste(bait_names, bait_regions, sep = "___")
    bait_ids_trans <- paste(bait_names_trans, bait_regions_trans, sep = "___")
  } else {
    bait_ids <- bait_names
    bait_ids_trans <- bait_names_trans
  }

  oe_ids <- oe_regions
  oe_ids_trans <- oe_regions_trans
  # Summarize counts
  bait_counts <- tapply(trans$count, bait_ids_trans, sum, na.rm = TRUE)
  oe_counts <- tapply(trans$count, oe_ids_trans, sum, na.rm = TRUE)

  # Assign summarized counts, using 0 for missing interactions
  x$bait.trans.count <- as.vector(bait_counts[match(bait_ids, names(bait_counts))])
  x$bait.trans.count[is.na(x$bait.trans.count)] <- 0

  x$target.trans.count <- as.vector(oe_counts[match(oe_ids, names(oe_counts))])
  x$target.trans.count[is.na(x$target.trans.count)] <- 0

  return(x)
})




#' Filter links for further analysis
#' @aliases filterLinks
#' @param x A linkSet object
#' @param filter_intra Whether to filter intra-chromosomal interactions
#' @param filter_unannotate Whether to filter unannotated interactions
#' @param distance The maximum distance between bait and other end
#' @return A linkSet object with filtered interactions
#' @examples
#' data(linkExample)
#' linkSet = c(linkExample,linkExample)
#' linkSet = countInteractions(linkSet)
#' linkSet = filterLinks(linkSet, filter_intra = FALSE, filter_unannotate = FALSE, distance = 100000)
#' @export
setMethod("filterLinks", "linkSet", function(x, filter_intra = TRUE,
                                            filter_unannotate = TRUE,
                                            distance = NULL) {
  # Ensure inter_type is annotated
  if (!.exist_inter(x)) {
    x <- annotateInter(x)
  }

  # Ensure distance is calculated
  if (is.null(x$distance)) {
    x <- pairdist(x)
  }

  # Initialize a logical vector for filtering
  keep <- rep(TRUE, length(x))

  if (filter_intra) {
    keep <- keep & (x$inter_type == "inter")
  }

  if (filter_unannotate) {
    bait_regions <- regionsBait(x)
    keep <- keep & (as.character(seqnames(bait_regions)) != "chrNULL")
  }

  if (!is.null(distance)) {
    keep <- keep & (!is.na(x$distance) & x$distance <= distance)
  }

  # Apply the filter
  filtered_x <- x[keep]

  # Clean unused regions
  filtered_x <- clean_unused_regions(filtered_x)

  return(filtered_x)
})

#' Cross gene enhancer
#' @aliases crossGeneEnhancer
#' @param x A linkSet object
#' @param score_threshold The minimum score to filter interactions
#' @return A linkSet object with filtered interactions
#' @examples
#' data(linkExample)
#' linkSet = c(linkExample,linkExample)
#' linkSet = countInteractions(linkSet)
#' linkSet = filterLinks(linkSet, filter_intra = FALSE, filter_unannotate = FALSE, distance = 100000)
#' linkSet = crossGeneEnhancer(linkSet, score_threshold = 10)
#' @export
setMethod("crossGeneEnhancer", "linkSet", function(x, score_threshold = NULL) {
  if (!is.null(score_threshold) & !"score" %in% colnames(mcols(x))) {
    warning("score column not found.")
    score_threshold <- NULL
  }
  else if (!is.null(score_threshold)) {
    x <- x[mcols(x)$score >= score_threshold]
  }

  # Get bait and oe regions
  baitName <- bait(x)
  oeName <- paste(oe(x))

  # Use data.table to handle large data efficiently
  baitDf <- data.table::as.data.table(list(oeName = oeName, baitName = baitName))
  baitDf <- data.table::setDT(baitDf)[, list(.N), by = list(oeName, baitName)]
  baitDf <- baitDf[get("N") > 0]

  enhancerDf <- data.table::setDT(baitDf)[, list(.N), by = list(oeName)]

  # Assign the Freq column to the mcol of x
  mcols(x)$crossFreq <- enhancerDf[["N"]][match(oeName, enhancerDf[["oeName"]])]

  return(x)
})


#' Order linkSet by mcols
#' @aliases orderLinks
#' @param x A linkSet object
#' @param by The column name to order by
#' @param decreasing Whether to sort in decreasing order
#' @return A linkSet object with ordered interactions
#' @examples
#' data(linkExample)
#' linkSet = c(linkExample,linkExample)
#' linkSet = countInteractions(linkSet)
#' linkSet = filterLinks(linkSet, filter_intra = FALSE, filter_unannotate = FALSE, distance = 100000)
#' linkSet = orderLinks(linkSet, by = "count", decreasing = TRUE)
#' @export
setMethod("orderLinks", "linkSet", function(x, by = "count", decreasing = TRUE) {
  if (!by %in% colnames(mcols(x))) {
    stop("Column '", by, "' not found in mcols of linkSet object.")
  }
  order_indices <- order(mcols(x)[[by]], decreasing = decreasing)
  return(x[order_indices])
})
