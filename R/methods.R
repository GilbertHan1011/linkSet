###############################################################
# Setting validity and show methods.

.check_inputs <- function(anchor1, anchor2, nameBait, regions,same.length=TRUE) {
  if (!all(is.finite(anchor1)) || !all(is.finite(anchor2))) {
    return("all anchor indices must be finite integers")
  }
  if (!all(anchor1 >= 1L) || !all(anchor2 >= 1L)) {
    return('all anchor indices must be positive integers')
  }
  nregs1 <- length(nameBait)
  nregs2 <- length(regions)
  if (!all(anchor1 <= nregs2) || !all(anchor2 <= nregs2)) {
    return("all anchor indices must refer to entries in 'regions'")
  }
  if (same.length && length(nameBait)!=length(anchor2)) {
    return("first and second anchor vectors have different lengths")
  }
  return(TRUE)
}
#
setValidity2("linkSet", function(object) {
  if (is.unsorted(regions(object))) { # Don't move into .check_inputs, as resorting comes after checking validity in various methods.
    return("'regions' should be sorted")
  }
  msg <- .check_inputs(anchor1(object), anchor2(object), bait(object),regions(object))
  if (is.character(msg)) { return(msg) }

  ### Length of anchors versus object is automatically checked by 'parallel_slot_names.'

  if (!is.null(names(object))) {
    if (length(names(object))!=length(object)) {
      stop("'NAMES' must be NULL or have length equal to that of the object")
    }
  }

  return(TRUE)
})


setMethod("parallel_slot_names", "linkSet", function(x) {
  base_slots <- callNextMethod() # Get the base slots from the parent class
  if (length(x@anchor1) == 0) {
    c("anchor2", "nameBait", "NAMES", base_slots)
  } else {
    c("anchor1", "anchor2", "nameBait", "NAMES", base_slots)
  }
})

# For coercion to an environment:
setMethod("parallelVectorNames", "linkSet", function(x) {
  c("anchor1", "anchor2","nameBait", "regions", "names")
})

#' Display detailed information about a linkSet object
#' @export
#' @rdname show
#' @examples
#' # Example usage of show method for linkSet object
#' gr1 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                ranges = IRanges(start = c(1000, 2000, 3000), width = 100),
#'                strand = "+", symbol = c("BRCA1", "TP53", "NONEXISTENT"))
#' gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                ranges = IRanges(start = c(5000, 6000, 7000), width = 100),
#'                strand = "+")
#' ls <- linkSet(gr1, gr2, specificCol = "symbol")
#' show(ls)
setMethod("show", "linkSet", function(object) {
  showLinkSet(object, margin="  ", print.seqinfo=TRUE, print.classinfo=TRUE, baitRegion=FALSE)
})

#' @export
#' @rdname show
#' @param object A linkSet object to be displayed.
#' @param margin A character string to be used as the margin for the display.
#' @param print.seqinfo A logical indicating whether to print sequence information.
#' @param print.classinfo A logical indicating whether to print class information.
#' @param baitRegion A logical indicating whether to display bait regions.
#' @examples
#' # Example usage of showLinkSet method for linkSet object
#' gr1 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                ranges = IRanges(start = c(1000, 2000, 3000), width = 100),
#'                strand = "+", symbol = c("BRCA1", "TP53", "NONEXISTENT"))
#' gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                ranges = IRanges(start = c(5000, 6000, 7000), width = 100),
#'                strand = "+")
#' ls <- linkSet(gr1, gr2, specificCol = "symbol")
#' showLinkSet(ls)
setMethod("showLinkSet", "linkSet",function(x, margin="", print.seqinfo=FALSE, 
                                          print.classinfo=FALSE, baitRegion=FALSE) {
  lx <- length(x)
  nr <- length(regions(x))
  nc <- .safeNMcols(x)
  cat(class(x), " object with ",
      lx, " ", ifelse(lx == 1L, "interaction", "interactions"), " and ",
      nc, " metadata ", ifelse(nc == 1L, "column", "columns"),
      ":\n", sep="")
  
  if (baitRegion && is.null(regionsBait(x))) {
    message("Please annotate bait first.")
    baitRegion <- FALSE
  }
  if (baitRegion) {
    out <- makePrettyMatrixForCompactPrinting(x, function(x) {
      .makeNakedMatFromGInteractions(x, baitRegion=TRUE)
    })
  } else {
    out <- makePrettyMatrixForCompactPrinting(x, .makeNakedMatFromGInteractions)
  }

  if (print.classinfo) {
    if (baitRegion) {
      .COL2CLASS <- c(bait = "character", seqnames_bait = "Rle", ranges_bait = "IRanges", 
                      "   " = "", seqnames_oe = "Rle", ranges_oe = "IRanges")
    } else {
      .COL2CLASS <- c(bait = "character", "   " = "", seqnames_oe = "Rle", ranges_oe = "IRanges")
    }
    extraColumnNames <- GenomicRanges:::extraColumnSlotNames(x)
    .COL2CLASS <- c(.COL2CLASS, getSlots(class(x))[extraColumnNames])
    classinfo <- makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
    classinfo[,"   "] <- ""
    stopifnot(identical(colnames(classinfo), colnames(out)))
    out <- rbind(classinfo, out)
  }

  if (nrow(out) != 0L) {
    rownames(out) <- paste0(margin, rownames(out))
  }
  print(out, quote=FALSE, right=TRUE, max=length(out))
  if (print.seqinfo) {
    cat(margin, "-------\n", sep="")
    ncr <-  .safeNMcols(regions(x))
    cat(margin, "regions: ", nr, " ranges and ", ncr, " metadata ", ifelse(ncr==1L, "column", "columns"), "\n", sep="")
    cat(margin, "seqinfo: ", summary(seqinfo(x)), "\n", sep="")
  }
})


.safeNMcols <- function(x) {
  #return column number safely
  nc <- ncol(mcols(x))
  if (is.null(nc)) { nc <- 0L }
  return(nc)
}

.makeNakedMatFromGInteractions <- function(x, baitRegion=FALSE) {
  lx <- length(x)
  nc <- .safeNMcols(x)
  
  if (baitRegion && !is.null(regionsBait(x))) {
    ans <- cbind(.pasteAnchor(anchors(x, type="bait"), append="bait"),
                 .pasteAnchor(regionsBait(x), append="bait"),
                 "   " = rep.int("---", lx),
                 .pasteAnchor(anchors(x, type="oe"), append="oe"))
  } else {
    ans <- cbind(.pasteAnchor(anchors(x, type="bait"), append="bait"),
                 "   " = rep.int("---", lx),
                 .pasteAnchor(anchors(x, type="oe"), append="oe"))
  }
  
  if (nc > 0L) {
    tmp <- do.call(data.frame, c(lapply(mcols(x), showAsCell), list(check.names=FALSE)))
    ans <- cbind(ans, `|` = rep.int("|", lx), as.matrix(tmp))
  }
  ans
}

#' Format anchor information for display
#' @importFrom GenomicRanges seqnames ranges
#' @importFrom S4Vectors showAsCell
#' @keywords internal
.pasteAnchor <- function(x, append) {
  if(is.character(x)){
    out <- as.matrix(x)
    colnames(out) <- "bait"
  } else{
    out <- cbind(as.character(seqnames(x)), showAsCell(ranges(x)))
    colnames(out) <- paste0(c("seqnames", "ranges"),"_", append)
  }
  out
}
###############################################################
# Constructors
.enforce_order <- function(anchor1, anchor2) {
    swap <- anchor2 < anchor1
    if (any(swap)) { 
        temp <- anchor1[swap]
        anchor1[swap] <- anchor2[swap]
        anchor2[swap] <- temp
    }
    return(list(anchor1=anchor1, anchor2=anchor2))
}


.resort_regions <- function(anchor1, anchor2, regions) {
  if (is.unsorted(regions)) {
    o <- order(regions)
    new.pos <- seq_along(o)
    new.pos[o] <- new.pos
    anchor1 <- new.pos[anchor1]
    anchor2 <- new.pos[anchor2]
    regions <- regions[o]
  }
  return(list(anchor1=anchor1, anchor2=anchor2, regions=regions))
}

#' @keywords internal
.new_LK <- function(anchor1, anchor2, nameBait, regions, metadata) {
  elementMetadata <- make_zero_col_DFrame(length(nameBait))

  # Checking odds and ends.
  anchor1 <- as.integer(anchor1)
  anchor2 <- as.integer(anchor2)
  if (is.null(nameBait)){
    nameBait <- paste(gr1)
  }
  msg <- .check_inputs(anchor1, anchor2, nameBait, regions)
  if (is.character(msg)) { stop(msg) }

  # out <- .resort_regions(anchor1, anchor2, regions)
  # anchor1 <- out$anchor1
  # anchor2 <- out$anchor2
  # regions <- out$regions

  cls <- "linkSet"
  #browser()
  new(cls,
      anchor1=anchor1,
      anchor2=anchor2,
      nameBait=nameBait,
      regions=regions,
      elementMetadata=elementMetadata,
      metadata=as.list(metadata))
}

#' @export
setMethod("linkSet", c("character", "GRanges","character_Or_missing"),
          function(anchor1, anchor2, specificCol,metadata=list(),  ...) {
            mcol2 <- mcols(anchor2)
            mcols(anchor2) <- NULL
            colnames(mcol2) <- sprintf("anchor2.%s", colnames(mcol2))
            extraCols <- DataFrame(...)
            if (ncol(extraCols) == 0L) {
              extraCols <- make_zero_col_DFrame(length(anchor1))
            }
            mcolBind <- cbind(extraCols, mcol2)
            nameBait <- anchor1
            anchor1 <- NULL

            collated <- .collate_GRanges(anchor2)
            regions <- collated$ranges
            anchor2 <- collated$indices[[1]]

            out <- .new_LK(anchor1=anchor1, anchor2=anchor2,
                           nameBait=nameBait,regions=regions,
                           metadata=metadata)
            mcols(out) <- mcolBind
            out
          }
)

.collate_GRanges <- function(...) {
  incoming <- list(...)
  obj.dex <- rep(factor(seq_along(incoming)), lengths(incoming))
  combined <- do.call(c, incoming)

  # Sorting and re-indexing.
  o <- order(combined)
  refdex <- integer(length(o))
  refdex[o] <- seq_along(combined)
  combined <- combined[o]

  # Removing duplicates and re-indexing.
  is.first <- !duplicated(combined)
  new.pos <- cumsum(is.first)
  combined <- combined[is.first]
  refdex <- new.pos[refdex]
  return(list(indices=split(refdex, obj.dex), ranges=combined))
}

#' @export
setMethod("linkSet", c("GRanges", "GRanges","character_Or_missing"),
          function(anchor1, anchor2, specificCol,metadata=list(),  ...) {
            # Stripping metadata and putting it somewhere else.
            mcol1 <- mcols(anchor1)
            mcols(anchor1) <- NULL
            colnames(mcol1) <- sprintf("anchor1.%s", colnames(mcol1))
            mcol2 <- mcols(anchor2)
            mcols(anchor2) <- NULL
            colnames(mcol2) <- sprintf("anchor2.%s", colnames(mcol2))

            # Additional Interaction-specific metadata
            extraCols <- DataFrame(...)
            if (ncol(extraCols) == 0L) {
              extraCols <- make_zero_col_DFrame(length(anchor1))
            }

            mcolBind <- cbind(extraCols, mcol1, mcol2)

            if (!missing(specificCol)){
              if (length(specificCol) > 1) {
                if (length(specificCol) == nrow(mcolBind)) {
                  nameBait <- specificCol
                } else {
                  warning("Length of specificCol does not match the number of rows in mcolBind. Using default naming.")
                  nameBait <- paste(anchor1)
                }
              } else {
                specificColName <- paste0("anchor1.",specificCol)
                if (specificColName %in% colnames(mcolBind)){
                  nameBait <- mcolBind[specificColName]
                  nameBait <- unlist(nameBait)
                } else{
                  warning(paste0("Can't find ", specificCol, "in metadata............"))
                  nameBait <- paste(anchor1)
                }
              }
            } else{
              nameBait <- paste(anchor1)
            }

            collated <- .collate_GRanges(anchor1, anchor2)
            regions <- collated$ranges
            anchor1 <- collated$indices[[1]]
            anchor2 <- collated$indices[[2]]

            # regionBait <- anchor1
            # anchor1 <- seq_along(regionBait)
            # regionOE <- anchor2
            # anchor2 <- seq_along(regionOE)

            out <- .new_LK(anchor1=anchor1, anchor2=anchor2,
                           nameBait=nameBait,
                           regions= regions,
                           metadata=metadata)
            mcols(out) <- mcolBind
            out
          }
)


#' Clean Unused Regions in a linkSet Object
#' 
#' This function removes unused regions from a linkSet object, ensuring that all regions are referenced by either anchor1 or anchor2.
#' 
#' @param x A linkSet object from which to remove unused regions.
#' @export
#' @examples
#' data(linkExample)
#' linkExample <- clean_unused_regions(linkExample)
setMethod("clean_unused_regions", "linkSet", function(x) {
    used_regions <- sort(unique(c(anchor1(x), anchor2(x))))
    new_regions <- regions(x)[used_regions]
    
    # Create a mapping from old indices to new indices
    index_map <- integer(length(regions(x)))
    index_map[used_regions] <- seq_along(used_regions)
    
    # Update anchor indices
    new_anchor1 <- index_map[anchor1(x)]
    new_anchor2 <- index_map[anchor2(x)]
    
    # Update the linkSet object
    unchecked_regions(x) <- new_regions
    unchecked_anchor1(x) <- new_anchor1
    unchecked_anchor2(x) <- new_anchor2
    
    return(x)
})

#== subset ==#
#' Subset linkSet object based on bait names
#' @rdname linkSet-subset-methods
#' @aliases subsetBait
#' @param x A linkSet object
#' @param subset A vector of bait names to keep
#'
#' @return A new linkSet object containing only the specified bait interactions
#' @export
#' @examples
#' data(linkExample)
#' subset_bait_names <- c("bait1", "bait2")
#' subsetted_linkSet <- subsetBait(linkExample, subset_bait_names)
setMethod("subsetBait", "linkSet", function(x, subset) {
  idx <- bait(x) %in% subset
  clean_unused_regions(x[idx])
})

#' Subset linkSet object based on bait regions
#' @rdname linkSet-subset-methods
#' @aliases subsetBaitRegion
#' @param x A linkSet object
#' @param subset A GRanges object specifying the regions to keep
#'
#' @return A new linkSet object containing only the interactions with bait regions overlapping the subset
#' @export
#' @examples
#' data(linkExample)
#' subset_bait_regions <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(100, 200), end = c(150, 250)))
#' subsetted_linkSet <- subsetBaitRegion(linkExample, subset_bait_regions)
setMethod("subsetBaitRegion", "linkSet", function(x, subset) {
  bait_regions <- regionsBait(x)
  if (is.null(bait_regions)) {
    stop("Bait regions are not available. Please annotate bait first.")
  }
  if (is.character(subset)) {
    subset <- .convert_to_grange(subset)
  }
  overlaps <- findOverlaps(bait_regions, subset)
  idx <- queryHits(overlaps)
  clean_unused_regions(x[idx])
})

#' Subset linkSet object based on other end (oe) regions
#' @rdname linkSet-subset-methods
#' @aliases subsetOE
#' @param x A linkSet object
#' @param subset A GRanges object specifying the regions to keep
#'
#' @return A new linkSet object containing only the interactions with oe regions overlapping the subset
#' @export
#' @examples
#' data(linkExample)
#' subset_oe_regions <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(300, 400), end = c(350, 450)))
#' subsetted_linkSet <- subsetOE(linkExample, subset_oe_regions)
setMethod("subsetOE", "linkSet", function(x, subset) {
  oe_regions <- oe(x)
  if (is.character(subset)) {
    subset <- .convert_to_grange(subset)
  }
  overlaps <- findOverlaps(oe_regions, subset)
  idx <- queryHits(overlaps)
  ls = x[idx]
  ls = clean_unused_regions(ls)
})
