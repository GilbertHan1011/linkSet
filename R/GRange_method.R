#######################################
# Internal functions for use in GRanges methods for paired regions.

.generate_regions <- function(x, FUN, args, other.args, region=region) {
    a1 <- anchor1(x)
    a2 <- anchor2(x)
    N1 <- length(a1)
    N2 <- length(a2)

    # Expanding all arguments to the requested length.
    arg1 <- arg2 <- args
    for (arg in names(args)) { 
        current <- args[[arg]]

        if (is.list(current)) { 
            c1 <- current[[1]]
            c2 <- current[[2]]            
        } else {
            c1 <- c2 <- current
        }

        c1 <- .expand_to_length(c1, N1)
        c2 <- .expand_to_length(c2, N2)
        arg1[[arg]] <- c1
        arg2[[arg]] <- c2
    }

    # Minimizing expansion by identifying unique modifications for each anchor region.
    reg <- regions(x)
    if (region == "bait" || region == "both") {
        reg1 <- .apply_unique_mods(reg, a1, FUN, arg1, other.args)
    } else {
        reg1 <- reg
    }
    if (region == "oe" || region == "both") {
        reg2 <- .apply_unique_mods(reg, a2, FUN, arg2, other.args)
    } else {
        reg2 <- reg
    }
    mod.reg <- .collate_GRanges(reg1, reg2)

    # Fiddling with the input/output ranges.
    unchecked_regions(x) <- mod.reg$ranges
    unchecked_anchor1(x) <- mod.reg$indices[[1]][a1]
    unchecked_anchor2(x) <- mod.reg$indices[[2]][a2]
    x <- clean_unused_regions(x)
    return(x)
}

.expand_to_length <- function(v, N) {
    tmp <- rep(v, length.out=N)
    tmp[] <- v # to trigger warning upon recycling a non-multiple vector.
    return(tmp)
}

.apply_unique_mods <- function(regions, index, FUN, args, other.args) {
    to.order <- c(list(index), args)
    o <- do.call(order, to.order)

    # Identifying unique modifications.
    mod <- integer(0) 
    if (length(o)) { 
        is.uniq <- FALSE
        for (element in to.order) {
            element <- element[o]
            is.uniq <- is.uniq | (element[-1]!=element[-length(o)])
        }
        mod <- o[c(TRUE, is.uniq)]
    }

    # Applying the modifications.
    mod.args <- lapply(args, "[", mod)
    to.mod <- index[mod]
    regions[to.mod] <- do.call(FUN, c(list(x=regions[to.mod]), mod.args, other.args))
    return(regions)
}

########################################
#' linkSet-GRange-Methods
#' @export
#' @rdname linkSet-GRange-Methods
#' @aliases trim
#' @description This man page documents intra range transformations of a [linkSet-class] object.
#' 
#' @param x A linkSet object
#' @param use.names A logical indicating whether to use names
#' @return A linkSet object
#' @author Gilbert Han
#' 
#' 
setMethod("trim", "linkSet", function(x, use.names=TRUE) {
    regions(x) <- trim(regions(x), use.names=use.names)
    return(x)
})

#' @export
#' @rdname linkSet-GRange-Methods
#' @aliases resize
#' 
setMethod("resize", "linkSet", function(x, width, fix="start", use.names=TRUE,...) {
    .generate_regions(x, FUN=resize, args=list(width=width, fix=fix), other.args=list(use.names=use.names, ...), region="both")
})

#' @export
#' @rdname linkSet-GRange-Methods
#' @aliases resizeRegions
setMethod("resizeRegions", "linkSet", function(x, width, fix="start", use.names=TRUE, region = "both",...) {
    .generate_regions(x, FUN=resize, args=list(width=width, fix=fix), other.args=list(use.names=use.names, ...), region=region)
})

#' @export
#' @rdname linkSet-GRange-Methods
#' @aliases narrow
setMethod("narrow", "linkSet", function(x, start=NA, end=NA, width=NA, use.names=TRUE) {
    .generate_regions(x, FUN=narrow, args=list(start=start, end=end, width=width), other.args=list(use.names=use.names), region="both")
})  

#' @export
#' @rdname linkSet-GRange-Methods
#' @aliases narrowRegions
setMethod("narrowRegions", "linkSet", function(x, start=NA, end=NA, width=NA, use.names=TRUE, region="both") {
    .generate_regions(x, FUN=narrow, args=list(start=start, end=end, width=width), other.args=list(use.names=use.names), region=region)
})

#' @export
#' @rdname linkSet-GRange-Methods
#' @aliases shift
setMethod("shift", "linkSet", function(x, shift=0L, use.names=TRUE) {
    .generate_regions(x, FUN=IRanges::shift, args=list(shift=shift), other.args=list(use.names=use.names), region="both")
})

#' @export
#' @rdname linkSet-GRange-Methods
#' @aliases shiftRegions
setMethod("shiftRegions", "linkSet", function(x, shift=0L, use.names=TRUE, region="both") {
    .generate_regions(x, FUN=IRanges::shift, args=list(shift=shift), other.args=list(use.names=use.names), region=region)
})

#' @export
#' @rdname linkSet-GRange-Methods
#' @aliases flank
setMethod("flank", "linkSet", function(x, width, start=TRUE, both=FALSE, use.names=TRUE, ignore.strand=FALSE) {
    .generate_regions(x, FUN=flank, args=list(width=width, start=start), 
                        other.args=list(use.names=use.names, ignore.strand=ignore.strand), region="both")
})

#' @export
#' @rdname linkSet-GRange-Methods
#' @aliases flankRegions
setMethod("flankRegions", "linkSet", function(x, width, start=TRUE, both=FALSE, use.names=TRUE, ignore.strand=FALSE,region="both") {
    .generate_regions(x, FUN=flank, args=list(width=width, start=start), 
                        other.args=list(use.names=use.names, ignore.strand=ignore.strand), region=region)
})

#' @export
#' @rdname linkSet-GRange-Methods
#' @aliases promoters
setMethod("promoters", "linkSet", function(x, upstream=2000, downstream=200, use.names=TRUE) {
    .generate_regions(x, FUN=promoters, args=list(upstream=upstream, downstream=downstream), 
                        other.args=list(use.names=use.names), region="both")
})

#' @export
#' @rdname linkSet-GRange-Methods
#' @aliases promoterRegions
setMethod("promoterRegions", "linkSet", function(x, upstream=2000, downstream=200, use.names=TRUE,region="both") {
    .generate_regions(x, FUN=promoters, args=list(upstream=upstream, downstream=downstream), 
                        other.args=list(use.names=use.names), region=region)
})


#' @export
#' @rdname linkSet-GRange-Methods
#' @aliases width
#' @examples
#' data(linkExample)
#' resize_bait <- resizeRegions(linkExample, width = 75, fix = "start", region = "bait")
#' resize_bait
#' 
#' narrow_bait <- narrowRegions(linkExample, start = 1, width = 5, region = "bait")
#' narrow_bait
#' 
#' shift_oe <- shiftRegions(linkExample, shift = 10, region = "oe")
#' shift_oe
#' 
#' flank_bait <- flankRegions(linkExample, width = 100, start = TRUE, both = FALSE, use.names = TRUE, ignore.strand = FALSE, region = "bait")
#' flank_bait
#' 
#' width(linkExample)
#' 
setMethod("width", "linkSet", function(x) {
    w <- width(regions(x))          
    list(bait=w[anchor1(x)], oe=w[anchor2(x)])
})

## reduce
#' Reduce a linkSet object
#'
#' This function reduces the bait and/or oe regions of a linkSet object and optionally counts interactions,
#' while maintaining the original length of the linkSet.
#'
#' @param x A linkSet object
#' @param reduceBait Logical, whether to reduce bait regions (default: TRUE)
#' @param reduceOE Logical, whether to reduce other end (oe) regions (default: TRUE)
#' @param countInteractions Logical, whether to count interactions after reducing (default: TRUE)
#' @param ... Additional arguments passed to GenomicRanges::reduce
#'
#' @return A reduced linkSet object with the same length as the input
#' @export
#'
#' @importFrom GenomicRanges reduce findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors subjectHits
#'
#' @examples
#' data(linkExample)
#' reduced_ls <- reduceRegions(linkExample, region = "both", countInteractions = TRUE)
#' reduced_ls
#'

setMethod("reduceRegions", "linkSet", function(x, region = "both", 
            countInteractions = TRUE, ...) {

  original_bait <- regionsBait(x)
  original_oe <- oe(x)

  if (region == "bait" || region == "both") {
    bait_reduced <- GenomicRanges::reduce(original_bait, ...)
    bait_overlaps <- findOverlaps(original_bait, bait_reduced)
    new_bait_regions <- bait_reduced[S4Vectors::subjectHits(bait_overlaps)]
    regionsBait(x) <- new_bait_regions
  }
  if (region == "oe" || region == "both") {
    oe_reduced <- GenomicRanges::reduce(original_oe, ...)
    oe_overlaps <- findOverlaps(original_oe, oe_reduced)
    new_oe <- oe_reduced[S4Vectors::subjectHits(oe_overlaps)]
    oe(x) <- new_oe
  }


  if (countInteractions) {
    x <- countInteractions(x)
  }
  return(x)
})

#' @export
#' @rdname linkSet-GRange-Methods
#' @aliases reduce
setMethod("reduce", "linkSet", function(x, ...) {
    reduceRegions(x, region = "both", countInteractions = TRUE, ...)
})