###############################################################
# Simple getters and setters (some of these are exported):

#' linkSet-accessors
#' @rdname linkSet-accessors
#' @aliases anchor1
#' @description
#' Methods to get and set fields in an linkSet object.
#' @return For the getters, values in various slots of x are returned, 
#' while for the setters, the slots of x are modified accordingly – see Details. 
#' @param x A linkSet object
#' @author Gilbert Han
#' @examples
#' data(linkExample)
#' anchor1(linkExample)
#' @export
setMethod("anchor1", "linkSet", function(x) { x@anchor1 })

#' 
#' @rdname linkSet-accessors
#' @aliases anchor2
#' @param x A linkSet object
#' @examples
#' data(linkExample)
#' anchor2(linkExample)
#' @export
setMethod("anchor2", "linkSet", function(x) { x@anchor2 })

#' 
#' @rdname linkSet-accessors
#' @aliases regions
#' @param x A linkSet object
#' @return A vector of the regions
#' @examples
#' data(linkExample)
#' regions(linkExample)
#' @export
setMethod("regions", "linkSet", function(x) {x@regions})


###############################################################
# Seqinfo getting and setting.
#' @rdname linkSet-accessors
#' @aliases seqinfo
#' @param x A linkSet object
#' @export
setMethod("seqinfo", "linkSet", function(x) {
  seqinfo(regions(x))
})

###############################################################
# Exported getters for anchors.
#' @rdname linkSet-accessors
#' @aliases anchorIds
#' @description
#' This method returns the anchor IDs of a linkSet object.
#' 
#' @param x A linkSet object.
#' @param type The type of anchor to return. Can be "both", "bait", or "oe".
#' @return A list of anchor IDs.
#' @examples
#' data(linkExample)
#' anchorIds(linkExample, type="both")
#' @export
setMethod("anchorIds", "linkSet", function(x, type="both") {
  type <- match.arg(type, c("both", "bait", "oe"))
  if (type=="both") {
    out <- list(bait=anchor1(x), oe=anchor2(x))
    names(out$bait) <- names(out$oe) <- names(x)
  } else if (type=="bait") {
    out <- anchor1(x)
    names(out) <- names(x)
  } else {
    out <- anchor2(x)
    names(out) <- names(x)
  }
  return(out)
})


#' @rdname linkSet-accessors
#' @aliases anchors
#' @description
#' This method returns the anchors of a linkSet object.
#' 
#' @param x A linkSet object.
#' @param type The type of anchor to return. Can be "both", "bait", or "oe".
#' @param id If TRUE, returns the anchor IDs instead of the anchors.
#' @return A list of anchors or anchor IDs.
#' @examples
#' data(linkExample)
#' anchors(linkExample, type="both", id=FALSE)
#' @export
setMethod("anchors", "linkSet", function(x, type="both", id=FALSE) {
  if (id) {
    return(anchorIds(x, type=type))
  }

  type <- match.arg(type, c("both", "bait", "oe"))
  if (type=="both") {
    out <- list(bait=x@nameBait, oe=regions(x)[anchor2(x)])
    names(out$bait) <- names(out$oe) <- names(x)
  } else if (type=="bait") {
    out <- x@nameBait
    names(out) <- names(x)
  } else {
    out <- regions(x)[anchor2(x)]
    names(out) <- names(x)
  }
  return(out)
})

# Defining some convenience methods.
#' @export
#' @rdname linkSet-accessors
#' @aliases first
#' @description
#' This method returns the bait anchors of a linkSet object.
#' 
#' @param x A linkSet object.
#' @return A GRanges object containing the bait anchors.
setMethod("first", "linkSet", function(x) { anchors(x, type="bait") })

#' @export
#' @rdname linkSet-accessors
#' @aliases second
#' @description
#' This method returns the other end (oe) anchors of a linkSet object.
#' 
#' @param x A linkSet object.
#' @return A GRanges object containing the oe anchors.
setMethod("second", "linkSet", function(x) { anchors(x, type="oe") })

#' @export
#' @rdname linkSet-accessors
#' @aliases bait
#' @description
#' This method is an alias for 'first' and returns the bait anchors of a linkSet object.
#' 
#' @param x A linkSet object.
#' @return A GRanges object containing the bait anchors.
setMethod("bait", "linkSet", function(x) { anchors(x, type="bait") })

#' @export
#' @rdname linkSet-accessors
#' @aliases oe
#' @description
#' This method is an alias for 'second' and returns the other end (oe) anchors of a linkSet object.
#' 
#' @param x A linkSet object.
#' @return A GRanges object containing the oe anchors.
setMethod("oe", "linkSet", function(x) { anchors(x, type="oe") })

#' @export
#' @rdname linkSet-accessors
#' @aliases regions
#' @description
#' This method returns the regions of a linkSet object.
#' 
#' @param x A linkSet object.
#' @return A GRanges object containing the regions.
setMethod("regions","linkSet",function(x){
  return(x@regions)
})

#' @export
#' @rdname linkSet-accessors
#' @aliases regionsBait
#' @description
#' This method returns the regions corresponding to the bait anchors of a linkSet object.
#' 
#' @param x A linkSet object.
#' @return A GRanges object containing the regions corresponding to the bait anchors.
setMethod("regionsBait", "linkSet", function(x) {
  if (length(x@anchor1) == 0) {
    return(NULL)
  }
  regions(x)[anchor1(x)]
})

#' @export
#' @rdname linkSet-accessors
#' @aliases bait
#' @description
#' This method replaces the bait anchors of a linkSet object with new values.
#' 
#' @param x A linkSet object.
#' @param value A GRanges object containing the new bait anchors.
#' @return The modified linkSet object with the new bait anchors.
setReplaceMethod("bait", "linkSet", function(x, value) {
  x@nameBait <- value
  return(x)
})

#' @export
#' @rdname linkSet-accessors
#' @aliases unchecked_regions
#' @description
#' This method replaces the regions of a linkSet object with new values.
#' 
#' @param x A linkSet object.
#' @param value A GRanges object containing the new regions.
#' @return The modified linkSet object with the new regions.
setReplaceMethod("unchecked_regions", "linkSet", function(x, value) {
    x@regions <- value
    return(x)        
})

#' @rdname linkSet-accessors
#' @aliases unchecked_anchor1
#' @description
#' This method replaces the anchor1 of a linkSet object with new values.
#' 
#' @param x A linkSet object.
#' @param value A vector containing the new anchor1 values.
#' @return The modified linkSet object with the new anchor1 values.
setReplaceMethod("unchecked_anchor1", "linkSet", function(x, value) {
    x@anchor1 <- value 
    return(x)        
})

#' @rdname linkSet-accessors
#' @aliases unchecked_anchor2
#' @description
#' This method replaces the anchor2 of a linkSet object with new values.
#' 
#' @param x A linkSet object.
#' @param value A vector containing the new anchor2 values.
#' @return The modified linkSet object with the new anchor2 values.
setReplaceMethod("unchecked_anchor2", "linkSet", function(x, value) {
    x@anchor2 <- value 
    return(x)        
})

#' @rdname linkSet-accessors
#' @aliases regions
#' @description
#' This method replaces the regions of a linkSet object with new values.
#' 
#' @param x A linkSet object.
#' @param value A GRanges object containing the new regions.
#' @return The modified linkSet object with the new regions.
setReplaceMethod("regions", "linkSet", function(x, value) {
  if (length(value)!=length(regions(x))) { 
      stop("assigned value must be of the same length as 'regions(x)'")
  }
  out <- .resort_regions(anchor1(x), anchor2(x), value)
  unchecked_anchor1(x) <- out$anchor1
  unchecked_anchor2(x) <- out$anchor2
  unchecked_regions(x) <- out$regions
  validObject(x)
  return(x)
})



#' @rdname linkSet-accessors
#' @aliases regionsBait
#' @description
#' This method replaces the regions corresponding to the bait anchors of a linkSet object with new values.
#' 
#' @param x A linkSet object.
#' @param value A GRanges object containing the new regions corresponding to the bait anchors.
#' @return The modified linkSet object with the new regions corresponding to the bait anchors.
setReplaceMethod("regionsBait", "linkSet", function(x, value) {
  if (!is(value, "GRanges")) {
    stop("The 'value' must be a GRanges object")
  }
  if (length(value) != length(oe(x))) {
    stop("The length of 'value' must be equal to the number of bait regions")
  }
  metadata <- mcols(value)
  mcols(value) <- NULL
  if (length(anchor1(x)) == 0) {
    x@regions <- c(regions(x),value)
    x@anchor1 <- (length(regions(x))-length(value)+1) : length(regions(x))
  } else {
    regions(x)[anchor1(x)] <- value
  }
  mcols(x) <- cbind(mcols(x),metadata)
  return(x)
})

#' @rdname linkSet-accessors
#' @aliases oe
#' @description
#' This method replaces the other end (oe) anchors of a linkSet object with new values.
#' 
#' @param x A linkSet object.
#' @param value A GRanges object containing the new other end anchors.
#' @return The modified linkSet object with the new other end anchors.
setReplaceMethod("oe", "linkSet", function(x, value) {
  if (!is(value, "GRanges")) {
    stop("The 'value' must be a GRanges object")
  }
  if (length(value) != length(regions(x)[anchor2(x)])) {
    stop("The length of 'value' must be equal to the number of other end regions")
  }
  metadata <- mcols(value)
  mcols(value) <- NULL
  regions(x)[anchor2(x)] <- value
  mcols(x) <- cbind(mcols(x),metadata)
  return(x)
})

###############################################################
# Defining some other getters and setters.
#' @rdname linkSet-accessors
#' @aliases $
#' @description
#' This method returns the metadata column of a linkSet object.
#' 
#' @param x A linkSet object.
#' @param name A character string specifying the name of the metadata column to return.
#' @return The value of the specified metadata column.
setMethod("$", "linkSet", function(x, name) {
    return(mcols(x)[[name]])
})

#' @rdname linkSet-accessors
#' @aliases $
#' @description
#' This method replaces the metadata column of a linkSet object with new values.
#' @importFrom S4Vectors mcols mcols<-
#' @param x A linkSet object.
#' @param name A character string specifying the name of the metadata column to replace.
#' @param value The new value to assign to the specified metadata column.
#' @return The modified linkSet object with the new metadata column value.
setReplaceMethod("$", "linkSet", function(x, name, value) {
    mcols(x)[[name]] <- value
    return(x)
})


###############################################################
# Name getting and setting.

setMethod("names", "linkSet", function(x) { 
    x@NAMES 
})

setReplaceMethod("names", "linkSet", function(x, value) {
    if (!is.null(value) && !is.character(value)) { value <- as.character(value) }                
    x@NAMES <- value
    validObject(x)
    return(x)
})
