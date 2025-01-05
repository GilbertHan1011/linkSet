#== make linkObj-------------------------
#' LinkSet object
#' @export
#' @aliases LinkSet
#' @description The linkSet object is a container for storing gene-enhancer interactions.
#' @details The linkSet object is a vectors of paired gene-enhancer interactions.
#' @slot nameBait A character vector of the bait names.
#' @slot anchor1 A integer vector of the first anchor indices.
#' @slot anchor2 A integer vector of the second anchor indices.
#' @slot regions A GenomicRanges object of the regions.
#' @slot NAMES A character vector of the region names.
#' @slot elementMetadata A DataFrame of the element metadata.
#' @seealso \code{\link{linkSet}}
#' @examples
#' showClass("linkSet")  # shows the known subclasses
#' 
#' set.seed(7000)
#' N <- 40
#' all.starts <- round(runif(N, 1, 100))
#' all.ends <- all.starts + round(runif(N, 5, 20))
#' all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))
#' genes = c(rep("SP7",4),rep("ASPN",10),rep("XBP1",6))
#' Np <- 20
#' all.anchor1 <- sample(N, Np)
#' gr1 <- all.regions[all.anchor1]
#' gr1$symbol <- genes
#' all.anchor2 <- setdiff(1:40,all.anchor1)
#' gr2 <- all.regions[all.anchor2]
#' x <- linkSet(gr1, gr2,specificCol = "symbol")
#' x
#' x2 <- linkSet(genes, gr2)
#' x2
#' 
setClass("linkSet",
         contains="Vector",
         representation(
           nameBait = "character",
           anchor1="integer",
           anchor2="integer",
           regions="GenomicRanges_OR_missing",
           NAMES="character_OR_NULL",
           elementMetadata="DataFrame"
         )
)

#' Character or Missing Class Union
#'
#' @name character_Or_missing-class
#' @docType class
#' @description A class union of character vectors and missing values used in linkSet
#' package for optional character arguments.
#'
#' @details This class is used internally by the linkSet package to handle optional
#' character arguments, particularly in the linkSet constructor and methods.
#'
#' @export
setClassUnion("character_Or_missing", c("character", "missing"))
# setClassUnion("integer_Or_missing", c("integer", "missing"))
