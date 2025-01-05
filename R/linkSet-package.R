#' @keywords internal
"_PACKAGE"

#' @name linkSet
#'
#' @aliases linkSet
#'
#' @title linkSet: Base Classes for Storing Genomic Link Data
#'
#' @description
#' The linkSet package provides tools for working with genomic link sets,
#' which represent connections between different genomic regions. This package
#' is designed for bioinformatics and genomic data analysis, offering various
#' methods to manipulate and analyze linkSet objects.
#'
#' @details
#' The main class provided by this package is the `linkSet` class, which is
#' designed to represent and analyze genomic interactions, particularly
#' focusing on gene-enhancer relationships. Key features include:
#'
#' \itemize{
#'   \item Representation of genomic interactions with two types of anchors:
#'         "bait" (typically genes) and "other end" (typically enhancers or
#'         other regulatory elements).
#'   \item Flexible input methods, supporting construction from various data types.
#'   \item Metadata storage for additional information about interactions.
#'   \item Integration with Bioconductor classes and tools.
#'   \item Methods for annotating promoters and distinguishing between inter-
#'         and intra-chromosomal interactions.
#' }
#'
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicRanges flank narrow promoters resize trim
#' @importFrom methods new setClass setGeneric setMethod show
#' @importFrom S4Vectors Rle DataFrame make_zero_col_DFrame mcols mcols<-
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom BiocGenerics start end width
#' @importFrom Organism.dplyr select
#' @importFrom foreach %dopar%
#' @importFrom rlang .data
#' @references
#' Add any relevant references here.
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/GilbertHan1011/linkSet}
#'   \item Report bugs at \url{https://github.com/GilbertHan1011/linkSet/issues/new}
#' }
#'
#' @examples
#' # Basic usage example
#' library(linkSet)
#' # Add a simple example here
#'
NULL
