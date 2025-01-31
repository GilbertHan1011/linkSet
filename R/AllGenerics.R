# getset.R
#' @export
setGeneric("anchors", function(x, ...) standardGeneric("anchors"))

#' @export
setGeneric("anchor1", function(x) standardGeneric("anchor1"))

#' @export
setGeneric("anchor2", function(x) standardGeneric("anchor2"))

#' @export
setGeneric("regions", function(x, ...) standardGeneric("regions"))

#' @export
setGeneric("bait", function(x) standardGeneric("bait"))

#' @export
setGeneric("oe", function(x) standardGeneric("oe"))

#' @export
setGeneric("anchorIds", function(x, ...) standardGeneric("anchorIds"))



## set functions for linkSet
#' @export
#' @rdname linkSet-accessors
setGeneric("bait<-", function(x, value) standardGeneric("bait<-"))

#' @export
#' @rdname linkSet-accessors
setGeneric("regions<-", function(x, value) standardGeneric("regions<-"))

#' @export
#' @rdname linkSet-accessors
setGeneric("anchor1<-", function(x, value) standardGeneric("anchor1<-"))

#' @export
#' @rdname linkSet-accessors
setGeneric("anchor2<-", function(x, value) standardGeneric("anchor2<-"))

#' @export
#' @rdname linkSet-accessors
setGeneric("unchecked_regions<-", function(x, value) standardGeneric("unchecked_regions<-"))
setGeneric("unchecked_anchor1<-", function(x, value) standardGeneric("unchecked_anchor1<-"))
setGeneric("unchecked_anchor2<-", function(x, value) standardGeneric("unchecked_anchor2<-"))

#' @export
setGeneric("subsetBait", function(x, ...) {
  standardGeneric("subsetBait")
})

#' @export
setGeneric("subsetBaitRegion", function(x, ...) {
  standardGeneric("subsetBaitRegion")
})

#' @export
setGeneric("subsetOE", function(x, ...) {
  standardGeneric("subsetOE")
})


# annotate.R
#' @export
setGeneric("annotatePromoter", function(x, genome = "hg38", keyType = "symbol",upstream = 500,overwrite = FALSE,...) standardGeneric("annotatePromoter"))

#' Execute Database Operation with Automatic Connection Management
#' 
#' @title Database Operation with Connection Management
#' @description Executes a database operation while managing the connection lifecycle automatically.
#' 
#' @param x The genome name or object to operate on
#' @param expr Expression to evaluate with database connection
#' @param ... Additional arguments
#' @return Result of the database operation
#' @examples
#' \dontrun{
#' requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' requireNamespace("org.Hs.eg.db")
#' result <- withTxDb("hg38", function(src) {
#'   genes <- Organism.dplyr::genes(src)
#'   return(genes)
#' })
#' print(result)
#' }
#' @export
setGeneric("withTxDb", function(x, expr, ...) {
    standardGeneric("withTxDb")
})


#' Set Other End Anchors for linkSet Object
#' 
#' @title Set Other End (OE) Anchors
#' @description Replace the other end (oe) anchors of a linkSet object with new values
#' 
#' @param x A linkSet object
#' @param value A GRanges object containing the new other end anchors
#' @return The modified linkSet object
#' @examples
#' # Create example data
#' gr1 <- GRanges("chr1", IRanges(1:3, width=1))
#' gr2 <- GRanges("chr1", IRanges(4:6, width=1))
#' ls <- linkSet(gr1, gr2)
#' 
#' # Create new other end anchors
#' new_oe <- GRanges("chr1", IRanges(7:9, width=1))
#' 
#' # Replace other end anchors
#' oe(ls) <- new_oe
#' @export
setGeneric("oe<-", function(x, value) standardGeneric("oe<-"))

#' Set Bait Regions for linkSet Object
#' 
#' @title Set Bait Regions
#' @description Replace the regions corresponding to the bait anchors of a linkSet object
#' 
#' @param x A linkSet object
#' @param value A GRanges object containing the new bait regions
#' @return The modified linkSet object
#' @examples
#' # Create example data
#' gr1 <- GRanges("chr1", IRanges(1:3, width=1))
#' gr2 <- GRanges("chr1", IRanges(4:6, width=1))
#' ls <- linkSet(gr1, gr2)
#' 
#' # Create new bait regions
#' new_bait <- GRanges("chr1", IRanges(7:9, width=1))
#' 
#' # Replace bait regions
#' regionsBait(ls) <- new_bait
#' @export
setGeneric("regionsBait<-", function(x, value) standardGeneric("regionsBait<-"))

# method.R
#' @export
setGeneric("linkSet", function(anchor1, anchor2, specificCol, ...) standardGeneric("linkSet"))

#' Display Detailed Information About a linkSet Object
#' 
#' @title Show linkSet Object Details
#' @description Displays detailed information about a linkSet object, including regions,
#' metadata, and optionally sequence information.
#' 
#' @param object A linkSet object to display
#' @param margin Character string for display margin (default: "")
#' @param print.seqinfo Logical, whether to print sequence information (default: FALSE)
#' @param print.classinfo Logical, whether to print class information (default: FALSE)
#' @param baitRegion Logical, whether to display bait regions (default: FALSE)
#' @param ... Additional arguments
#' 
#' @return None (invisible NULL)
#' @examples
#' gr1 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                ranges = IRanges(start = c(1000, 2000, 3000), width = 100),
#'                strand = "+", symbol = c("BRCA1", "TP53", "NONEXISTENT"))
#' gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                ranges = IRanges(start = c(5000, 6000, 7000), width = 100),
#'                strand = "+")
#' ls <- linkSet(gr1, gr2, specificCol = "symbol")
#' showLinkSet(ls)
#' @export
setGeneric("showLinkSet", function(object, margin="", print.seqinfo=FALSE, print.classinfo=FALSE, baitRegion=FALSE,...) standardGeneric("showLinkSet"))

#' @export
setGeneric("regionsBait", function(x) standardGeneric("regionsBait"))

#' Clean Unused Regions in a linkSet Object
#' 
#' This function removes unused regions from a linkSet object, ensuring that all regions 
#' are referenced by either anchor1 or anchor2.
#' 
#' @param x A linkSet object from which to remove unused regions
#' @return A linkSet object with unused regions removed
#' @examples
#' data(linkExample)
#' linkExample <- clean_unused_regions(linkExample)
#' @export
setGeneric("clean_unused_regions", function(x) standardGeneric("clean_unused_regions"))

# distance.R
#' @export
setGeneric("annotateInter", function(x) standardGeneric("annotateInter"))



#' @export
setGeneric("pairdist", function(x, type="mid") standardGeneric("pairdist"))

#' @export
setGeneric("diagnoseLinkSet", function(x) standardGeneric("diagnoseLinkSet"))



# convert.R
#' @export
setGeneric("Convert", function(x, ...) {
  standardGeneric("Convert")
})

#' Convert GInteractions with bait range and oe ranges to linkSet
#' 
#' @title Convert GInteractions to linkSet with bait annotations
#' @param x A GInteractions object
#' @param geneGr A GRanges object representing genes
#' @param peakGr A GRanges object representing peaks
#' @param ... Additional arguments
#' @return A linkSet object
#' @export
setGeneric("baitGInteractions", function(x, geneGr, peakGr,...) {
  standardGeneric("baitGInteractions")
})

#' Convert linkSet object to GInteractions
#' 
#' @title Convert to GInteractions
#' @param x A linkset object
#' @return A GInteractions object
#' @examples
#' data(linkExample)
#' gi <- as.GInteractions(linkExample)
#' gi
#' @export
setGeneric("as.GInteractions", function(x){
  standardGeneric("as.GInteractions")
})

#' Export linkSet to interBed Format
#' 
#' @title Export linkSet to interBed format
#' @description Exports a linkSet object to a tab-delimited interBed format file
#' 
#' @param x A linkSet object
#' @param outfile Output file path
#' @return None. The function writes to the specified file.
#' @examples
#' data(linkExample)
#' tmpfile <- tempfile(fileext = ".txt")
#' exportInterBed(linkExample, tmpfile)
#' cat(readLines(tmpfile), sep = "\n")
#' @export
setGeneric("exportInterBed", function(x, outfile){
  standardGeneric("exportInterBed")
})

#' Export linkSet to WashU Format
#' 
#' @title Export linkSet to WashU browser format
#' @description Exports a linkSet object to a tab-delimited format compatible with the WashU genome browser
#' 
#' @param x A linkSet object
#' @param outfile Output file path
#' @return None. The function writes to the specified file.
#' @examples
#' data(linkExample)
#' tmpfile <- tempfile(fileext = ".txt")
#' exportWashU(linkExample, tmpfile)
#' cat(readLines(tmpfile), sep = "\n")
#' @export
setGeneric("exportWashU", function(x, outfile){
  standardGeneric("exportWashU")
})

# grange_method.R

#' @export
setGeneric("resizeRegions", function(x,width = 1000, fix = "center", use.names = TRUE,region="both", ...) standardGeneric("resizeRegions"))

#' @export
setGeneric("narrowRegions", function(x,  start = 1L, end = 1000L, width = 1000L, use.names = TRUE,region="both",...) standardGeneric("narrowRegions"))

#' @export
setGeneric("shiftRegions", function(x, width = 1000, shift = 0L, use.names = TRUE,region="both",...) standardGeneric("shiftRegions"))

#' @export
setGeneric("flankRegions", function(x,width = 1000,start=TRUE, both=FALSE, use.names=TRUE,ignore.strand=FALSE,region="both", ...) standardGeneric("flankRegions"))

#' @export
setGeneric("promoterRegions", function(x,upstream=2000, downstream=200, use.names=TRUE,region="both", ...) standardGeneric("promoterRegions"))



# count.R

#' Count Bait and Other End Interactions
#'
#' This function takes a linkSet object and counts the number of interactions 
#' for each bait and other end.
#'
#' @param x A linkSet object
#' @param baitRegions Whether to count bait regions (default: TRUE)
#' @return A linkSet object with counts for each unique interaction
#' @examples
#' data(linkExample)
#' linkSet = c(linkExample,linkExample)
#' linkSet = countInteractions(linkSet)
#' linkSet
#' @aliases countBaitOe
#' @export
setGeneric("countInteractions", function(x, baitRegions = TRUE) standardGeneric("countInteractions"))

#' Count Interaction Interactibility
#' 
#' @title Count bait and oe interactibility
#' @description This function calculates the number of trans interactions for each bait and oe. 
#' The word "interactibility" can refer to https://doi.org/10.1038%2Fnature11279.
#' 
#' @param x A linkSet object
#' @param baitRegions Whether to count bait regions (default: TRUE)
#' @return A linkSet object with counts for each unique interaction
#' @examples
#' data(linkExample)
#' linkSet = c(linkExample,linkExample)
#' linkSet = countInteractions(linkSet)
#' linkSet = countInteractibility(linkSet)
#' @export
setGeneric("countInteractibility", function(x, baitRegions = TRUE) standardGeneric("countInteractibility"))

#' Reduce a linkSet Object
#' 
#' @title Reduce Regions in a linkSet Object
#' @description This function reduces the bait and/or oe regions of a linkSet object 
#' and optionally counts interactions, while maintaining the original length of the linkSet.
#' 
#' @param x A linkSet object
#' @param region Character, specifying which regions to reduce: "both", "bait", or "oe" (default: "both")
#' @param countInteractions Logical, whether to count interactions after reducing (default: TRUE)
#' @param ... Additional arguments passed to GenomicRanges::reduce
#' 
#' @return A reduced linkSet object with the same length as the input
#' @examples
#' data(linkExample)
#' reduced_ls <- reduceRegions(linkExample, region = "both", countInteractions = TRUE)
#' reduced_ls
#' @export
setGeneric("reduceRegions", function(x, region = "both", countInteractions = TRUE, ...) {
  standardGeneric("reduceRegions")
})

#' @export
setGeneric("filterLinks", function(x, filter_intra = TRUE, filter_unannotate = TRUE, distance = NULL) {
  standardGeneric("filterLinks")
})

#' @export
setGeneric("crossGeneEnhancer", function(x, score_threshold = NULL) {
  standardGeneric("crossGeneEnhancer")
})
#' @export
setGeneric("orderLinks", function(x, by = "count", decreasing = TRUE) {
  standardGeneric("orderLinks")
})

# plot.R
#' Add Genome Links to Coverage Plot
#' 
#' @title Add Genome Links to Coverage Plot
#' @description Creates a visualization of genomic links for a linkSet object
#' 
#' @param linkSet A linkSet object
#' @param score.col Column name containing score information (default: "count")
#' @param score.threshold Score threshold for filtering links (default: NULL)
#' @param score.color Color vector for score visualization (default: c("grey70", "#56B1F7", "#132B43"))
#' @param scale.range Scale factor for link height (default: 10)
#' @param plot.space Top and bottom margin (default: 0.1)
#' @param plot.height Relative height of link to coverage plot (default: 0.2)
#' @param arrow.size Size of arrow heads (default: 0.05)
#' @param remove_x_axis Whether to remove x-axis (default: FALSE)
#' @param link_plot_on_top Whether to plot links above coverage (default: FALSE)
#' @param extend.base Base pair extension range (default: 10000)
#' @param show.rect Whether to show rectangle borders (default: FALSE)
#' @param x.range Range for x-axis (default: NULL)
#' @param log.scale Whether to use log scale for scores (default: TRUE)
#' 
#' @return A ggplot layer object
#' @export
setGeneric("geom_linkset", 
    function(linkSet, 
             score.col = "count", 
             score.threshold = NULL, 
             score.color = c("grey70", "#56B1F7", "#132B43"),
             scale.range = 10, 
             plot.space = 0.1, 
             plot.height = 0.2, 
             arrow.size = 0.05, 
             remove_x_axis = FALSE,
             link_plot_on_top = FALSE,
             extend.base = 10000,
             show.rect = FALSE, 
             x.range = NULL, 
             log.scale = TRUE) {
        standardGeneric("geom_linkset")
    }
)

#' @export
setGeneric("plot_genomic_ranges", function(linkset, showBait = NULL, showOE = NULL,x.range = NULL,
                                            score.col = "count",
                                            show.rect = TRUE,
                                            extend.base = 10000,
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
                                            arrow.size = 0.05, remove_x_axis = FALSE,
                                            plot.height = 0.4, plot.space = 0.1,
                                            log.scale = TRUE) {
  standardGeneric("plot_genomic_ranges")
})


# statical.R
#' Run ChICANE Analysis on linkSet Object
#' 
#' @title Run ChICANE Analysis
#' @description This function adapts the \code{chicane} function from the \code{ChICANE} 
#' package to work with the \code{linkSet} object format. It runs the full method for 
#' detecting significant interactions in capture Hi-C experiments.
#' 
#' @param linkSet A linkSet object containing interaction data
#' @param ... Additional arguments passed to methods
#' 
#' @return A linkSet object with additional columns:
#' \itemize{
#'   \item{expected}{The expected number of reads linking fragments under the fitted model}
#'   \item{p.value}{P-value for test of observed vs expected read counts}
#'   \item{q.value}{FDR-corrected p-value}
#' }
#' @export
setGeneric("run_chicane", function(linkSet, ...) {
    standardGeneric("run_chicane")
})