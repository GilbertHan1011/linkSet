% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{show,linkSet-method}
\alias{show,linkSet-method}
\alias{showLinkSet,linkSet-method}
\title{Display detailed information about a linkSet object}
\usage{
\S4method{show}{linkSet}(object)

\S4method{showLinkSet}{linkSet}(
  x,
  margin = "",
  print.seqinfo = FALSE,
  print.classinfo = FALSE,
  baitRegion = FALSE
)
}
\description{
Display detailed information about a linkSet object
}
\examples{
# Example usage of show method for linkSet object
gr1 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
               ranges = IRanges(start = c(1000, 2000, 3000), width = 100),
               strand = "+", symbol = c("BRCA1", "TP53", "NONEXISTENT"))
gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
               ranges = IRanges(start = c(5000, 6000, 7000), width = 100),
               strand = "+")
ls <- linkSet(gr1, gr2, specificCol = "symbol")
show(ls)
# Example usage of showLinkSet method for linkSet object
gr1 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
               ranges = IRanges(start = c(1000, 2000, 3000), width = 100),
               strand = "+", symbol = c("BRCA1", "TP53", "NONEXISTENT"))
gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
               ranges = IRanges(start = c(5000, 6000, 7000), width = 100),
               strand = "+")
ls <- linkSet(gr1, gr2, specificCol = "symbol")
showLinkSet(ls)
}
