% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{withTxDb,character,function-method}
\alias{withTxDb,character,function-method}
\title{Execute Database Operation with Automatic Connection Management}
\usage{
\S4method{withTxDb}{character,function}(x, expr, ...)
}
\arguments{
\item{x}{The genome name or object to operate on}

\item{expr}{Expression to evaluate with database connection}

\item{...}{Additional arguments}
}
\value{
Result of the database operation
}
\description{
Execute Database Operation with Automatic Connection Management
}
\examples{
requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")
requireNamespace("org.Hs.eg.db")
result <- withTxDb("hg38", function(src) {
  genes <- Organism.dplyr::genes(src)
  return(genes)
})
print(result)
}
