% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formatConverter.R
\name{baitGInteractions,GInteractions,GRanges,GRanges-method}
\alias{baitGInteractions,GInteractions,GRanges,GRanges-method}
\title{Convert GInteractions with bait range and oe ranges to linkSet}
\usage{
\S4method{baitGInteractions}{GInteractions,GRanges,GRanges}(x, geneGr, peakGr, geneSymbol = NULL)
}
\arguments{
\item{geneGr}{A GRanges object representing genes}

\item{peakGr}{A GRanges object representing peaks}

\item{geneSymbol}{A character vector with same length as geneGr or column name in mcols(geneGr) for gene symbols}

\item{gi}{A GInteractions object}
}
\value{
A linkSet object
}
\description{
Convert GInteractions with bait range and oe ranges to linkSet
}
