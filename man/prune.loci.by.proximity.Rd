% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GREA.functions.R
\name{prune.loci.by.proximity}
\alias{prune.loci.by.proximity}
\title{prune.loci.by.proximity
\code{prune.loci.by.proximity} Given a set of \code{positions}(genomic coordinates), will return a data.frame with the genes within the input \code{positions}, and if pValue is provided than will trim the original \code{positions} leaving only the SNP with the lowest pValue.}
\usage{
prune.loci.by.proximity(positions, window = "2MB",
  position.colnames = c("chr", "pos"), prune = FALSE, pvalCol = NULL)
}
\arguments{
\item{positions}{A data.frame with genomic coordinates.}

\item{position.colnames}{Character vector defining the colnames used for defining the genomic coordinates within the \code{positions} data.frame.}

\item{gene.window}{Numeric, which will deifine the genomic window that will be used to get neighbouring genes given a certain genomic position defined by \code{positions}.}

\item{gene.info}{A data.frame with genomic coordinates for genes. For reference please look at data(gene.info)}

\item{pVal.col}{Character, if provided, then an extra column to the returned data-frame will be added with the pValue of the SNP associated.}
}
\value{
A dataframe
}
\description{
prune.loci.by.proximity
\code{prune.loci.by.proximity} Given a set of \code{positions}(genomic coordinates), will return a data.frame with the genes within the input \code{positions}, and if pValue is provided than will trim the original \code{positions} leaving only the SNP with the lowest pValue.
}
