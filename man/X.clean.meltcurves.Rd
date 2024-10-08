% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/X.clean.meltcurves.R
\name{X.clean.meltcurves}
\alias{X.clean.meltcurves}
\title{X.clean.meltcurves}
\usage{
X.clean.meltcurves(
  data = NULL,
  keratins = TRUE,
  serum = TRUE,
  trypsin = TRUE,
  average.duplicate.proteins = TRUE,
  remove.nas = TRUE
)
}
\arguments{
\item{data}{Data frame with the loaded mass spec data}

\item{serum}{Logical: Should the serum proteins be removed?}

\item{trypsin}{Logical: Should trypsin be removed?}

\item{average.duplicate.proteins}{Logical: Should the duplicate proteins be averaged?}

\item{remove.nas}{Logical: Should the rows with NA values be removed?}

\item{keratin}{Logical: Should the keratins be removed?}
}
\value{
List of two elements. $data is the data frame with contaminants removed and $removed is the data frame listing the removed contaminants.
}
\description{
Cleans the loaded mass spec data
}
\examples{
MCdata.clean <- X.clean.meltcurves(MCdata.raw)
}
