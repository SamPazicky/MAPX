% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/X.AI.R
\name{X.AI}
\alias{X.AI}
\title{X.AI}
\usage{
X.AI(
  data = NULL,
  complexes = list(),
  quality.data = NULL,
  index = "semirandom",
  noPCs = 0,
  trials = 500,
  cluster.cut = 0.65,
  weights = NULL,
  y.trans = TRUE
)
}
\arguments{
\item{data}{Data frame with features reduce by PCA. Columns are 'protein', and PC columns starting with 'PC'. Additional columns
'condition' and 'replicate' can be used.}

\item{complexes}{Named list: A list where each element is vector with proteins and the name of that element is the name of the protein complex.}

\item{quality.data}{Data.frame with a column for protein and a quality indicator.}

\item{index}{Character string: Either 'random' or 'semirandom'. 'random' will compare within-complex Euclidean distances with completely random
distances. 'semirandom' will always use PC coordinates of a each subunit and compare them to random coordinates.}

\item{noPCs}{Integer: How many principal components should be used for the assembly index calculation?}

\item{trials}{Integer: How many random Euclidean distances should be calculated for the statistics?}

\item{cluster.cut}{Numeric: At what height should the clustered data from each complex be cut? The largest cluster will be used to calculate
the complex centroid, distance from which is further calculated and gives}

\item{weights}{Numeric vector: Weights for averaging principal components for Euclidean distance calculations.}

\item{y.trans}{Logical: Should the assembly index be log2 transformed? Default is TRUE.}
}
\value{
PCA data and different supporting information.
}
\description{
Calculate assembly indeces from complex features reduce into a PCA plot.
}
\examples{
all.fitdata.frame <- data.frame()
for(tp in timepoints) {
  for(rep in reps) {
    all.fitdata.frame <- bind_rows(all.fitdata.frame,
                                   all.fitdata[[tp]][[rep]] \%>\% mutate(replicate=as.character(rep)) \%>\% mutate(condition=tp)
    )
  }
}
quality.data <- all.fitdata.frame \%>\%
  dplyr::select(protein,condition,replicate,R2) \%>\%
  rename(quality=R2)

data.AIs <- X.AI(data.PCAs$data,complexes, quality.data)

}
