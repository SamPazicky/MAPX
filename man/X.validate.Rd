% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/X.validate.R
\name{X.validate}
\alias{X.validate}
\title{X.validate}
\usage{
X.validate(
  data = NULL,
  model = NULL,
  predict.type = "prob",
  labels.col = NA,
  eval.metric = "prc",
  plot = FALSE,
  black = FALSE,
  kappa.cutoff = 0.5,
  kappa.weight = "unweighted",
  F1.cutoff = 0.5,
  evaluate = TRUE
)
}
\arguments{
\item{data}{Data frame with measured parameters.}

\item{model}{The model.}

\item{predict.type}{The prediction type of the model (typically 'prob' or 'probability')}

\item{labels.col}{Character string: 'data' column name with labels (1 for complex-forming, 0 for others.)}

\item{eval.metric}{Character string: Method for predictor evaluation. 'roc' for area under the receiver-operator curve, 'prc' for area under the
precision-recall curve, 'kappa' for Cohen's kappa and 'F1' for F1 score.}

\item{plot}{Logical: Should the plots be plotted and outputted?}

\item{black}{Logical: Should the plotted curve be black or color gradient? Default is FALSE.}

\item{kappa.cutoff}{Numeric: What is the score cutoff for evaluating based on Cohen's kappa? Default is 0.5.}

\item{kappa.weight}{Character string: Same as argument weight in irr:kappa2}

\item{F1.cutoff}{Numeric: What is the score cutoff for evaluating based on F1? Default is 0.5.}

\item{evaluate}{Logical: Should the model be evaluated?}
}
\value{
A list with three elements. $eval.metric is a number corresponding to the chosen evaluation method. $curvedata are data that allow plotting
a plot corresponding to the ealuation metric. $plot is the plot of the corresponding metric.
}
\description{
Validate models.
}
\examples{
evaluation <- X.validate(data,model,"score","complex","prc")
}
