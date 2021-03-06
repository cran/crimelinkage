% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/naivebayes.R
\name{predict.naiveBayes}
\alias{predict.naiveBayes}
\title{Generate prediction (sum of log bayes factors) from a \code{naiveBayes} object}
\usage{
\method{predict}{naiveBayes}(object, newdata, components = FALSE,
  vars = NULL, ...)
}
\arguments{
\item{object}{a naive bayes object from \code{\link{naiveBayes}}}

\item{newdata}{data frame of new predictors, column names must match NB names}

\item{components}{(logical) return the log bayes factors from each component
  or return the sum of log bayes factors}

\item{vars}{the names or column numbers of specific predictors. If NULL, then
  all predictors will be used}

\item{\ldots}{not currently used}
}
\value{
BF if \code{components = FALSE}, the sum of log bayes factors, if
   \code{components = TRUE} the component bayes factors (useful for plotting).

   It will give a warning, but still produce output if X is missing predictors.
   The output in this situation will be based on the predictors that are in X.
}
\description{
This does not include the log prior odds, so will be off by a
   constant.
}
\examples{
# See vignette: "Statistical Methods for Crime Series Linkage" for usage.
}
\seealso{
\code{\link{naiveBayes}}, \code{\link{plot.naiveBayes}}
}

