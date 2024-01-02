\name{assessmentGroup}
\alias{assessmentGroup}
\title{
Assessment Functions
}
\description{
Calculates the specified assessment function for each participant and each condition.
}
\usage{
assessmentGroup(inData, stopping.rule=c("OR", "AND"), correct=c(TRUE, FALSE), 
  fast=c(TRUE, FALSE), detection=TRUE, plotAt=TRUE, ...)
}
\arguments{
  \item{inData}{Data collected from a Double Factorial Paradigm experiment in standard form.}
  \item{OR}{Indicates whether to compare performance to an OR or AND processing baseline.  Provided for backwards compatibility for package version < 2. }
  \item{stopping.rule}{Indicates whether to use OR or AND processing baseline to calculate individual assessment functions.}
  \item{plotAt}{Indicates whether or not to generate plots of the assessment functions.}
  \item{correct}{Indicates whether to assess performance on correct trials.}
  \item{fast}{Indicates whether to use cumulative distribution functions or survivor functions to assess performance.}
  \item{detection}{Indicates whether to use a detection task baseline or a discrimination task baseline.}
  \item{...}{Arguments to be passed to plot function.}
}
\details{
For the details of the assessment functions, see \code{\link{assessment}}.
}
\value{
  A list containing the following components:
  \item{overview}{A data frame indicating whether the OR and AND capacity coefficients significantly above baseline (super), below baseline (limited) or neither (unlimited) both at the individual level and at the group level in each condition.  NA is returned for any participant that had performance below the accuracy cutoff in a condition.}
  \item{Ct.fn}{Matrix with each row giving the values of the of the estimated capacity coefficient for one participant in one condition for values of times.  The rows match the ordering of statistic.}
  \item{Ct.var}{Matrix with each row giving the values of the of the variance of the estimated capacity coefficient for one participant in one condition for values of times.  Only returned if ratio=FALSE.}
  \item{capacity}{A list with the returned values from the capacity function for each participant and each condition.}
  \item{times}{Times at which the matrix capacity coefficients are calculated in Ct.fn matrix.}
}
\references{
Townsend, J.T. and Altieri, N. (2012).  An accuracy-response time capacity assessment function that measures performance against standard parallel predictions. \emph{Psychological Review, 3}, 500-516.

Donkin, C, Little, D.R. and Houpt (2013). Assessing the effects of caution on the capacity of information processing.  \emph{Manuscript submitted for publication.}
}
\author{
Joe Houpt <joseph.houpt@utsa.edu>
}
\seealso{
\code{\link{assessment}}
}
\examples{
\dontrun{
data(dots)
assessmentGroup(subset(dots, Condition=="OR"), 
  stopping.rule="OR", correct=TRUE, fast=FALSE, 
  detection=TRUE)
assessmentGroup(subset(dots, Condition=="AND"), 
  stopping.rule="AND", correct=TRUE, fast=TRUE, )
}
}
\keyword{ sft }
