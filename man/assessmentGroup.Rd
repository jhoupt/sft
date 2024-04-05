\name{assessmentGroup}
\alias{assessmentGroup}
\title{
Assessment Functions
}
\description{
Calculates the specified assessment function for each participant and each condition.
}
\usage{
assessmentGroup(inData, stopping.rule=c("OR", "AND", "STST"), correct=c(TRUE, FALSE), 
  fast=c(TRUE, FALSE), detection=TRUE, plotAt=TRUE, ...)
}
\arguments{
  \item{inData}{Data collected from a Double Factorial Paradigm experiment in standard form.}
  \item{stopping.rule}{Indicates whether to use OR, AND, or single-target-self-terminating (STST) processing baseline to calculate individual assessment functions.}
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
  \item{overview}{A data frame indicating order of subject and condition for the assessment functions.}
  \item{At.fn}{Matrix with each row giving the values of the of the estimated assessment function for one participant in one condition for values of times.  The rows match the ordering of statistic.}
  \item{assessment}{A list with the returned values from the assessment function for each participant and each condition.}
  \item{times}{Times at which the assessment functions are calculated in At.fn matrix.}
}
\references{
Townsend, J.T. and Altieri, N. (2012).  An accuracy-response time capacity assessment function that measures performance against standard parallel predictions. \emph{Psychological Review, 3}, 500-516.

Donkin, C., Little, D. R., and Houpt, J. W. (2014). Assessing the speed-accuracy trade-off effect on the capacity of information processing. \emph{Journal of Experimental Psychology: Human Perception and Performance, 40(3),} 1183.
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
