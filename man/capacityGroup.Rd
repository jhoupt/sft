\name{capacityGroup}
\alias{capacityGroup}
\title{
Capacity Analysis
}
\description{
Performs workload capacity analysis on each participant and each condition.  Plots each capacity coefficient individually and returns the results of the nonparametric null-hypothesis test for unlimited capacity independent parallel performance.
}
\usage{
capacityGroup(inData, acc.cutoff=.9, ratio=TRUE, OR=NULL, 
  stopping.rule=c("OR", "AND", "STST"), plotCt=TRUE, ...)
}
\arguments{
  \item{inData}{Data collected from a Double Factorial Paradigm experiment in standard form.}
  \item{acc.cutoff}{Minimum accuracy for each stimulus category used in calculating the capacity coefficient.}
  \item{OR}{Indicates whether to compare performance to an OR or AND processing baseline.  Provided for backwards compatibility for package version < 2. }
  \item{stopping.rule}{Indicates whether to use OR, AND or Single Target Self Terminating (STST) processing baseline to calculate individual capacity functions.}
  \item{ratio}{Indicates whether to return the standard ratio capacity coefficient or, if FALSE, the difference form.}
  \item{plotCt}{Indicates whether or not to generate plots of the capacity coefficients.}
  \item{...}{Arguments to be passed to plot function.}
}
\details{
For the details of the capacity coefficients, see \code{\link{capacity.or}}, \code{\link{capacity.and}} and \code{\link{capacity.stst}}.  If accuracy in any of the stimulus categories used to calculate a capacity coefficient falls below the cutoff, NA is retuned for that value in both the statistic and the Ct matrix.
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
Townsend, J.T. & Wenger, M.J. (2004). A theory of interactive parallel processing: New capacity measures and predictions for a response time inequality series. \emph{Psychological Review, 111}, 1003--1035.

Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Houpt, J.W. & Townsend, J.T. (2012). Statistical Measures for Workload Capacity Analysis.  \emph{Journal of Mathematical Psychology, 56}, 341-355.

Houpt, J.W., Blaha, L.M., McIntire, J.P., Havig, P.R. and Townsend, J.T. (2013). Systems Factorial Technology with R. \emph{Behavior Research Methods}.
}
\author{
Joe Houpt <joseph.houpt@utsa.edu>
}
\seealso{
\code{\link{capacity.and}}
\code{\link{capacity.or}}
\code{\link{capacity.stst}}
\code{\link{ucip.test}}
}
\examples{
\dontrun{
data(dots)
capacityGroup(subset(dots, Condition=="OR"), 
  stopping.rule="OR")
capacityGroup(subset(dots, Condition=="AND"), 
  stopping.rule="AND")
}
}
\keyword{ sft }
