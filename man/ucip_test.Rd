\name{ucip.test}
\alias{ucip.test}
\title{
A Statistical Test for Super or Limited Capacity
}
\description{
A nonparametric test for capacity values significantly different than those predicted by the estimated unlimited capacity, independent parallel model.
}
\usage{
ucip.test(RT, CR=NULL, OR=NULL, stopping.rule=c("OR","AND","STST"))
}
\arguments{
  \item{RT}{A list of response time arrays.  The first array in the list is assumed to be the exhaustive condition.}
  \item{CR}{A list of correct/incorrect indicator arrays.  If NULL, assumes all are correct.}
  \item{OR}{Indicates whether to compare performance to an OR or AND processing baseline.  Provided for backwards compatibility for package version < 2. }
  \item{stopping.rule}{Indicates whether to compare performance to an OR, AND or Single Target Self Terminating (STST) processing baseline.}
}
\details{
The test is based on the Nelson-Aalen formulation of the log-rank test.  The function takes a weighted difference between estimated unlimited capacity, independent parallel performance, based on a participants single source response times, and the participants true performance when all sources are present.
}
\value{
A list of class "htest" containing:
  \item{statistic}{Z-score of a null-hypothesis test for UCIP performance.}
  \item{p.val}{p-value of a two-tailed null-hypothesis test for UCIP performance.}
 \item{alternative}{A description of the alternative hypothesis.}
 \item{method}{A string indicating whether the ART or ANOVA was used.}
 \item{data.name}{A string indicating which data were used for which input.}
}
\references{
Houpt, J.W. & Townsend, J.T. (2012). Statistical Measures for Workload Capacity Analysis.  \emph{Journal of Mathematical Psychology, 56}, 341-355.
}
\author{
Joe Houpt <joseph.houpt@utsa.edu>
}
\seealso{
\code{\link{capacity.or}}
\code{\link{capacity.and}}
\code{\link{estimateUCIPor}}
\code{\link{estimateUCIPand}}
\code{\link{estimateNAH}}
\code{\link{estimateNAK}}
}
\examples{
rate1 <- .35
rate2 <- .3
RT.pa <- rexp(100, rate1)
RT.ap <- rexp(100, rate2)

CR.pa <- runif(100) < .98
CR.ap <- runif(100) < .98
CR.pp <- runif(100) < .96
CRlist <- list(CR.pp, CR.pa, CR.ap)

#  OR Processing
RT.pp.limited <- pmin( rexp(100, .5*rate1), rexp(100, .5*rate2))
RT.pp.unlimited <- pmin( rexp(100, rate1), rexp(100, rate2))
RT.pp.super <- pmin( rexp(100, 2*rate1), rexp(100, 2*rate2))
z.limited   <- ucip.test(RT=list(RT.pp.limited, RT.pa, RT.ap), CR=CRlist, stopping.rule="OR")
z.unlimited <- ucip.test(RT=list(RT.pp.unlimited, RT.pa, RT.ap), CR=CRlist, stopping.rule="OR")
z.super     <- ucip.test(RT=list(RT.pp.super, RT.pa, RT.ap), CR=CRlist, stopping.rule="OR")

#  AND Processing
RT.pp.limited <- pmax( rexp(100, .5*rate1), rexp(100, .5*rate2))
RT.pp.unlimited <- pmax( rexp(100, rate1), rexp(100, rate2))
RT.pp.super <- pmax( rexp(100, 2*rate1), rexp(100, 2*rate2))
z.limited   <- ucip.test(RT=list(RT.pp.limited, RT.pa, RT.ap), CR=CRlist, stopping.rule="AND")
z.unlimited <- ucip.test(RT=list(RT.pp.unlimited, RT.pa, RT.ap), CR=CRlist, stopping.rule="AND")
z.super     <- ucip.test(RT=list(RT.pp.super, RT.pa, RT.ap), CR=CRlist, stopping.rule="AND")
}
\keyword{ ~sft }
