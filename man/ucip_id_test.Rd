\name{ucip.id.test}
\alias{ucip.id.test}
\title{
A Statistical Test for Super or Limited Capacity
}
\description{
A nonparametric test for capacity values significantly different than those predicted by the estimated unlimited capacity, independent parallel model.
}
\usage{
ucip.id.test(dt.rt, nt.rt, st.rts, dt.cr=NULL, nt.cr=NULL, st.crs=NULL)
}
\arguments{
  \item{dt.rt}{Response times from trials on which all signals indicate targets.}
  \item{nt.rt}{Response times from trials on which all signals are either absent or are non-targets.}
  \item{st.rts}{A list of arrays of response times from with each list containing response times for trials on which a distinct signal is the only signal present or that is indicating the target response.}
  \item{dt.cr}{A vector of indicators from trials on which all signals indicate targets indicating whether the participant correctly responded to the corresponding trial. If NULL, assumes all are correct.}
  \item{nt.cr}{A vector of indicators from trials on which all signals are either absent or are non-targets indicating whether the participant correctly responded to the corresponding trial. If NULL, assumes all are correct.}
  \item{st.crs}{A list of vectors of indicators from with each list containing response times for trials on which a distinct signal is the only signal present or that is indicating the target response indicating whether the participant correctly responded to the corresponding trial. If NULL, assumes all are correct.}
}
\details{
The test is based on the Nelson-Aalen formulation of the log-rank test.  The function takes a weighted difference between estimated unlimited capacity, independent parallel performance, based on a participants single source response times, and the participants true performance when all sources are present.  Use this statistic with caution. It has not been thorougly tested nor subjected to peer review. 
}
\value{
A list of class "htest" containing:
  \item{statistic}{Z-score of a null-hypothesis test for UCIP performance.}
  \item{p.val}{p-value of a two-tailed null-hypothesis test for UCIP performance.}
 \item{alternative}{A description of the alternative hypothesis.}
 \item{method}{A string indicating whether the ART or ANOVA was used.}
}
\references{
Houpt, J.W. & Townsend, J.T. (2012). Statistical Measures for Workload Capacity Analysis.  \emph{Journal of Mathematical Psychology, 56}, 341-355.

Howard, Z. L., Garrett, P., Little, D. R., Townsend, J. T., & Eidels, A. (2021). A show about nothing: No-signal processes in systems factorial technology. Psychological Review, 128(1), 187-201. doi:https://doi.org/10.1037/rev0000256 
}
\author{
Joe Houpt <joseph.houpt@utsa.edu>
}
\seealso{
\code{\link{capacity.or}}
\code{\link{capacity.and}}
\code{\link{capacity.id}}
\code{\link{estimateUCIPor}}
\code{\link{estimateUCIPand}}
\code{\link{estimateNAH}}
\code{\link{estimateNAK}}
}
\examples{
rate1p <- .35
rate1a <- .25
rate2p <- .35
rate2a <- .25
RT.pa <- pmax(rexp(100, rate1p), rexp(100, rate2a))
RT.ap <- pmax(rexp(100, rate2p), rexp(100, rate1a))
RT.nt <- pmax(rexp(100, rate1a), rexp(100, rate2a))
RT.pp.limited <- pmax( rexp(100, .5*rate1p), rexp(100, .5*rate2p))
RT.pp.unlimited <- pmax( rexp(100, rate1p), rexp(100, rate2p))
RT.pp.super <- pmax( rexp(100, 2*rate1p), rexp(100, 2*rate2p))

z.limited   <- ucip.id.test(dt.rt=RT.pp.limited,   nt.rt=RT.nt, st.rts=list(RT.pa, RT.ap)) 
z.unlimited <- ucip.id.test(dt.rt=RT.pp.unlimited, nt.rt=RT.nt, st.rts=list(RT.pa, RT.ap)) 
z.super     <- ucip.id.test(dt.rt=RT.pp.super,     nt.rt=RT.nt, st.rts=list(RT.pa, RT.ap)) 
}
\keyword{ ~sft }
