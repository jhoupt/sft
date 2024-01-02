\name{capacity.id}
\alias{capacity.id}
\title{
Capacity Coefficient for Full Identification (ID) Exhaustive Processing
}
\description{
Calculates the capacity coefficient for exhaustive processing accounting for processing differences in no-responses to single targets.  The motivation for this version of the AND capacity coefficient is described in Howard et. al (2020).
}
\usage{
capacity.id(dual.target.rt, no.target.rt, single.target.rt, dual.target.cr, no.target.cr, single.target.crs, ratio=TRUE)
}
\arguments{
  \item{dual.target.rt}{Response times from trials on which all signals indicate targets.}
  \item{no.target.rt}{Response times from trials on which all signals are either absent or are non-targets.}
  \item{single.target.rts}{A list of arrays of response times from with each list containing response times for trials on which a distinct signal is the only signal present or that is indicating the target response.}
  \item{dual.target.cr}{A vector of indicators from trials on which all signals indicate targets indicating whether the participant correctly responded to the corresponding trial. If NULL, assumes all are correct.}
  \item{no.target.cr}{A vector of indicators from trials on which all signals are either absent or are non-targets indicating whether the participant correctly responded to the corresponding trial. If NULL, assumes all are correct.}
  \item{single.target.crs}{A list of vectors of indicators from with each list containing response times for trials on which a distinct signal is the only signal present or that is indicating the target response indicating whether the participant correctly responded to the corresponding trial. If NULL, assumes all are correct.}
  \item{ratio}{Indicates whether to return the standard ratio capacity coefficient or, if FALSE, the difference form.}
}
\details{
The identification-AND capacity coefficient compares performance on task to an unlimited-capacity, independent, parallel (UCIP) model using cumulative reverse hazard functions.  Suppose \eqn{K_i(t)}{Ki(t)} is the cumulative reverse hazard function for response times when process \eqn{i}{i} indicates a target-present response but all other signals imply a non-target response, \eqn{K_{\rm AND}(t)}{Kand(t)} is the cumulative reverse hazard function for response times when all target processes must completed together, and \eqn{K_{\rm NT}(t)}{Knt(t)} si teh cumulative reverse hazard function for response times when all non-target response processes are completed together.  Then the ID capacity coefficient is given by, 
\deqn{C_{\rm ID}(t)=\frac{\sum_i K_i(t)}{K_{\rm AND}(t) + K_{\rm ID}}.}{Cand(t)=sum[Ki(t)]/[Kand(t) + Knt(t)].}
The numerator is the estimated cumulative reverse hazard function for the UCIP model, based on the response times for each process in isolation and the denominator is the actual performance with a correction for the non-target response processes present in the single-target conditions. 

\eqn{C_{\rm ID}(t) <1}{Cid(t)<1} implies worse performance than the UCIP model.  This indicates that either there are limited processing resources, there is inhibition among the subprocesses, or the items are not processed in parallel (e.g., the items may be processed serially).

\eqn{C_{\rm ID}(t) >1}{Cid(t)>1} implies better performance than the UCIP model.  This indicates that either there are more processing resources available per process when there are more processes, that there is facilitation among the subprocesses, or the items are not processed in parallel (e.g., the items may be processed coactively).

The difference form of the capacity coefficient (returned if ratio=FALSE) is given by,
\deqn{C_{\rm ID}(t)=K_{\rm ID}(t) + K_{\rm NT}(t) - \sum_i K_i(t).}{Cid(t)=Kand(t) + Knt(t) - sum[Ki(t)].}  Negative values indicate worse than UCIP performance and positive values indicate better than UCIP performance.
}
\value{
  \item{Ct}{An object of class approxfun representing the estimated ID capacity coefficient.}
  \item{Var}{An object of class approxfun representing the variance of the estimated ID capacity coefficient.  Only returned if ratio=FALSE.}
  \item{Ctest}{A list with class "htest" that is returned from \code{\link{ucip.test}} and contains the statistic and p-value.}
}
\references{
Townsend, J.T. & Wenger, M.J. (2004). A theory of interactive parallel processing: New capacity measures and predictions for a response time inequality series. \emph{Psychological Review, 111}, 1003--1035.

Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Howard, Z. L., Garrett, P., Little, D. R., Townsend, J. T., & Eidels, A. (2021). A show about nothing: No-signal processes in systems factorial technology. Psychological Review, 128(1), 187-201. doi:https://doi.org/10.1037/rev0000256 

Houpt, J.W. & Townsend, J.T. (2012). Statistical Measures for Workload Capacity Analysis.  \emph{Journal of Mathematical Psychology, 56}, 341-355.

Houpt, J.W., Blaha, L.M., McIntire, J.P., Havig, P.R. and Townsend, J.T. (2013). Systems Factorial Technology with R. \emph{Behavior Research Methods}.
}
\author{
Joe Houpt <joseph.houpt@utsa.edu>
}
\seealso{
\code{\link{ucip.test}}
\code{\link{capacityGroup}}
\code{\link{capacity.and}}
\code{\link{estimateUCIPand}}
\code{\link{estimateNAK}}
\code{\link{approxfun}}
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
tvec <- sort(unique(c(RT.pa, RT.ap, RT.pp.limited, RT.pp.unlimited, RT.pp.super)))

cap.limited <- capacity.id(dual.target.rt=RT.pp.limited, no.target.rt=RT.nt, 
			   single.target.rts=list(RT.pa, RT.ap))
print(cap.limited$Ctest)
cap.unlimited <- capacity.id(dual.target.rt=RT.pp.unlimited, no.target.rt=RT.nt, 
			     single.target.rts=list(RT.pa, RT.ap))
cap.super <- capacity.id(dual.target.rt=RT.pp.super, no.target.rt=RT.nt, 
			 single.target.rts=list(RT.pa, RT.ap))

matplot(tvec, cbind(cap.limited$Ct(tvec), cap.unlimited$Ct(tvec), cap.super$Ct(tvec)),
  type='l', lty=1, ylim=c(0,3), col=2:4, main="Example Capacity Functions", xlab="Time", 
  ylab="C(t)")
abline(1,0)
legend('topright', c("Limited", "Unlimited", "Super"), lty=1, col=2:4)

}
\keyword{ sft }
