\name{capacity.stst}
\alias{capacity.stst}
\title{
Capacity Coefficient for Single-Target Self-Terminating (STST) Processing
}
\description{
Calculates the Capacity Coefficient for Single-Target Self-Terminating (STST) Processing
}
\usage{
capacity.stst(RT, CR=NULL, ratio=TRUE)
}
\arguments{
  \item{RT}{A list of response time arrays.  The first array in the list is assumed to be the single-target among N distractors condition.}
  \item{CR}{A list of correct/incorrect indicator arrays.  If NULL, assumes all are correct.}
  \item{ratio}{Indicates whether to return the standard ratio capacity coefficient or, if FALSE, the difference form.}
}
\details{
The STST capacity coefficient compares performance on task to an unlimited-capacity, independent, parallel (UCIP) model using cumulative reverse hazard functions.  Suppose \eqn{K_{i,1}(t)}{Ki,1(t)} is the cumulative reverse hazard function for response times when single-target process \eqn{i}{i} is completed in isolation and \eqn{K_{i,n}(t)}{Ki,n(t)} is the cumulative reverse hazard function for response times when the single-target \eqn{i}{i} is processed among \eqn{n}{n} other processes, all completed together.  Then the STST capacity coefficient is given by, 
\deqn{C_{\rm STST}(t)=\frac{K_{i,1}(t)}{K_{i,n}(t)}.}{Cstst(t)=[Ki,1(t)]/[Ki,n(t)].}
The numerator is the estimated cumulative reverse hazard function for the UCIP model, based on the response times for the \eqn{i}{i} process in isolation and the denominator is the actual performance on the \eqn{i}{i} process among \eqn{n}{n} distractors or other active channels. 

\eqn{C_{\rm STST}(t) <1}{Cstst(t)<1} implies worse performance than the UCIP model.  This indicates that either there are limited processing resources, there is inhibition among the subprocesses, or the items are not processed in parallel (e.g., the items may be processed serially).

\eqn{C_{\rm STST}(t) >1}{Cstst(t)>1} implies better performance than the UCIP model.  This indicates that either there are more processing resources available per process when there are more processes, that there is facilitation among the subprocesses, or the items are not processed in parallel (e.g., the items may be processed coactively).

The difference form of the capacity coefficient (returned if ratio=FALSE) is given by,
\deqn{C_{\rm STST}(t)=K_{i,n}(t) - K_{i,1}(t).}{Cstst(t)=K1,n(t) - Ki,1(t).}  Negative values indicate worse than UCIP performance and positive values indicate better than UCIP performance.
}
\value{
  \item{Ct}{An object of class approxfun representing the estimated STST capacity coefficient.}
  \item{Var}{An object of class approxfun representing the variance of the estimated STST capacity coefficient.  Only returned if ratio=FALSE.}
  \item{Ctest}{A list with class "htest" that is returned from \code{\link{ucip.test}} and contains the statistic and p-value.}
}
\references{
Blaha, L.M. & Townsend, J.T. (under review). On the capacity of single-target self-terminating processes.

Houpt, J.W. & Townsend, J.T. (2012). Statistical measures for workload capacity analysis.  \emph{Journal of Mathematical Psychology, 56}, 341-355.

Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Townsend, J.T. & Wenger, M.J. (2004). A theory of interactive parallel processing: New capacity measures and predictions for a response time inequality series. \emph{Psychological Review, 111}, 1003--1035.
}

\author{
Leslie Blaha <leslie.blaha@us.af.mil>

Joe Houpt <joseph.houpt@utsa.edu>
}
\seealso{
\code{\link{ucip.test}}
\code{\link{capacityGroup}}
\code{\link{capacity.or}}
\code{\link{capacity.and}}
\code{\link{estimateNAK}}
\code{\link{approxfun}}
}
\examples{
rate1 <- .35
RT.pa <- rexp(100, rate1)
RT.pp.limited <- rexp(100, .5*rate1)
RT.pp.unlimited <- rexp(100, rate1)
RT.pp.super <- rexp(100, 2*rate1)
tvec <- sort(unique(c(RT.pa, RT.pp.limited, RT.pp.unlimited, RT.pp.super)))

cap.limited <- capacity.stst(RT=list(RT.pp.limited, RT.pa))
print(cap.limited$Ctest)
cap.unlimited <- capacity.stst(RT=list(RT.pp.unlimited, RT.pa))
cap.super <- capacity.stst(RT=list(RT.pp.super, RT.pa))

matplot(tvec, cbind(cap.limited$Ct(tvec), cap.unlimited$Ct(tvec), cap.super$Ct(tvec)),
  type='l', lty=1, ylim=c(0,5), col=2:4, main="Example Capacity Functions", xlab="Time", 
  ylab="C(t)")
abline(1,0)
legend('topright', c("Limited", "Unlimited", "Super"), lty=1, col=2:4, bty="n")

}
\keyword{ sft }
