\name{capacity.or}
\alias{capacity.or}
\title{
Capacity Coefficient for First-Terminating (OR) Processing
}
\description{
Calculates the Capacity Coefficient for First-Terminating (OR) Processing
}
\usage{
capacity.or(RT, CR=NULL, ratio=TRUE)
}
\arguments{
  \item{RT}{A list of response time arrays.  The first array in the list is assumed to be the exhaustive condition.}
  \item{CR}{A list of correct/incorrect indicator arrays.  If NULL, assumes all are correct.}
  \item{ratio}{Indicates whether to return the standard ratio capacity coefficient or, if FALSE, the difference form.}
}
\details{
The OR capacity coefficient compares performance on task to an unlimited-capacity, independent, parallel (UCIP) model using cumulative hazard functions.  Suppose \eqn{H_i(t)}{Hi(t)} is the cumulative hazard function for response times when process \eqn{i}{i} is completed in isolation and \eqn{H_i(t)}{Hi(t)} is the cumulative hazard function for response times when all processes occur together and a response is made as soon as any of the processes finish.  Then the OR capacity coefficient is given by, 
\deqn{C_{\rm OR}(t)=\frac{H_{\rm OR}(t)}{\sum_i H_i(t)}.}{Cor(t)=sum[Hor(t)]/[Hi(t)].}
The denominator is the estimated cumulative hazard function for the UCIP model, based on the response times for each process in isolation and the numerator is the actual performance. 

\eqn{C_{\rm OR}(t) <1}{Cor(t)<1} implies worse performance than the UCIP model.  This indicates that either there are limited processing resources, there is inhibition among the subprocesses, or the items are not processed in parallel (e.g., the items may be processed serially).

\eqn{C_{\rm OR}(t) >1}{Cor(t)>1} implies better performance than the UCIP model.  This indicates that either there are more processing resources available per process when there are more processes, that there is facilitation among the subprocesses, or the items are not processed in parallel (e.g., the items may be processed coactively).

The difference form of the capacity coefficient is given by,
\deqn{C_{\rm OR}(t)=H_{\rm OR}(t) - \sum_i H_i(t).}{Cor(t)=Hor(t) - sum[Hi(t)].}  Negative values indicate worse than UCIP performance and positive values indicate better than UCIP performance.
}
\value{
  \item{Ct}{An object of class approxfun representing the OR capacity coefficient.}
  \item{Var}{An object of class approxfun representing the variance of the estimated OR capacity coefficient.  Only returned if ratio=FALSE.}
  \item{Ctest}{A list with class "htest" that is returned from \code{\link{ucip.test}} and contains the statistic and p-value.}
}
\references{
Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Houpt, J.W. & Townsend, J.T. (2012). Statistical Measures for Workload Capacity Analysis.  \emph{Journal of Mathematical Psychology, 56}, 341-355.

Houpt, J.W., Blaha, L.M., McIntire, J.P., Havig, P.R. and Townsend, J.T. (2013). Systems Factorial Technology with R. \emph{Behavior Research Methods}.
}
\author{
Joe Houpt <joseph.houpt@wright.edu>
}
\seealso{
\code{\link{ucip.test}}
\code{\link{capacityGroup}}
\code{\link{capacity.and}}
\code{\link{estimateUCIPor}}
\code{\link{estimateNAH}}
\code{\link{approxfun}}
}
\examples{
rate1 <- .35
rate2 <- .3
RT.pa <- rexp(100, rate1)
RT.ap <- rexp(100, rate2)
RT.pp.limited <- pmin( rexp(100, .5*rate1), rexp(100, .5*rate2))
RT.pp.unlimited <- pmin( rexp(100, rate1), rexp(100, rate2))
RT.pp.super <- pmin( rexp(100, 2*rate1), rexp(100, 2*rate2))
tvec <- sort(unique(c(RT.pa, RT.ap, RT.pp.limited, RT.pp.unlimited, RT.pp.super)))

cap.limited <- capacity.or(RT=list(RT.pp.limited, RT.pa, RT.ap))
print(cap.limited$Ctest)
cap.unlimited <- capacity.or(RT=list(RT.pp.unlimited, RT.pa, RT.ap))
cap.super <- capacity.or(list(RT=RT.pp.super, RT.pa, RT.ap))

matplot(tvec, cbind(cap.limited$Ct(tvec), cap.unlimited$Ct(tvec), cap.super$Ct(tvec)),
  type='l', lty=1, ylim=c(0,3), col=2:4, main="Example Capacity Functions", xlab="Time", 
  ylab="C(t)")
abline(1,0)
legend('topright', c("Limited", "Unlimited", "Super"), lty=1, col=2:4)
}
\keyword{ sft }
