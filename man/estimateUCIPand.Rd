\name{estimateUCIPand}
\alias{estimateUCIPand}
\title{
UCIP Performance on AND Tasks
}
\description{
Estimates the reverse cumulative hazard function of an unlimited capacity, independent, parallel process on an AND task.
}
\usage{
estimateUCIPand(RT, CR)
}
\arguments{
  \item{RT}{A list of arrays of response times.  Each list is used to estimate the response time distribution of a separate channel.}
  \item{CR}{A list of arrays of correct (1) or incorrect (0) indicators corresponding to each element of the list RT.}
}
\details{
This function concerns the processing time of an unlimited capacity, independent, parallel (UCIP) system.  This means that the completion time for each processing channel does not vary based on the presence of other processes.  Thus, the performance on tasks with a single process can be used to estimate performance of the UCIP model with multiple processes occurring. 

For example, in a two channel UCIP system the probability that both processes have finished (AND processing) is the product of the probabilities of that each channel has finished.  
\deqn{P(T_{ab} \le t) = P(T_a \le t)P(T_b \le t)}{% 
      P(Tab <= t) = P(Ta <= t)P(Tb <= t)}

We are interested in the cumulative reverse hazard function, which is the natural log of the cumulative distribution function.  Because the log of a product is the sum of the logs, this gives us the following equality for the two channel AND process. 
\deqn{K_{ab}(t) = K_a(t) + K_b(t)}{% 
      Kab(t) = Ka(t) + Kb(t)}

In general, the cumulative reverse hazard function of a UCIP AND process is estimated by the sum of the cumulative reverse hazard functions of each sub-process.
\deqn{K_{\rm UCIP}(t) = \sum_{i=1}^n K_i(t)}{% 
      Kucip(t) = sum( i = 1 .. n) Ki(t)}

The cumulative reverse hazard functions of the sub-processes are estimated using the Nelson-Aalen estimator.  The Nelson-Aalen estimator is a Gaussian martingale, so the estimate of the UCIP performance is also a Gaussian martingale and the variance of the estimator can be estimated with the sum of variance estimates for each sub-process.
}

\value{
  \item{K }{A function of class "stepfun" that returns the Nelson-Aalen estimator of the cumulative reverse hazard function of a UCIP model on an exhaustive (AND) task.}
  \item{Var }{A function of class "stepfun" that returns the estimated variance of the Nelson-Aalen estimator of the cumulative reverse hazard function of a UCIP model on an exhaustive (AND) task.}
}
\references{
Townsend, J.T. & Wenger, M.J. (2004). A theory of interactive parallel processing: New capacity measures and predictions for a response time inequality series. \emph{Psychological Review, 111}, 1003-1035.

Houpt, J.W. & Townsend, J.T. (2012). Statistical Measures for Workload Capacity Analysis.  \emph{Journal of Mathematical Psychology, 56}, 341-355.
}
\author{
Joe Houpt <joseph.houpt@utsa.edu>
}

\seealso{
\code{\link{estimateNAK}}
}
\examples{
# Channel completion times and accuracy
rt1 <- rexp(100, rate=.5)
cr1 <- runif(100) < .90
rt2 <- rexp(100, rate=.4)
cr2 <- runif(100) < .95
Kucip = estimateUCIPand(list(rt1, rt2), list(cr1, cr2))


# Plot the estimated UCIP cumulative reverse hazard function
plot(Kucip$K, do.p=FALSE, 
  main="Estimated UCIP Cumulative Reverse Hazard Function\n
    X~max(X1,X2)    X1~Exp(.5)    X2~Exp(.4)", 
  xlab="X", ylab="K_UCIP(x)")
# Plot 95% Confidence intervals
times <- seq(0,10, length.out=100)
lines(times, Kucip$K(times) + sqrt(Kucip$Var(times))*qnorm(1-.05/2), lty=2)
lines(times, Kucip$K(times) - sqrt(Kucip$Var(times))*qnorm(1-.05/2), lty=2)
# Plot true UCIP cumulative reverse hazard function
lines(times[-1], log(pexp(times[-1], .5)) + log(pexp(times[-1], .4)), col='red')

}
\keyword{ sft }
