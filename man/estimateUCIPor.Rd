\name{estimateUCIPor}
\alias{estimateUCIPor}
\title{
UCIP Performance on OR Tasks
}
\description{
Estimates the cumulative hazard function of an unlimited capacity, independent, parallel process on an OR task.
}
\usage{
estimateUCIPor(RT, CR)
}
\arguments{
  \item{RT}{A list of arrays of response times.  Each list is used to estimate the response time distribution of a separate channel.}
  \item{CR}{A list of arrays of correct (1) or incorrect (0) indicators corresponding to each element of the list RT.}
}

\details{
This function concerns the processing time of an unlimited capacity, independent, parallel (UCIP) system.  This means that the completion time for each processing channel does not vary based on the presence of other processes.  Thus, the performance on tasks with a single process can be used to estimate performance of the UCIP model with multiple processes occurring. 

For example, in a two channel UCIP system the probability that no process has finished (OR processing) is the product of the probabilities of that each channel has not finished.  
\deqn{P(T_{ab} > t) = P(T_a > t)P(T_b > t)}{% 
      P(Tab > t) > P(Ta > t)P(Tb > t)}

We are interested in the cumulative hazard function, which is the natural log of the suvivor function (which is one minus the cumulative distribution function).  Because the log of a product is the sum of the logs, this gives us the following equality for the two channel OR process.
\deqn{H_{ab}(t) = H_a(t) + H_b(t)}{% 
      Hab(t) = Ha(t) + Hb(t)}

In general, the cumulative hazard function of a UCIP OR process is estimated by the sum of the cumulative hazard functions of each sub-process.
\deqn{H_{\rm UCIP}(t) = \sum_{i=1}^n H_i(t)}{% 
      Hucip(t) = sum( i = 1 .. n) Hi(t)}

The cumulative hazard functions of the sub-processes are estimated using the Nelson-Aalen estimator.  The Nelson-Aalen estimator is a Gaussian martingale, so the estimate of the UCIP performance is also a Gaussian martingale and the variance of the estimator can be estimated with the sum of variance estimates for each sub-process.
}
\value{
  \item{H }{A function of class "stepfun" that returns the Nelson-Aalen estimator of the cumulative hazard function of a UCIP model on a first-terminating (OR) task.}
  \item{Var }{A function of class "stepfun" that returns the estimated variance of the Nelson-Aalen estimator of the cumulative reverse hazard function of a UCIP model on a first-terminating (OR) task.}
}
\references{
Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Houpt, J.W. & Townsend, J.T. (2012). Statistical Measures for Workload Capacity Analysis.  \emph{Journal of Mathematical Psychology, 56}, 341-355.
}
\author{
Joe Houpt <joseph.houpt@utsa.edu>
}

\seealso{
\code{\link{estimateNAH}}
}
\examples{
# Channel completion times and accuracy
rt1 <- rexp(100, rate=.5)
cr1 <- runif(100) < .90
rt2 <- rexp(100, rate=.4)
cr2 <- runif(100) < .95
Hucip = estimateUCIPor(list(rt1, rt2), list(cr1, cr2))


# Plot the estimated UCIP cumulative hazard function
plot(Hucip$H, do.p=FALSE, 
  main="Estimated UCIP Cumulative Hazard Function\n
    X~min(X1,X2)    X1~Exp(.5)    X2~Exp(.4)", 
  xlab="X", ylab="H_UCIP(t)")
# Plot 95% Confidence intervals
times <- seq(0,10, length.out=100)
lines(times, Hucip$H(times) + sqrt(Hucip$Var(times))*qnorm(1-.05/2), lty=2)
lines(times, Hucip$H(times) - sqrt(Hucip$Var(times))*qnorm(1-.05/2), lty=2)
#Plot true UCIP cumulative hazard function
abline(0,.9, col='red')
}
\keyword{ sft }
