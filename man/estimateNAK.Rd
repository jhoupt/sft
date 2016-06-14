\name{estimateNAK}
\alias{estimateNAK}
\title{Neslon-Aalen Estimator of the Reverse Cumulative Hazard Function}
\description{
Computes the Nelson-Aalen estimator of a reverse cumulative hazard function.
}
\usage{
estimateNAK(RT, CR)
}
\arguments{
  \item{RT}{A vector of times at which an event occurs (e.g., a vector of response times).}
  \item{CR}{A vector of status indicators, 1=normal, 0=censored.  For response time data, this corresponds to 1=correct, 0=incorrect.}
}
\details{
The Nelson-Aalen estimator of the cumulative reverse hazard function is a step function with jumps at each event time.  The jump size is given by the number of events that have occurred up to and including the event.  If G(t) is the number events that have occurred up to and including t, then the N-A estimator of the cumulative reverse hazard function is given by:
\deqn{K(t) = -\sum_{s \in \{{\rm Event Times} > t\}} \frac{1}{G(s)}}{% 
      K(t) = -sum(s in {Event Times > t}) 1/G(s)}
}
\value{
  \item{H }{A function of class "stepfun" that returns the Nelson-Aalen estimator of the cumulative hazard function.}
  \item{Var }{A function of class "stepfun" that returns estimated variance of the Nelson-Aalen estimator of the cumulative hazard function.}
}
\references{
Aalen, O. O., Borgan, O., & Gjessing, H. K. (2008). \emph{Survival and event history analysis: A process point of view}. New York: Springer.

Houpt, J.W. & Townsend, J.T. (2012). Statistical Measures for Workload Capacity Analysis.  \emph{Journal of Mathematical Psychology, 56}, 341-355.
}
\author{
Joe Houpt <joseph.houpt@wright.edu>
}

\seealso{
\code{\link{estimateNAH}}
\code{\link{stepfun}}
}

\examples{
x <- rexp(50, rate=.5)
censoring <- runif(50) < .90
K.NA <- estimateNAK(x, censoring)

# Plot the estimated cumulative reverse hazard function
plot(K.NA$K, 
  main="Cumulative Reverse Hazard Function\n X ~ Exp(.5)    n=50", 
  xlab="X", ylab="K(x)")

# Plot 95% Confidence intervals
times <- seq(0,10, length.out=100)
lines(times, K.NA$K(times) + sqrt(K.NA$Var(times))*qnorm(1-.05/2), lty=2)
lines(times, K.NA$K(times) - sqrt(K.NA$Var(times))*qnorm(1-.05/2), lty=2)

# Plot the true cumulative reverse hazard function
lines(times, log(pexp(times, .5)), col='red')
}

\keyword{ survival }
\keyword{ sft }
