
\name{estimateNAH}
\alias{estimateNAH}
\title{Neslon-Aalen Estimator of the Cumulative Hazard Function}
\description{
Computes the Nelson-Aalen estimator of a cumulative hazard function.
}
\usage{
estimateNAH(RT, CR)
}
\arguments{
  \item{RT}{A vector of times at which an event occurs (e.g., a vector of response times).}
  \item{CR}{A vector of status indicators, 1=normal, 0=censored.  For response time data, this corresponds to 1=correct, 0=incorrect.}
}
\details{
The Nelson-Aalen estimator of the cumulative hazard function is a step function with jumps at each event time.  The jump size is given by the number at risk up until immediately before the event.  If Y(t) is the number at risk immediately before t, then the N-A estimator is given by:
\deqn{H(t) = \sum_{s \in \{{\rm Event Times} < t\}} \frac{1}{Y(s)}}{% 
      H(t) = sum(s in {Event Times < t}) 1/Y(s)}

}
\value{
  \item{H }{A function of class "stepfun" that returns the Nelson-Aalen estimator of the cumulative hazard function.}
  \item{Var }{A function of class "stepfun" that returns the estimated variance of the Nelson-Aalen estimator of the cumulative hazard function.}
}
\references{
Aalen, O. O., Borgan, O., & Gjessing, H. K. (2008). \emph{Survival and event history analysis: A process point of view}. New York: Springer.
}
\author{
Joe Houpt <joseph.houpt@wright.edu>
}

\seealso{
\code{\link{estimateNAK}}
\code{\link{stepfun}}
}

\examples{
x <- rexp(50, rate=.5)
censoring <- runif(50) < .90
H.NA <- estimateNAH(x, censoring)

# Plot the estimated cumulative hazard function
plot(H.NA$H, 
  main="Cumulative Hazard Function\n X ~ Exp(.5)    n=50", 
  xlab="X", ylab="H(x)")

# Plot 95% Confidence intervals
times <- seq(0,10, length.out=100)
lines(times, H.NA$H(times) + sqrt(H.NA$Var(times))*qnorm(1-.05/2), lty=2)
lines(times, H.NA$H(times) - sqrt(H.NA$Var(times))*qnorm(1-.05/2), lty=2)

# Plot the true cumulative hazard function
abline(0,.5, col='red')
}

\keyword{ survival }
\keyword{ sft }
