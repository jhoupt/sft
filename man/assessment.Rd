\name{assessment}
\alias{assessment}
\title{
Workload Assessment Functions
}
\description{
Calculates the Workload Assessment Functions
}
\usage{
assessment(RT, CR, stopping.rule=c("OR","AND", "STST"), correct=c(TRUE, FALSE), 
  fast=c(TRUE, FALSE), detection=TRUE) 
}
\arguments{
  \item{RT}{A list of response time arrays.  The first array in the list is assumed to be the exhaustive condition.}
  \item{CR}{A list of correct/incorrect indicator arrays.  If NULL, assumes all are correct.}
  \item{stopping.rule}{Indicates whether to compare performance to an OR, AND, or STST processing baseline.}
  \item{correct}{Indicates whether to assess performance on correct trials.}
  \item{fast}{Indicates whether to use cumulative distribution functions or survivor functions to assess performance.}
  \item{detection}{Indicates whether to use a detection task baseline or a discrimination task baseline.}
}
\details{
The assessment functions are a generalization of the workload capacity functions that account for incorrect responses. Townsend & Altieri (2012) derived four different assessment functions each for AND and OR tasks to compare performance with two targets with the performance of an unlimited-capacity, independent, parallel (UCIP) model.  The correct assessment functions assess performance on correct trials and the incorrect assessment functions assess performance on the trials with incorrect responses.  The fast assessment functions use the cumulative distribution functions, similar to the AND capacity coefficient, and the slow assessment functions use the survivor functions, similar to the OR capacity coefficient.

In an OR task, the detection model assumes that the response will be correct if it is correct on either source, i.e., if either source is detected.  In discrimination OR tasks, a participant may respond based on whichever source finishes first.  Hence, the response will be incorrect if the first to finish is incorrect even if the second source would have been correct.  This results in a slightly different baseline for performance assessment.  See Donkin, Little and Houpt (2013) for details.
}
\value{
  \item{A}{An object of class stepfun representing the estimated assessment function.}
}
\references{
Townsend, J.T. and Altieri, N. (2012).  An accuracy-response time capacity assessment function that measures performance against standard parallel predictions. \emph{Psychological Review, 3}, 500-516.

Donkin, C, Little, D.R. and Houpt (2013). Assessing the effects of caution on the capacity of information processing.  \emph{Manuscript submitted for publication.}
}
\author{
Joe Houpt <joseph.houpt@utsa.edu>
}
\seealso{
\code{\link{capacity.or}}
\code{\link{capacity.and}}
\code{\link{stepfun}}
}
\examples{
c1c.12 <- rexp(10000, .015)
c1i.12 <- rexp(10000, .01)
c1c    <- rexp(10000, .015)
c1i    <- rexp(10000, .01)

c2c.12 <- rexp(10000, .014)
c2i.12 <- rexp(10000, .01)
c2c    <- rexp(10000, .014)
c2i    <- rexp(10000, .01)

RT.1 <- pmin(c1c, c1i)
CR.1 <- c1c < c1i
RT.2 <- pmin(c2c, c2i)
CR.2 <- c2c < c2i

c1Correct <- c1c.12 < c1i.12
c2Correct <- c2c.12 < c2i.12

# OR Detection
CR.12 <- c1Correct | c2Correct
RT.12 <- rep(NA, 10000)
RT.12[c1Correct & c2Correct] <- pmin(c1c.12, c2c.12)[c1Correct & c2Correct]
RT.12[c1Correct & !c2Correct] <- c1c.12[c1Correct & !c2Correct]
RT.12[!c1Correct & c2Correct] <- c2c.12[!c1Correct & c2Correct]
RT.12[!c1Correct & !c2Correct] <- pmax(c1i.12, c2i.12)[!c1Correct & !c2Correct]

RT <- list(RT.12, RT.1, RT.2)
CR <- list(CR.12, CR.1, CR.2)
a.or.cf <- assessment(RT, CR, stopping.rule="OR", correct=TRUE, fast=TRUE, detection=TRUE)
a.or.cs <- assessment(RT, CR, stopping.rule="OR", correct=TRUE, fast=FALSE, detection=TRUE)
a.or.if <- assessment(RT, CR, stopping.rule="OR", correct=FALSE, fast=TRUE, detection=TRUE)
a.or.is <- assessment(RT, CR, stopping.rule="OR", correct=FALSE, fast=FALSE, detection=TRUE)

par(mfrow=c(2,2))
plot(a.or.cf, ylim=c(0,2))
plot(a.or.cs, ylim=c(0,2))
plot(a.or.if, ylim=c(0,2))
plot(a.or.is, ylim=c(0,2))


# AND 
CR.12 <- c1Correct & c2Correct
RT.12 <- rep(NA, 10000)
RT.12[CR.12] <- pmax(c1c.12, c2c.12)[CR.12]
RT.12[c1Correct & !c2Correct] <- c2i.12[c1Correct & !c2Correct]
RT.12[!c1Correct & c2Correct] <- c1i.12[!c1Correct & c2Correct]
RT.12[!c1Correct & !c2Correct] <- pmin(c1i.12, c2i.12)[!c1Correct & !c2Correct]

RT <- list(RT.12, RT.1, RT.2)
CR <- list(CR.12, CR.1, CR.2)
a.and.cf <- assessment(RT, CR, stopping.rule="AND", correct=TRUE, fast=TRUE, detection=TRUE)
a.and.cs <- assessment(RT, CR, stopping.rule="AND", correct=TRUE, fast=FALSE, detection=TRUE)
a.and.if <- assessment(RT, CR, stopping.rule="AND", correct=FALSE, fast=TRUE, detection=TRUE)
a.and.is <- assessment(RT, CR, stopping.rule="AND", correct=FALSE, fast=FALSE, detection=TRUE)

par(mfrow=c(2,2))
plot(a.and.cf, ylim=c(0,2))
plot(a.and.cs, ylim=c(0,2))
plot(a.and.if, ylim=c(0,2))
plot(a.and.is, ylim=c(0,2))
}
\keyword{ sft }
