\name{mic.test}
\alias{mic.test}
\title{
Test of the Mean Interaction Contrast
}
\description{
Performs either an Adjusted Rank Transform or ANOVA test for an interaction at the mean level.
}
\usage{
mic.test(HH, HL, LH, LL, method=c("art", "anova"))
}
\arguments{
  \item{HH}{Response times from the High--High condition.}
  \item{HL}{Response times from the High--Low condition.}
  \item{LH}{Response times from the Low--High condition.}
  \item{LL}{Response times from the Low--Low condition.}
  \item{method}{If "art", use the adjusted rank transform test.  If "anova" use ANOVA.}
}
\details{
  The mean interaction contrast (MIC) indicates the architecture of a process.  Serial processes result in MIC equal to zero.  Parallel-OR and Coactive process have a positive MIC.  Parallel-AND process have a negative MIC.  A test for a significant MIC can be done with a nonparametric adjusted rank transform test (described below) or an ANOVA.

  The Adjusted Rank Transform is a nonparametric test for an interaction between two discrete variables.  The test is carried out by first subtracting the mean effect of the salience level on each channel.  Suppose, 
\eqn{m_{H,\cdot} =} E[RT; Level of Channel 1 is Fast],
\eqn{m_{L,\cdot} =} E[RT; Level of Channel 1 is Slow],
\eqn{m_{\cdot, H} =} E[RT; Level of Channel 2 is Fast],
\eqn{m_{\cdot, L} =} E[RT; Level of Channel 2 is Slow].
Then for each response time from the fast--fast condition, \eqn{m_{H, \cdot}} and \eqn{m_{\cdot,H}} are subtracted.  Likewise, for each of the other conditions, the appropriate \eqn{m} is subtracted.  Next, each mean subtracted response time is replaced with its rank across all conditions (e.g., the fastest time of all conditions would be replaced with a 1).  In this implementation, tied response times are assigned using the average rank.  Finally, a standard ANOVA on the ranks is done on the ranks and the \eqn{p}-value of that test is returned.  This test was recommended by Sawilowsky (1990) based on a survey of a number of nonparametric tests for interactions.  He credits Reinach (1965) for first developing the test.
}
\value{
A list of class "htest" containing:
 \item{statistic}{The value of the MIC.}
 \item{p.value}{The \eqn{p}-value of the ART or ANOVA test.}
 \item{alternative}{A description of the alternative hypothesis.}
 \item{method}{A string indicating that the Houpt-Townsend statistic was used.}
 \item{data.name}{A string indicating which data were used for which input.}
}
\references{
Reinach, S.G. (1965).  A nonparametric analysis for a multiway classification with one element per cell.  \emph{South African Journal of Agricultural Science, 8}, 941--960.

Sawilowsky, S.S. (1990).  Nonparametric tests of interaction in experimental design.  \emph{Review of Educational Research, 60}, 91--126.

Houpt, J.W. & Townsend, J.T. (2010). The statistical properties of the survivor interaction contrast. \emph{Journal of Mathematical Psychology, 54}, 446-453.
}
\author{
Joe Houpt <joseph.houpt@wright.edu>
}
\examples{
T1.h <- rweibull(300, shape=2 , scale=400 )
T1.l <- rweibull(300, shape=2 , scale=800 )
T2.h <- rweibull(300, shape=2 , scale=400 )
T2.l <- rweibull(300, shape=2 , scale=800 )

Serial.hh <- T1.h + T2.h
Serial.hl <- T1.h + T2.l
Serial.lh <- T1.l + T2.h
Serial.ll <- T1.l + T2.l
mic.test(HH=Serial.hh, HL=Serial.hl, LH=Serial.lh, LL=Serial.ll)

Parallel.hh <- pmax(T1.h, T2.h)
Parallel.hl <- pmax(T1.h, T2.l)
Parallel.lh <- pmax(T1.l, T2.h)
Parallel.ll <- pmax(T1.l, T2.l)
mic.test(HH=Parallel.hh, HL=Parallel.hl, LH=Parallel.lh, LL=Parallel.ll, method="art")

}
\keyword{ sft }
