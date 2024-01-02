\name{sic.test}
\alias{sic.test}
\title{
Statistical test of the SIC.
}
\description{
Function to test for statistical significance of the positive and negative parts of a SIC.
}
    
\usage{
sic.test(HH, HL, LH, LL, method="ks")
}
\arguments{
  \item{HH}{Response times from the High--High condition.}
  \item{HL}{Response times from the High--Low condition.}
  \item{LH}{Response times from the Low--High condition.}
  \item{LL}{Response times from the Low--Low condition.}
  \item{method}{Which type of hypothesis test to use for SIC form.}
}
\details{
     SIC(t) = (S_LL - S_LH) - (S_HL - S_HH)

This function performs a statistical analysis to determine whether the positive and negative parts of the SIC are significantly different from zero. Currently the only statistical test is based on the generalization of the two-sample Kolmogorov-Smirnov test described in Houpt & Townsend, 2010.  This test performs two separate null-hypothesis tests:  One test for whether the largest positive value of the SIC is significantly different from zero and one test for whether the largest negative value is significantly different from zero.
}
\value{
 \item{positive}{A list of class "htest" containing the statistic and \eqn{p}-value along with descriptions of the alternative hypothesis, method and data names for the test of a significant positive portion of the SIC.}
 \item{negative}{A list of class "htest" containing the statistic and \eqn{p}-value along with descriptions of the alternative hypothesis, method and data names for the test of a significant negative portion of the SIC.}
}
\references{
Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Houpt, J.W. & Townsend, J.T. (2010). The statistical properties of the survivor interaction contrast. \emph{Journal of Mathematical Psychology, 54}, 446-453.
}
\author{
Joe Houpt <joseph.houpt@utsa.edu>
}
\seealso{
\code{\link{stepfun}}
\code{\link{sicGroup}}
\code{\link{sic}}
\code{\link{mic.test}}
}
\examples{
T1.h <- rexp(50, .2)
T1.l <- rexp(50, .1)
T2.h <- rexp(50, .21)
T2.l <- rexp(50, .11)

SerialAND.hh <- T1.h + T2.h
SerialAND.hl <- T1.h + T2.l
SerialAND.lh <- T1.l + T2.h
SerialAND.ll <- T1.l + T2.l
sic.test(HH=SerialAND.hh, HL=SerialAND.hl, LH=SerialAND.lh, LL=SerialAND.ll)

p1 <- runif(200) < .3
SerialOR.hh <- p1[1:50]    * T1.h + (1-p1[1:50]   )*T2.h
SerialOR.hl <- p1[51:100]  * T1.h + (1-p1[51:100] )*T2.l
SerialOR.lh <- p1[101:150] * T1.l + (1-p1[101:150])*T2.h
SerialOR.ll <- p1[151:200] * T1.l + (1-p1[151:200])*T2.l
sic.test(HH=SerialOR.hh, HL=SerialOR.hl, LH=SerialOR.lh, LL=SerialOR.ll)

ParallelAND.hh <- pmax(T1.h, T2.h)
ParallelAND.hl <- pmax(T1.h, T2.l)
ParallelAND.lh <- pmax(T1.l, T2.h)
ParallelAND.ll <- pmax(T1.l, T2.l)
sic.test(HH=ParallelAND.hh, HL=ParallelAND.hl, LH=ParallelAND.lh, LL=ParallelAND.ll)

ParallelOR.hh <- pmin(T1.h, T2.h)
ParallelOR.hl <- pmin(T1.h, T2.l)
ParallelOR.lh <- pmin(T1.l, T2.h)
ParallelOR.ll <- pmin(T1.l, T2.l)
sic.test(HH=ParallelOR.hh, HL=ParallelOR.hl, LH=ParallelOR.lh, LL=ParallelOR.ll)
}
\keyword{ sft }
