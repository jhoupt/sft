\name{siDominance}
\alias{siDominance}
\title{
Dominance Test for Selective Influence
}
\description{
Function to test for the survivor function ordering predicted by the selective influence of the salience manipulation.
}
    
\usage{
siDominance(HH, HL, LH, LL, method="ks")
}
\arguments{
  \item{HH}{Response times from the High--High condition.}
  \item{HL}{Response times from the High--Low condition.}
  \item{LH}{Response times from the Low--High condition.}
  \item{LL}{Response times from the Low--Low condition.}
  \item{method}{Which type of hypothesis test to use for testing stochastic dominance relations, either as series of KS tests ("ks") or the dominance test based on Dirichlet process priors ("dp"). DP not yet implemented.}
}
\details{
For an SIC function to distinguish among the processing types, the salience manipulation on each channel must selectively influence its respective channel (although see Eidels, Houpt, Altieri, Pei & Townsend, 2010 for SIC prediction from interactive parallel models).  Although the selective influence assumption cannot be directly tested, one implication is that the distribution the HH response times stochastically dominates the HL and LH distributions which each in turn stochastically dominate the LL response time distribution.  This implication is automatically tested in this function.  The KS dominance test uses eight two-sample Kolmogorov-Smirnov tests:  HH < HL, HH < LH, HL < LL, LH < LL should be significant while HH > HL, HH > LH, HL > LL, LH > LL should not.  The DP uses four tests to determine which relation has the highest Bayes factor assuming a Dirichlet process prior for each of (HH, HL), (HH, LH), (HL, LL) and (LH, LL).  See Heathcote, Brown, Wagenmakers & Eidels, 2010, for more details.

}
\value{
 A data frame with the first column indicating which ordering was tested, the second column indicating the test statistic and the third indicating the p-value for that value of the statistic.
}
\references{
Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Houpt, J.W. & Townsend, J.T. (2010). The statistical properties of the survivor interaction contrast. \emph{Journal of Mathematical Psychology, 54}, 446-453.

Dzhafarov, E.N., Schweickert, R., & Sung, K. (2004). Mental architectures with selectively influenced but stochastically interdependent components. \emph{Journal of Mathematical Psychology}, 48, 51-64.
}
\author{
Joe Houpt <joseph.houpt@wright.edu>
}
\seealso{
\code{\link{ks.test}}
\code{\link{sic}}
\code{\link{sicGroup}}
\code{\link{mic.test}}
}
\examples{
T1.h <- rexp(50, .2)
T1.l <- rexp(50, .1)
T2.h <- rexp(50, .21)
T2.l <- rexp(50, .11)

HH <- T1.h + T2.h
HL <- T1.h + T2.l
LH <- T1.l + T2.h
LL <- T1.l + T2.l
siDominance(HH, HL, LH, LL)
}
\keyword{ sft }
