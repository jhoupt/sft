\name{sic}
\alias{sic}
\title{
Calculate the Survivor Interaction Contrast
}
\description{
Function to calculate survivor interaction contrast and associated measures.
}
    
\usage{
sic(HH, HL, LH, LL, domtest="ks", sictest="ks", mictest=c("art", "anova"))
}
\arguments{
  \item{HH}{Response times from the High--High condition.}
  \item{HL}{Response times from the High--Low condition.}
  \item{LH}{Response times from the Low--High condition.}
  \item{LL}{Response times from the Low--Low condition.}
  \item{sictest}{Which type of hypothesis test to use for SIC form.}
  \item{domtest}{Which type of hypothesis test to use for testing stochastic dominance relations, either as series of KS tests ("ks") or the dominance test based on Dirichlet process priors ("dp"). DP not yet implemented.}
  \item{mictest}{Which type of hypothesis test to use for the MIC, either adjusted rank transform or ANOVA. }
}
\details{
     SIC(t) = (S_LL - S_LH) - (S_HL - S_HH)

This function calculates the Survivor Interaction Contrast (SIC; Townsend & Nozawa, 1995).  The SIC indicates the architecture and stopping-rule of the underlying information processing system.  An entirely positive SIC indicates parallel first-terminating processing.  An entirely negative SIC indicates parallel exhaustive processing.  An SIC that is always zero indicates serial first-terminating processing.  An SIC that is first positive then negative indicates either serial exhaustive or coactive processing.  To distinguish between these two possibilities, an additional test of the mean interaction contrast (MIC) is used; coactive processing leads to a positive MIC while serial processing leads to an MIC of zero.

For the SIC function to distinguish among the processing types, the salience manipulation on each channel must selectively influence its respective channel (although see Eidels, Houpt, Altieri, Pei & Townsend, 2010 for SIC prediction from interactive parallel models).  Although the selective influence assumption cannot be directly tested, one implication is that the distribution the HH response times stochastically dominates the HL and LH distributions which each in turn stochastically dominate the LL response time distribution.  This implication is automatically tested in this function.  The KS dominance test uses eight two-sample Kolmogorov-Smirnov tests:  HH < HL, HH < LH, HL < LL, LH < LL should be significant while HH > HL, HH > LH, HL > LL, LH > LL should not.  The DP uses four tests to determine which relation has the highest Bayes factor assuming a Dirichlet process prior for each of (HH, HL), (HH, LH), (HL, LL) and (LH, LL).  See Heathcote, Brown, Wagenmakers & Eidels, 2010, for more details.

This function also performs a statistical analysis to determine whether the positive and negative parts of the SIC are significantly different from zero. Currently the only statistical test is based on the generalization of the two-sample Kolmogorov-Smirnov test described in Houpt & Townsend, 2010.  This test performs two separate null-hypothesis tests:  One test for whether the largest positive value of the SIC is significantly different from zero and one test for whether the largest negative value is significantly different from zero.
}
\value{
 \item{SIC}{An object of class stepfun representing the SIC.}
 \item{Dominance}{A data frame with the first column indicating which ordering was tested, the second column indicating the test statistic and the third indicating the \eqn{p}-value for that value of the statistic.}
 \item{Dvals}{A Matrix containing the values of the test statistic and the associated \eqn{p}-values.}
 \item{MIC}{Results of an adjusted rank transform test of the mean interaction contrast.}
 \item{N}{The scaling factor used for the KS test of the SIC form.}
}
\references{
Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Houpt, J.W. & Townsend, J.T. (2010). The statistical properties of the survivor interaction contrast. \emph{Journal of Mathematical Psychology, 54}, 446-453.

Houpt, J.W., Blaha, L.M., McIntire, J.P., Havig, P.R. and Townsend, J.T. (2013). Systems Factorial Technology with R. \emph{Behavior Research Methods}.
}
\author{
Joe Houpt <joseph.houpt@utsa.edu>
}
\seealso{
\code{\link{stepfun}}
\code{\link{sicGroup}}
\code{\link{mic.test}}
\code{\link{sic.test}}
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
SerialAND.sic <- sic(HH=SerialAND.hh, HL=SerialAND.hl, LH=SerialAND.lh, 
  LL=SerialAND.ll)
print(SerialAND.sic$Dvals)
plot(SerialAND.sic$SIC, do.p=FALSE, ylim=c(-1,1))

p1 <- runif(200) < .3
SerialOR.hh <- p1[1:50]    * T1.h + (1-p1[1:50]   )*T2.h
SerialOR.hl <- p1[51:100]  * T1.h + (1-p1[51:100] )*T2.l
SerialOR.lh <- p1[101:150] * T1.l + (1-p1[101:150])*T2.h
SerialOR.ll <- p1[151:200] * T1.l + (1-p1[151:200])*T2.l
SerialOR.sic <- sic(HH=SerialOR.hh, HL=SerialOR.hl, LH=SerialOR.lh, LL=SerialOR.ll)
print(SerialOR.sic$Dvals)
plot(SerialOR.sic$SIC, do.p=FALSE, ylim=c(-1,1))

ParallelAND.hh <- pmax(T1.h, T2.h)
ParallelAND.hl <- pmax(T1.h, T2.l)
ParallelAND.lh <- pmax(T1.l, T2.h)
ParallelAND.ll <- pmax(T1.l, T2.l)
ParallelAND.sic <- sic(HH=ParallelAND.hh, HL=ParallelAND.hl, LH=ParallelAND.lh, 
  LL=ParallelAND.ll)
print(ParallelAND.sic$Dvals)
plot(ParallelAND.sic$SIC, do.p=FALSE, ylim=c(-1,1))

ParallelOR.hh <- pmin(T1.h, T2.h)
ParallelOR.hl <- pmin(T1.h, T2.l)
ParallelOR.lh <- pmin(T1.l, T2.h)
ParallelOR.ll <- pmin(T1.l, T2.l)
ParallelOR.sic <- sic(HH=ParallelOR.hh, HL=ParallelOR.hl, LH=ParallelOR.lh, 
  LL=ParallelOR.ll)
print(ParallelOR.sic$Dvals)
plot(ParallelOR.sic$SIC, do.p=FALSE, ylim=c(-1,1))
}
\keyword{ sft }
