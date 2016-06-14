\name{sicGroup}
\alias{sicGroup}
\title{
SIC Analysis for a Group
}
\description{
Calculates the SIC for each individual in each condition of a DFP experiment.  The function will plot each individuals SIC and return the results of the test for stochastic dominance and the statistical test of SIC form.
}
\usage{
sicGroup(inData, sictest = "ks", domtest = "ks", plotSIC=TRUE, ...)
}
\arguments{
  \item{inData}{Data collected from a Double Factorial Paradigm experiment in standard form.}
  \item{sictest}{Which type of hypothesis test to use for SIC form. "ks" is the only test currently implemented.}
  \item{domtest}{Which type of hypothesis test to use for testing stochastic dominance relations, either as series of KS tests ("ks") or the dominance test based on Dirichlet process priors ("dp").  DP not yet implemented.}
  \item{plotSIC}{Indicates whether or not to generate plots of the survivor interaction contrasts.}
  \item{...}{Arguments to be passed to plot function.}
}
\details{
%  TO DO:  Add details of sftData format
See the help page for the \code{\link{sic}} function for details of the survivor interaction contrast.
}
\value{
  \item{overview}{Data frame summarizing the test outcomes for each participant and condition.}
  \item{Subject}{The participant identifier from inData.}
  \item{Condition}{The condition identifier from inData.}
  \item{Selective.Influence}{The results of the survivor function dominance test for selective influence.  Pass indicates HH < HL, LH and LL > LH, HL, but not HL, LH < HH and not LH, HL > LL,  where A < B indicates that A is significantly faster than B at the level of the distribution.  Ambiguous means neither HL, LH < HH, nor LH, HL > LL, but at least one of HH < HL, LH or LL >HL, LH did not hold.  Fail means that at least one of HL, LH < HH or HL, LH > LL.}
  \item{Positive.SIC}{Indicates whehter the SIC is significantly positive at any time.}
  \item{Negative.SIC}{Indicates whehter the SIC is significantly negative at any time.}
  \item{MIC}{Indicates whether or not the MIC is significantly non-zero.}
  \item{Model}{Indicates which model would predict the pattern of data, assuming selective influence.}
  \item{SICfn}{Matrix with each row giving the values of the of the estimated SIC for one participant in one condition for values of times.  The rows match the ordering of statistic.}
  \item{sic}{List with each element giving the result applying sic() to an individual in a condition.  sic has the same ordering as overview.}
  \item{times}{Times at which the SICs in SICfn are calculated.}
}
\references{
Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Houpt, J.W. & Townsend, J.T. (2010). The statistical properties of the survivor interaction contrast. \emph{Journal of Mathematical Psychology, 54}, 446-453.

Heathcote, A., Brown, S.D., Wagenmakers, E-J. & Eidels, A. (2010) Distribution-free tests of stochastic dominance for small samples. \emph{Journal of Mathematical Psychology, 54}, 454-463. 

Houpt, J.W., Blaha, L.M., McIntire, J.P., Havig, P.R. and Townsend, J.T. (2013). Systems Factorial Technology with R. \emph{Behavior Research Methods}.
}
\author{
Joe Houpt <joseph.houpt@wright.edu>
}
\seealso{
\code{\link{sic}}
\code{\link{capacityGroup}}
}
\examples{
\dontrun{
data(dots)
sicGroup(dots)
}
}
\keyword{ ~sft }
