\name{fPCAassessment}
\alias{fPCAassessment}
\title{Functional Principal Components Analysis for the Assessment Functions}
\description{Calculates the principle functions and scores for the workload assessment measure of performance by each individual in each condition.}
\usage{
fPCAassessment(sftData, dimensions, stopping.rule=c("OR", "AND", "STST"), 
               correct=c(TRUE,FALSE), fast=c(TRUE,FALSE), detection=TRUE, 
               register=c("median","mean","none"), plotPCs=FALSE, ...) 
}
\arguments{
  \item{sftData}{Data collected from a Double Factorial Paradigm experiment in standard form.}
  \item{dimensions}{The number of principal functions with which to represent the data.}
  \item{stopping.rule}{Indicates whether to use OR, AND or Single Target Self Terminating (STST) processing baseline to calculate individual assessment functions.}
  \item{correct}{Indicates whether to assess performance on correct trials.}
  \item{fast}{Indicates whether to use cumulative distribution functions or survivor functions to assess performance.}
  \item{detection}{Indicates whether to use a detection task baseline or a discrimination task baseline.}
  \item{register}{Indicates value to use for registering the assessment data.}
  \item{plotPCs}{Indicates whether or not to generate plots of the principal functions.}
  \item{...}{Arguments to be passed to plot function.}
}
\details{
Functional principal components analysis (fPCA) is an extension of standard principal components analysis to infinite dimensional (function) spaces.  Just as in standard principal components analysis, fPCA is a method for finding a basis set of lower dimensionality than the original space to represent the data.  However, in place of basis vectors, fPCA has basis functions.  Each function in the original dataset can then be represented by a linear combination of those bases, so that given the bases, the each datum is represented by a vector of its coefficients (or scores) in that linear combination.  

The assessment coefficient is a function across time, so the differences among assessment coefficients from different participants and/or conditions may be quite informative.  fPCA gives a well motivated method for representing those differences in a concise way.  The factor scores can be used to examine differences among assessment coefficients, accounting for variation across the entire function.

This function implements the steps outlines in Burns, Houpt, Townsend and Endres (2013) applied to the assessment functions defined in Townsend and Altieri (2012) and Donkin, Little, and Houpt (2013).  First, the data are shifted by subtracting the median response time within each condition for each participant, but across both single target and multiple target trials, so that the assessment curves will be registered. Second, each assessment coefficient is calculated with the shifted response times.  Next, the mean assessment coefficient is subtracted from each assessment coefficient, then the representation of the resulting assessment coefficients are translated to a b-spline basis.  The fPCA procedure extracts the basis function from the bspline space that accounts for the largest variation across the assessment coefficients, then the next basis function which must be orthogonal to the first but explains the most amount of variation in the assessment coefficients given that constraint and so on until the indicated number of basis have been extracted.  Once the assessment functions are represented in the reduced space, a varimax rotation is applied.

The assessment functions can be registered to the mean or median response time across all levels of workload but within each participant and condition, or the analyses can be performed without registration. 

For details on fPCA for the assessment coefficient, see Burns, Houpt, Townsend and Endres (2013).  For details on fPCA in general using R, see Ramsay, Hooker and Graves (2009).
}
\value{
  \item{Scores}{Data frame containing the Loading values for each participant and condition.}
  \item{MeanAT}{Object of class approxfun representing the mean At function.}
  \item{PF}{List of objects of class approxfun representing the principal functions.}
  \item{medianRT}{Size of shift used to register each assessment curve (median RT).}
}
\references{
Burns, D.M., Houpt, J.W., Townsend, J.T. & Endres, M.J. (2013).  Functional principal components analysis of workload assessment functions.  \emph{Behavior Research Methods}

Donkin, C, Little, D.R. and Houpt (2013). Assessing the effects of caution on the capacity of information processing.  \emph{Manuscript submitted for publication.}

Ramsay, J., Hooker, J. & Graves, S. (2009).  Functional Data Analysis with R and MATLAB.  New York, NY:  Springer.

Townsend, J.T. and Altieri, N. (2012).  An accuracy-response time capacity assessment function that measures performance against standard parallel predictions. \emph{Psychological Review, 3}, 500-516.

}
\author{
Joe Houpt <joseph.houpt@wright.edu>
}
\seealso{
\code{\link{assessment}}
\code{\link{fda}}
}
\examples{
\dontrun{
data(dots)
fPCAassessment(dots, dimensions=2, stopping.rule="OR", register="median",
               correct=TRUE, fast=FALSE, detection=TRUE, plotPCs=TRUE)
}
}
\keyword{ sft }
