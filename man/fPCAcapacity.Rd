\name{fPCAcapacity}
\alias{fPCAcapacity}
\title{Functional Principal Components Analysis for the Capacity Coefficient}
\description{Calculates the principle functions and scores for the workload capacity measure of performance by each individual in each condition.}
\usage{
fPCAcapacity(sftData, dimensions, acc.cutoff=.75, OR=NULL, 
  stopping.rule=c("OR","AND","STST"), ratio=TRUE, 
  register=c("median","mean","none"), plotPCs=FALSE, ...)
}
\arguments{
  \item{sftData}{Data collected from a Double Factorial Paradigm experiment in standard form.}
  \item{dimensions}{The number of principal functions with which to represent the data.}
  \item{acc.cutoff}{Minimum accuracy needed to be included in the analysis.}
  \item{OR}{Indicates whether to compare performance to an OR or AND processing baseline.  Provided for backwards compatibility for package version < 2. }
  \item{stopping.rule}{Indicates whether to use OR, AND or Single Target Self Terminating (STST) processing baseline to calculate individual capacity functions.}
  \item{ratio}{Whether to use ratio Ct or difference Ct.}
  \item{register}{Indicates value to use for registering the capacity data.}
  \item{plotPCs}{Indicates whether or not to generate plots of the principal functions.}
  \item{...}{Arguments to be passed to plot function.}
}
\details{
Functional principal components analysis (fPCA) is an extension of standard principal components analysis to infinite dimensional (function) spaces.  Just as in standard principal components analysis, fPCA is a method for finding a basis set of lower dimensionality than the original space to represent the data.  However, in place of basis vectors, fPCA has basis functions.  Each function in the original dataset can then be represented by a linear combination of those bases, so that given the bases, the each datum is represented by a vector of its coefficients (or scores) in that linear combination.  

The capacity coefficient is a function across time, so the differences among capacity coefficients from different participants and/or conditions may be quite informative.  fPCA gives a well motivated method for representing those differences in a concise way.  The factor scores can be used to examine differences among capacity coefficients, accounting for variation across the entire function.

This function implements the steps outlines in Burns, Houpt, Townsend and Endres (2013).  First, the data are shifted by subtracting the median response time within each condition for each participant, but across both single target and multiple target trials, so that the capacity curves will be registered. Second, each capacity coefficient is calculated with the shifted response times.  Next, the mean capacity coefficient is subtracted from each capacity coefficient, then the representation of the resulting capacity coefficients are translated to a b-spline basis.  The fPCA proceedure extracts the basis function from the bspline space that accounts for the largest variation across the capacity coefficients, then the next basis function which must be orthogonal to the first but explains the most amount of variation in the capacity coefficients given that constraint and so on until the indicated number of basis have been extracted.  Once the capacity functions are represented in the reduced space, a varimax rotation is applied.

The capacity functions can be registered to the mean or median response time across all levels of workload but within each participant and condition, or the analyses can be performed without registration. 

For details on fPCA for the capacity coefficient, see Burns, Houpt, Townsend and Endres (2013).  For details on fPCA in general using R, see Ramsay, Hooker and Graves (2009).
}
\value{
  \item{Scores}{Data frame containing the Loading values for each participant and condition.}
  \item{MeanCT}{Object of class approxfun representing the mean Ct function.}
  \item{PF}{List of objects of class approxfun representing the principal functions.}
  \item{medianRT}{Size of shift used to register each capacity curve (median RT).}
}
\references{
Burns, D.M., Houpt, J.W., Townsend, J.T. & Endres, M.J. (2013).  Functional principal components analysis of workload capacity functions.  \emph{Behavior Research Methods}

Ramsay, J., Hooker, J. & Graves, S. (2009).  Functional Data Analysis with R and MATLAB.  New York, NY:  Springer.
}
\author{
Joe Houpt <joseph.houpt@wright.edu>
Devin Burns <devburns@indiana.edu>
}
\seealso{
\code{\link{capacity.and}}
\code{\link{capacity.or}}
\code{\link{fda}}
}
\examples{
\dontrun{
data(dots)
fPCAcapacity(dots, dimensions=2,stopping.rule="OR", 
  plotPCs=TRUE)
}
}
\keyword{ sft }
