\name{estimate.bounds}
\alias{estimate.bounds}
\title{
Bounds on Response Time Cumulative Distribution Functions for Parallel Processing Models
}
\description{
Calculates the bounds on the range of cumulative distribution functions for response time data for parallel processing models under specified stopping rules (OR, AND, or Single-Target Self-Terminating).
}
\usage{
estimate.bounds(RT, CR = NULL, stopping.rule = c("OR","AND","STST"), 
  assume.ID=FALSE, numchannels=NULL, unified.space=FALSE)
}
\arguments{
  \item{RT}{A list of numeric response time arrays for the individual processing channels}
  \item{CR}{A list of correct/incorrect indicator arrays.  If NULL, assumes all are correct.}
  \item{stopping.rule}{A character string specifying the stopping rule for the parallel processing model; must be one of "OR", "AND", "STST". If NULL, then "OR" is the default model.}
  \item{assume.ID}{A logical indicating whether the individual channel distributions are assumed to be Identically Distributed (ID). If FALSE, non-ID distributions are estimated from the RT data. If TRUE, only the first array in RT is used, together with numchannels, to estimate the distributions.}
  \item{numchannels}{Number of channels in the parallel processing model when all channels are active. If NULL,number channels will be estimated equal to length of RT.}
  \item{unified.space}{A logical indicating whether the unified capacity space version of the bounds should be estimated.}
}
\details{
The \emph{estimate.bounds} function uses the response times from individual channels processing in isolation to estimae the response time distributions for an \eqn{n}{n}-channel parallel model. The input argument RT must be a list of numeric arrays, containing either one array for each of the \eqn{n}{n} channels to be estimated (so length(RT)=n), or it can have length(RT)=1 and the bounds can be found under an assumption that the \eqn{n}{n} channels are identically distributed to the data in RT. For this latter case, assume.ID=TRUE and numchannels=n (where \eqn{n \ge 2}{n \ge 2}) must be specified.

Standard unlimited capacity parallel processing models for \eqn{n}{n} simultaneously operating channels can produce a range of behavior, which is bounded by various functions derived from the probability distributions on each of the \eqn{i}{i}, for \eqn{i=1,\ldots,n}{i=1,...,n}, channels operating in isolation. These bounds depend on the stopping rule under which is the parallel model is assumed to be operating (OR, AND, or STST).



\emph{stopping.rule="OR"}

Let \eqn{F_{n}(t)=P[RT \le t]=P[min(T_{i}) \le t]}{F_n(t)=P[RT \le t]=P[min(T_i) \le t]} for \eqn{i=1,...,n}{i=1,...,n}, denote the cumulative distribution of response times under a minimum time (logical OR) stopping rule.
The general bounds for n-channel parallel processing under an OR stopping rule are:
\deqn{max_{i}F^{i}_{n-1}(t) \le F_{n}(t) \le min_{i,j}[F^{i}_{n-1}(t) + F^{j}_{n-1}(t) - F^{ij}_{n-2}(t)]}{max_{i}F^{i}_{n-1}(t) \le F_{n}(t) \le min_{i,j}[F^{i}_{n-1}(t) + F^{j}_{n-1}(t) - F^{ij}_{n-2}(t)].}

Under the assumption or conditions that the individual channels are identically distributed, this inequality chain simplifies to
\deqn{F_{n-1}(t) \le F_{n}(t) \le [2*F_{n-1}(t) - F_{n-2}(t)]}{F_{n-1}(t) \le F_{n}(t) \le [2*F_{n-1}(t) - F_{n-2}(t)].}

When the model under scrutiny has only \eqn{n=2}{n=2} channels, the inequality chain takes the form:
\deqn{F^{i}_{1}(t) \le F_{2}(t) \le [F^{i}_{1}(t) + F^{j}_{1}(t)]}{F^{i}_{1}(t) \le F_{2}(t) \le [F^{i}_{1}(t) + F^{j}_{1}(t)].}

\emph{stopping.rule="AND"}

Let \eqn{G_{n}(t)=P[RT \le t]=P[max(T_i) \le t]}{G_n(t)=P[RT \le t]=P[max(T_i) \le t]}, for \eqn{i=1,...,n}{i=1,...,n}, denote the cumulative distribution of repsonse times under a maximum time (logical AND, ehxaustive) stopping rule.
The general bounds for n-channel parallel processing under an AND stopping rule are:
\deqn{max_{i,j}[G^{i}_{n-1}(t) + G^{j}_{n-1}(t) - G^{ij}_{n-2}(t)] \le G_{n}(t) \le min_{i}G^{i}_{n-1}(t)}{max_{i,j}[G^{i}_{n-1}(t) + G^{j}_{n-1}(t) - G^{ij}_{n-2}(t)] \le G_{n}(t) \le min_{i}G^{i}_{n-1}(t).}

Under the assumption or conditions that the individual channels are identically distributed, this inequality chain simplifies to
\deqn{[2*G_{n-1}(t) - G_{n-2}(t)] \le G_{n}(t) \le G_{n-1}(t)}{[2*G_{n-1}(t) - G_{n-2}(t)] \le G_{n}(t) \le G_{n-1}(t).}

When the model under scrutiny has only \eqn{n=2}{n=2} channels, the inequality chain takes the form:
\deqn{[G^{i}_{1}(t) + G^{j}_{1}(t) - 1] \le G_{2}(t) \le G^{i}_{1}(t)}{[G^{i}_{1}(t) + G^{j}_{1}(t) - 1] \le G_{2}(t) \le G^{i}_{1}(t).}

\emph{stopping.rule="STST"}

Let \eqn{F_{k}(t)=P[RT \le t]}{F_k(t)=P[RT \le t]} denote the cumulative distribution of repsonse times under a single-target self-terminating (STST) stopping rule, where the target of interest is on processing channel \eqn{k}{k} among \eqn{n}{n} active channels.
The general bounds for n-channel parallel processing under an STST stopping rule are:
\deqn{\prod_{i=1}^{n} F_{1}(t) \le F_{k}(t) \le \sum_{i=1}^{n} F_{1}(t)}{\prod_{i=1}^{n} F_{1}(t) \le F_{k}(t) \le \sum_{i=1}^{n} F_{1}(t).}

Under the assumption or conditions that the individual channels are identically distributed, this inequality chain simplifies to
\deqn{[F_{1}(t)]^{n} \le F_{k}(t) \le n*F_{1}(t)}{[F_{1}(t)]^{n} \le F_{k}(t) \le n*F_{1}(t).}

When the model under scrutiny has only \eqn{n=2}{n=2} channels, the inequality chain takes the form:
\deqn{[F^{i}_{1}(t) * F^{j}_{1}(t)] \le F_{k}(t) \le [F^{i}_{1}(t) + F^{j}_{1}(t)]}{[F^{i}_{1}(t) * F^{j}_{1}(t)] \le F_{k}(t) \le [F^{i}_{1}(t) + F^{j}_{1}(t)].}
Note that in this case, \eqn{k=i}{k=i} or \eqn{k=j}{k=j}, but this may not be specifiable \emph{a priori} depending on experimental design.

Across all stopping rule conditions, violation of the upper bound indicates performance that is faster than can be predicted by an unlimited capacity parallel model. This may arise from positive (facilitatory) crosstalk between parallel channels, super capacity parallel processing, or some form of co-active architecture in the measured human response time data.

Violation of the lower bound indicates performance that is slower than predictd by an unlimited capacity parallel model. This may arise from negative (inhibitory) crosstalk between parallel channels, fixed capacity or limited capacity processing, or some form of serial architecture in the measured human response time data.

}
\value{
A list containing the following components:
  \item{Upper.Bound}{An object of class "approxfun"" representing the estimated upper bound on the cumulative distribution function for an unlimited capacity parallel model.}
  \item{Lower.Bound}{A object of class "approxfun" representing the estimated lower bound on the cumulative distribution function for an unlimited capacity parallel model.}

}
\references{
Blaha, L.M. & Townsend, J.T. (under review). On the capacity of single-target self-terminating processes.

Colonius, H. & Vorberg, D. (1994). Distribution inequalities for parallel models with unlimited capacity. \emph{Journal of Mathematical Psychology, 38}, 35-58. 

Grice, G.R., Canham, L., & Gwynne, J.W. (1984). Absense of a redundant-signals effect in a raction time task with divided attention. \emph{Perception & Psychophysics, 36}, 565-570.

Miller, J. (1982). Divided attention: Evidence for coactivation with redundant signals. \emph{Cognitive Psychology, 14}, 247-279.

Townsend, J.T. & Eidels, A. (2011). Workload capacity spaces: a unified methodology for response time measures of efficiency as workload is varied. \emph{Psychonomic Bulletin & Review, 18}, 659-681.

Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Townsend, J.T. & Wenger, M.J. (2004). A theory of interactive parallel processing: New capacity measures and predictions for a response time inequality series. \emph{Psychological Review, 111}, 1003--1035.
}

\author{
Leslie Blaha <leslie.blaha@us.af.mil>

Joe Houpt <joseph.houpt@wright.edu>
}
\seealso{
\code{\link{ucip.test}}
\code{\link{capacity.or}}
\code{\link{capacity.and}}
\code{\link{capacity.stst}}
\code{\link{approxfun}}
}
\examples{
#randomly generated data
rate1 <- .35
rate2 <- .3
rate3 <- .4
RT.paa <- rexp(100, rate1)
RT.apa <- rexp(100, rate2)
RT.aap <- rexp(100, rate3)
RT.or <- pmin(rexp(100, rate1), rexp(100, rate2), rexp(100, rate3))
RT.and <- pmax(rexp(100, rate1), rexp(100, rate2), rexp(100, rate3))
tvec <- sort(unique(c(RT.paa, RT.apa, RT.aap, RT.or, RT.and)))

or.bounds <- estimate.bounds(RT=list(RT.paa, RT.apa, RT.aap), CR=NULL, assume.ID=FALSE, 
  unified.space=FALSE)
and.bounds <- estimate.bounds(RT=list(RT.paa, RT.apa, RT.aap))

\dontrun{
#plot the or bounds together with a parallel OR model
matplot(tvec, 
  cbind(or.bounds$Upper.Bound(tvec), or.bounds$Lower.Bound(tvec), ecdf(RT.or)(tvec)),
  type='l', lty=1, ylim=c(0,1), col=2:4, main="Example OR Bounds", xlab="Time", 
  ylab="P(T<t)")
abline(1,0)
legend('topright', c("Upper Bound", "Lower Bound", "Parallel OR Model"), 
  lty=1, col=2:4, bty="n")

#using the dots data set in sft package
data(dots)
attach(dots)
RT.A <- dots[Subject=='S1' & Condition=='OR' & Channel1==2 & Channel2==0, 'RT']
RT.B <- dots[Subject=='S1' & Condition=='OR' & Channel1==0 & Channel2==2, 'RT']
RT.AB <- dots[Subject=='S1' & Condition=='OR' & Channel1==2 & Channel2==2, 'RT']
tvec <- sort(unique(c(RT.A, RT.B, RT.AB)))
Cor.A <- dots[Subject=='S1' & Condition=='OR' & Channel1==2 & Channel2==0, 'Correct']
Cor.B <- dots[Subject=='S1' & Condition=='OR' & Channel1==0 & Channel2==2, 'Correct']
Cor.AB <- dots[Subject=='S1' & Condition=='OR' & Channel1==2 & Channel2==2, 'Correct']
capacity <- capacity.or(list(RT.AB,RT.A,RT.B), list(Cor.AB,Cor.A,Cor.B), ratio=TRUE)
bounds <- estimate.bounds(list(RT.A,RT.B), list(Cor.A,Cor.B), unified.space=TRUE)

#plot unified capacity coefficient space
plot(tvec, capacity$Ct(tvec), type="l", lty=1, col="red", lwd=2)
lines(tvec, bounds$Upper.Bound(tvec), lty=2, col="blue", lwd=2)
lines(tvec, bounds$Lower.Bound(tvec), lty=4, col="blue", lwd=2)
abline(h=1, col="black", lty=1)
}
}
\keyword{ sft }
