\name{dots}
\alias{dots}
\docType{data}
\title{RT and and Accuracy from a Simple Detection Task}
\description{
Data from a simple Double Factorial Paradigm task.
}
\usage{data(dots)}
\format{
  A data frame with 57600 observations on the following 6 variables.
  \describe{
    \item{\code{Subject}}{A character vector indicating the participant ID.}
    \item{\code{Condition}}{A character vector indicating whether participants could respond as soon as they detected either dot (OR) or both dots (AND).}
    \item{\code{Correct}}{A logical vector indicating whether or not the participant responded correctly.}
    \item{\code{RT}}{A numeric vector indicating the response time on a given trial.}
    \item{\code{Channel1}}{A numeric vector indicating the stimulus level for the upper dot.  0:  Absent; 1:  Low contrast (slow); 2:  High contrast (fast).}
    \item{\code{Channel2}}{A numeric vector indicating the stimulus level for the lower dot.  0:  Absent; 1:  Low contrast (slow); 2:  High contrast (fast).}
  }
}
\details{
These data include response time and accuracy from nine participants that completed two versions of a Double Factorial Paradigm task.  Stimuli were either two dots, one above fixation and one below, a single dot above fixation, a single dot below fixation, or a blank screen.  Each dot could be presented either high or low contrast when present.  In the OR task, participants were instructed to  respond 'yes' whenever they saw either dot and 'no' otherwise.  In the AND task, participants were instructed to respond 'yes' only when both dots were present and 'no' otherwise.  See Eidels et al. (2012) or Houpt & Townsend (2010) for a more thorough description of the task.
}
\source{
Eidels, A., Townsend, J. T., Hughes, H. C., & Perry, L. A. (2012). Complementary relationship between response times, response accuracy, and task requirements in a parallel processing system. \emph{Journal Cognitive Psychology}. Manuscript
(submitted for publication).
}
\references{
Houpt, J.W. & Townsend, J.T. (2010). The statistical properties of the survivor interaction contrast. \emph{Journal of Mathematical Psychology, 54}, 446-453.
}
\examples{
data(dots)
summary(dots)
\dontrun{
sicGroup(dots)
capacityGroup(dots)
}
}
\keyword{datasets}
