generateData <- function(nSubjects, nTrials, arch, stoprule, highlowratio=3) {
  nu.l <- 0.1
  nu.h <- nu.l * highlowratio
  alpha <- 30
  diff <- 1
  sop <- .5

  rData <- matrix(NA, 4 * nSubjects * nTrials, 6)
  rData[,1] <- rep(1:nSubjects, each=4*nTrials)
  rData[,2] <- 1
  rData[,3] <- rep(c(rep(2,2*nTrials), rep(1, 2*nTrials)), nSubjects)
  rData[,4] <- rep(c(rep(2,nTrials), rep(1, nTrials)), 2*nSubjects)

  allrt <- numeric()
  for( subj in 1:nSubjects) {
    if(arch=="coactive") {
      rt.hh <- rinvGauss(nTrials, nu=alpha/(nu.h+nu.h), lambda=.5*(alpha/diff)^2)
      rt.hl <- rinvGauss(nTrials, nu=alpha/(nu.h+nu.l), lambda=.5*(alpha/diff)^2)
      rt.lh <- rinvGauss(nTrials, nu=alpha/(nu.l+nu.h), lambda=.5*(alpha/diff)^2)
      rt.ll <- rinvGauss(nTrials, nu=alpha/(nu.l+nu.l), lambda=.5*(alpha/diff)^2)

    } else {
      # Channel1 High
      x1h1 <- rinvGauss(nTrials, nu=alpha/nu.h, lambda=(alpha/diff)^2)
      x1h2 <- rinvGauss(nTrials, nu=alpha/nu.h, lambda=(alpha/diff)^2)

      # Channel1 Low
      x1l1 <- rinvGauss(nTrials, nu=alpha/nu.l, lambda=(alpha/diff)^2)
      x1l2 <- rinvGauss(nTrials, nu=alpha/nu.l, lambda=(alpha/diff)^2)

      # Channel2 High
      x2h1 <- rinvGauss(nTrials, nu=alpha/nu.h, lambda=(alpha/diff)^2)
      x2h2 <- rinvGauss(nTrials, nu=alpha/nu.h, lambda=(alpha/diff)^2)

      # Channel2 Low
      x2l1 <- rinvGauss(nTrials, nu=alpha/nu.l, lambda=(alpha/diff)^2)
      x2l2 <- rinvGauss(nTrials, nu=alpha/nu.l, lambda=(alpha/diff)^2)

      if (arch == "parallel") {
        if (stoprule == "and") {
          rt.hh <- pmax(x1h1, x2h1)
          rt.hl <- pmax(x1h2, x2l1)
          rt.lh <- pmax(x1l1, x2h2)
          rt.ll <- pmax(x1l2, x2l2)
        }else if (stoprule == "or") {
          rt.hh <- pmin(x1h1, x2h1)
          rt.hl <- pmin(x1h2, x2l1)
          rt.lh <- pmin(x1l1, x2h2)
          rt.ll <- pmin(x1l2, x2l2)
        } else {
          cat("Unknown stopping rule!\n")
          return(NULL)
        }
      } else if (arch == "serial") {
        if (stoprule == "and") {
          rt.hh <- x1h1 + x2h1
          rt.hl <- x1h2 + x2l1
          rt.lh <- x1l1 + x2h2
          rt.ll <- x1l2 + x2l2
        }else if (stoprule == "or") {
          oneFirst <- runif(nTrials) < .5
          rt.hh <- oneFirst * x1h1 + (1-oneFirst)*x2h1
          rt.hl <- oneFirst * x1h2 + (1-oneFirst)*x2l1
          rt.lh <- oneFirst * x1l1 + (1-oneFirst)*x2h2
          rt.ll <- oneFirst * x1l2 + (1-oneFirst)*x2l2
        } else {
          cat("Unknown stopping rule!\n")
          return(NULL)
        }
      } else {
        cat("Unknown architecture!\n")
        return(NULL)
      }
    }
    allrt <- c(allrt, rt.hh, rt.hl, rt.lh, rt.ll)
  }
  rData[,5] <- allrt
  rData[,6] <- TRUE
  rData <- as.data.frame(rData)
  rData[,2] <- paste(arch,stoprule, sep="-")
  names(rData) <- c("Subject", "Condition", "Channel1", "Channel2", "RT", "Correct")
  rData$Subject <- as.factor(rData$Subject)
  rData$Condition <- as.factor(rData$Condition)
  return(rData)
}

