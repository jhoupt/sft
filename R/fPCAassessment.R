fPCAassessment<- function(sftData, dimensions, 
                          stopping.rule=c("OR", "AND", "STST"), 
                          correct=c(TRUE,FALSE), fast=c(TRUE,FALSE), 
                          detection=TRUE,
                          register=c("median","mean","none"), 
                          plotPCs=FALSE, ...) {
  subjects <- sort(unique(sftData$Subject))
  subjects <- factor(subjects)
  n.subjects <- length(subjects)

  conditions <- sort(unique(sftData$Condition))
  conditions <- factor(conditions)
  n.conditions <- length(conditions)

  subj.out <- character()
  cond.out <- character()

  channels <- grep("Channel", names(sftData), value=TRUE)
  n.channels <- length(channels)
  if(n.channels < 2) {
    stop("Not enough channels for assessment analysis.")
  }

  register <- match.arg(register, c("mean","median","none"))

  rule <- match.arg(stopping.rule, c("OR","AND","STST"))
  if(rule == "OR") {
    capacity <- capacity.or
  } else if (rule == "AND") {
    capacity <- capacity.and 
  } else if (rule == "STST"){
    capacity <- capacity.stst
  } else  {
    stop("Please choose a valid stopping rule for fPCAassessment.")
  }

  if (rule!="STST") {
    # Currently only does present versus absent
    #  To be implemented:  separate tests for each factorial
    #  salience condition;  Negative numbers for distractor
    for ( ch in channels ) {
      if(is.factor(sftData[,ch])) {
        sftData[,ch] <- as.numeric(levels(sftData[,ch]))[sftData[,ch]]
      }
      sftData <- subset(sftData, sftData[,ch] >=0)
    }
  }

  if(rule=="STST") { n.channels <- 2 }

  times <- seq(quantile(sftData$RT,.001), quantile(sftData$RT,.999), 
              length.out=1000)

  midpoint <- 0# floor(length(times)/2)
  at.all.mat <- numeric()
  subj.vec <- c()
  cond.vec <- c()

  allRT <- numeric()
  registervals <- numeric()

  #good <- logical()
  if (rule=="STST") {
    RTlist <- vector("list", n.channels)
    CRlist <- vector("list", n.channels)
  } else {
    RTlist <- vector("list", n.channels+1)
    CRlist <- vector("list", n.channels+1)
  }

  ltyvec <- numeric()
  colvec <- numeric()
  condLegend <- levels(conditions)

  # Calculate capacity for each participant in each condition
  for ( cn in 1:n.conditions ) {
    cond <- levels(conditions)[cn]
    cond.subjects <- with(sftData, factor(Subject[Condition==cond]))
    n.cond.subjects <- length(levels(cond.subjects))
    if (n.cond.subjects == 0) { 
      next
    }
    for (sn in 1:n.cond.subjects) {
      subj <- levels(cond.subjects)[sn]
      subj.vec <- c(subj.vec, subj)
      cond.vec <- c(cond.vec, cond)

      

      if (rule == "STST") {
        # Target with distractors 
        use.channel <- sftData$Subject == subj & sftData$Condition == cond &
                       (apply(sftData[,channels]>0, 1, sum)==1) & 
                       (apply(sftData[,channels]<0, 1, sum)>0)
        RTlist[[1]] <- with(sftData, RT[use.channel & 
                            RT < quantile(RT[use.channel], .975)])
        CRlist[[1]] <- with(sftData, Correct[use.channel & 
                            RT < quantile(RT[use.channel], .975)])

        # Isolated Target 
        use.channel <- sftData$Subject == subj & sftData$Condition == cond &
                       (apply(sftData[,channels] >= 0, 1, all)) & 
                       (apply(sftData[,channels] != 0, 1, sum) == 1)
        RTlist[[2]] <- with(sftData, RT[use.channel & 
                            RT < quantile(RT[use.channel], .975)])
        CRlist[[2]] <- with(sftData, Correct[use.channel & 
                            RT < quantile(RT[use.channel], .975)])

      } else {
        # Target Response Times
        use.channel <- sftData$Subject == subj & sftData$Condition == cond &
                       apply(sftData[,channels] > 0, 1, all)
        RTlist[[1]] <- with(sftData, RT[use.channel & 
                            RT < quantile(RT[use.channel], .975)])
        CRlist[[1]] <- with(sftData, Correct[use.channel & 
                            RT < quantile(RT[use.channel], .975)])

        # Single Target Response Times
        for (ch in 1:n.channels) {
          use.channel <- sftData$Subject == subj & 
                         sftData$Condition == cond &
                         sftData[, channels[ch]] > 0 & 
                         apply(as.matrix(sftData[, channels[-ch]]==0), 1, 
                               all)
          RTlist[[ch+1]] <- with(sftData, RT[use.channel & 
                                 RT < quantile(RT[use.channel], .975)])
          CRlist[[ch+1]] <- with(sftData, Correct[use.channel & 
                                 RT < quantile(RT[use.channel], .975)])
        }
      }

# Track the amount of offset for each capacity function (for registration)
      if (register == "median") {
        registervals <- c(registervals, 
                          mean(median(RTlist[[1]]), 
                               median(c(RTlist[2:n.conditions], 
                                        recursive = TRUE))))
        shift.n <- midpoint - max( which(times < tail(registervals,1)))
      } else if (register == "mean") {
        registervals <- c(registervals, mean(c(RTlist,recursive=TRUE)) )
        shift.n <- midpoint - max( which(times < tail(registervals,1)))
      } else {
        shift.n <- 0
      }


      at.out <- assessment(RT=RTlist, CR=CRlist, stopping.rule=rule, 
                           correct=correct, fast=fast, detection=detection)

      subj.out <- c(subj.out, subj)
      cond.out <- c(cond.out, cond)
      ltyvec <- c(ltyvec, sn)
      colvec <- c(colvec, cn)

      t.min <- max(c(lapply(RTlist, quantile, probs=c(.01)), 
                     recursive=TRUE), na.rm=TRUE)
      t.max <- min(c(lapply(RTlist, quantile, probs=c(.99)), 
                     recursive=TRUE), na.rm=TRUE)
      at <- at.out(times)
      at[times < t.min] <- NA
      at[times > t.max] <- NA

      if (register == "none") { 
        at.all.mat <- rbind(at.all.mat, at)
      } else {
        at.all.mat <- rbind(at.all.mat, shift(at, shift.n))
      }
    }
  }

  if(register != "none") {
    times <- times - mean(registervals)
  }
  
  t.min <- min(times[!apply(is.na(at.all.mat), 2, all)])
  t.max <- max(times[!apply(is.na(at.all.mat), 2, all)])
  t.good <- times[times >= t.min & times <= t.max]
  at.good.mat <- at.all.mat[, times >= t.min & times <= t.max]

  k <- dim(at.good.mat)[1]

# Replace NA values in each function with the average assessment value 
# across functions.
  at.good.mn <- apply(at.good.mat, 2, mean, na.rm=TRUE)
  for (i in 1:k) {
    at.good.mat[i, is.na(at.good.mat[i,])] <- 
      at.good.mn[is.na(at.good.mat[i,])]
  }

# Subtract mean (across participants and conditions) assessment value
  at.good.mn <- apply(at.good.mat, 2, mean)
  at.good.mat <- at.good.mat - matrix(at.good.mn, k, length(t.good), 
                                      byrow=TRUE)
  
  basis <- create.bspline.basis(rangeval=c(t.min,t.max), 
                                nbasis= k - 1, norder=4)
  at.good.fd <- smooth.basis(t.good, t(at.good.mat), basis)
  pca.str.good <- pca.fd(at.good.fd$fd, dimensions)
  if (dimensions > 1) { 
    pca.str.good.varmx <- varmx.pca.fd(pca.str.good)
  } else {
    pca.str.good.varmx <- pca.str.good
  }
 
  if(plotPCs) {
    values <- pca.str.good$values
    dev.new()
    par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0))
    plot(1:5, values[1:5]/sum(values),
         xlim=c(1, 5), ylim=c(0,1), pch=19,
         main="Scree Plot", xlab="Eigenfunction", 
         ylab="Variance Accounted For")
    lines(1:5, values[1:5]/sum(values))
  }
  
  #harm.mat <- eval.fd(t.good, pca.str.good$harmonics)
  #fac.mult <- apply(abs(pca.str.good$scores), 2, mean)

  harm.mat.v <- eval.fd(t.good, pca.str.good.varmx$harmonics)
  fac.mult.v <- apply(abs(pca.str.good.varmx$scores), 2, mean)

  #scoreout <- data.frame(subj.vec, cond.vec)
  #for (i in 1:dimensions) {
  #  scoreout[[i + 2]] <- pca.str.good$scores[,i]
  #}
  #names(scoreout) <- c("Subject", "Condition", 
  #                     paste("D", 1:dimensions, sep=""))

  scoreout.v <- data.frame(subj.vec, cond.vec)
  for (i in 1:dimensions) {
    scoreout.v[[i + 2]] <- pca.str.good.varmx$scores[, i]
  }
  names(scoreout.v) <- c("Subject", "Condition",
                         paste("D", 1:dimensions, sep=""))

  pf.list <- vector("list", length=dimensions)
  for (ifac in 1:dimensions) {
    pf.list[[ifac]] <- approxfun(t.good,harm.mat.v[,ifac])
  }

  if(plotPCs) {
    y.upper1 <- max(at.good.mn+ 
                   matrix(rep(fac.mult.v, length(t.good)), 
                          length(t.good), dimensions, byrow=T) * harm.mat.v)
    y.upper1 <- 1.1 * max(y.upper1, max(at.good.mn))
    y.lower1 <- min(at.good.mn+ 
                   matrix(rep(fac.mult.v, length(t.good)), 
                          length(t.good), dimensions, byrow=T) * harm.mat.v)
    y.lower1 <- 0.9 * min(y.lower1, min(at.good.mn))
    ylim1=c(y.lower1, y.upper1)

    y.upper2 <- max(matrix(rep(fac.mult.v, length(t.good)), 
                          length(t.good), dimensions, byrow=T) * harm.mat.v)
    y.lower2 <- min(matrix(rep(fac.mult.v, length(t.good)), 
                          length(t.good), dimensions, byrow=T) * harm.mat.v)
    ylim2=c(.9*y.lower2, 1.1*y.upper2)
    dev.new()
    par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25, 0), 
        mfrow=c(dimensions, 3))
    for (ifac in 1:dimensions) {
      mainstr <- paste("PC", ifac, "-", 
                       floor(100*pca.str.good.varmx$varprop[ifac]), "%")
      fn.i <- at.good.mn + fac.mult.v[ifac]* harm.mat.v[,ifac]
      plot(t.good, fn.i, type='l', lty=2, main="", xlab="Time (Adjusted)", 
           ylab="", ylim=ylim, xlim=c(t.min, t.max))
      lines(t.good, at.good.mn)
      abline(0,0, col=grey(.4))
      mtext(mainstr, side=2, line=1)

      if(ifac == 1) {
        mtext("Component Function", side=3, line=.5)
        legend("topright", c("Component", "Mean"), lty=c(2,1))
      }

      plot(t.good, fn.i - at.good.mn, type='l', main="", 
           xlab="Time (Adjusted)", ylab="", ylim=ylim, xlim=c(t.min, t.max))
      abline(0,0, col=grey(.4))

      if(ifac ==1) {
        mtext("Component - Mean", side=3, line=.5)
      }

      plot(scoreout.v$Subject, scoreout.v[[ifac+2]], type="n", 
           xaxt='n', ylab="", xlab="Subject")
      axis(1, at=1:10, labels=rep("",10), las=0, cex=.1, tck=-.02)
      mtext(side=1, 1:10, at=1:10, line=.05, cex=.7)
      text(scoreout.v$Subject, scoreout.v[[ifac + 2]], 
           labels=scoreout.v$Condition, col=colvec)

      if(ifac==1) {
        mtext("Score", side=3, line=.5)
      }
    }
  }
  
  return(list(Scores=scoreout.v, MeanAT=approxfun(t.good, at.good.mn), 
              PF=pf.list, medianRT=registervals))
}
