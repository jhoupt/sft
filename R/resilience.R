resilience <- function(RT, CR=NULL, ratio=TRUE, rho=0) {
    times <- sort(unique(c(RT, recursive=TRUE))) 

    Hab <- estimateNAH(RT[[1]], CR[[1]])
    Hay <- estimateNAH(RT[[2]], CR[[2]])
    Hxb <- estimateNAH(RT[[3]], CR[[3]])

    rtest <- resilience.test(RT, CR, rho=0)
    if(ratio) {
      R <- Hab$H(times)/ (Hay$H(times) + Hxb$H(times))
      R[is.nan(R)] <- NA
      R[is.infinite(R)] <- NA
      R <- approxfun(times, R)
      return(list(Rt=R, Rtest=rtest))
    }
    else {
      R <- Hab$H(times) - (Hay$H(times) + Hxb$H(times))
      R[is.nan(R)] <- NA
      R[is.infinite(R)] <- NA
      Var.R <- Hab$Var(times) + Hay$Var(times) + Hxb$Var(times)
      Var.R <- approxfun(times, Var.R)
      return(list(Rt=R, Var=Var.R, Rtest=rtest))
    }
}


resilience.test <- function(RT, CR=NULL, rho=0) {
    METHOD <- "Resilience test"
    DNAME <- deparse(substitute(RT))
    ALTERNATIVE <- "response times are different than those predicted by the UCIP-OR model"

    RTab <- RT[[1]]
    RTay <- RT[[2]] 
    RTxb <- RT[[3]]
    allRT <- c(RT, recursive=TRUE)

    if ( is.null(CR) ) {
      CRab <- rep(1, length(RTab))
      CRay <- rep(1, length(RTay))
      CRxb <- rep(1, length(RTxb))
      allCR <- rep(1, length(allRT))
    } else {
      DNAME <- paste(DNAME, "and", deparse(substitute(CR)))
      CRab <- as.numeric(CR[[1]])
      CRay <- as.numeric(CR[[2]])
      CRxb <- as.numeric(CR[[3]])
      allCR <- c(CR, recursive=TRUE)
    }

    index <- numeric()
    for ( i in 1:3 ) {
      index <- c(index, rep(i, length(RT[[i]])) )
    }

    RTmat <- cbind( allRT, allCR, index)
    RT.sort <- sort(RTmat[,1], index.return=TRUE)
    cr.s <- RTmat[RT.sort$ix,2]
    cond.s <- RTmat[RT.sort$ix,3]
    tvec <- RT.sort$x

    Yab <- rep(0, length(tvec))
    Yay <- rep(0, length(tvec))
    Yxb <- rep(0, length(tvec))
    Jab <- rep(0, length(tvec))
    Jay <- rep(0, length(tvec))
    Jxb <- rep(0, length(tvec))
    for (i in 1:length(tvec)) { 
      s <- tvec[i]

      Jab[i] <- sum( RTab[CRab==1]==s )
      Jay[i] <- sum( RTay[CRay==1]==s )
      Jxb[i] <- sum( RTxb[CRxb==1]==s )

      Yab[i] <- sum(RTab >= s)
      Yay[i] <- sum(RTay >= s)
      Yxb[i] <- sum(RTxb >= s)
    }
    St <- (1-(Jab + Jay + Jxb) / (Yab + Yay + Yxb) )
    St <- cumprod(St)
    Lt <- St ^ rho * ( Yab * (Yay + Yxb ) / (Yab + Yay + Yxb) )

    numer <- 0 
    numer <- numer + sum( (Lt / Yab)[cond.s==1 & cr.s==1 & Yab > 0 & Yay > 0 & Yxb > 0]  )
    numer <- numer - sum( (Lt / Yay)[cond.s==2 & cr.s==1 & Yab > 0 & Yay > 0 & Yxb > 0]  ) 
    numer <- numer - sum( (Lt / Yxb)[cond.s==3 & cr.s==1 & Yab > 0 & Yay > 0 & Yxb > 0]  ) 
    
    denom <- 0 
    denom <- denom + sum( (Lt / Yab)[cond.s==1 & cr.s==1 & Yab > 0 & Yay > 0 & Yxb > 0]^2 )
    denom <- denom + sum( (Lt / Yay)[cond.s==2 & cr.s==1 & Yab > 0 & Yay > 0 & Yxb > 0]^2 )
    denom <- denom + sum( (Lt / Yxb)[cond.s==3 & cr.s==1 & Yab > 0 & Yay > 0 & Yxb > 0]^2 )
    denom <- sqrt(denom)

    STATISTIC <- numer/denom
    # names(STATISTIC) = "Rz"

    pval <- 2*min(pnorm(numer/denom),1-pnorm(numer/denom))
    rval <- list(statistic=STATISTIC, p.value=pval, alternative=ALTERNATIVE,
              method=METHOD, data.name=DNAME)
    class(rval) <- "htest"
    return(rval)
}

conflict.contrast <- function(RTay.H, RTay.L, RTxb.H, RTxb.L, CRay.H, CRay.L, CRxb.H, CRxb.L, rho) {
    #times <- sort(unique(c(RT, recursive=TRUE))) 
    times <- sort(unique(c(RTay.H, RTay.L, RTxb.H, RTxb.L)))

    Hay.H <- estimateNAH(RTay.H, CRay.H)
    Hay.L <- estimateNAH(RTay.L, CRay.L)
    Hxb.H <- estimateNAH(RTxb.H, CRxb.H)
    Hxb.L <- estimateNAH(RTxb.L, CRxb.L)

    cctest <- conflict.contrast.test(RTay.H, RTay.L, RTxb.H, RTxb.L, CRay.H, CRay.L, CRxb.H, CRxb.L, rho) 

    CC <- (Hay.L$H(times) - Hay.H$H(times)) + (Hxb.L$H(times) - Hxb.H$H(times))
    CC <- approxfun(times, CC)
    Var.CC <- (Hay.H$H(times) + Hay.L$H(times)) + (Hxb.H$H(times) + Hxb.L$H(times))
    Var.CC <- approxfun(times, Var.CC)
    return(list(CCt=CC, Var=Var.CC, CCtest=cctest))
}

conflict.contrast.test <- function(RTay.H, RTay.L, RTxb.H, RTxb.L, CRay.H, CRay.L, CRxb.H, CRxb.L, rho) {
    METHOD <- "Conflict Contrast test"
    ALTERNATIVE <- "The contrast conflict function is nonzero."
    RT <- list(RTay.H, RTay.L, RTxb.H, RTxb.L)
    DNAME <- deparse(substitute(RT))
    
    if ( is.null(CRay.H) ) {
      CRay.H <- rep(1, length(RTay.H))
    } else {
      DNAME <- paste(DNAME, deparse(CRay.H))
    }
    if ( is.null(CRay.L) ) {
      CRay.L <- rep(1, length(RTay.L))
    } else {
      DNAME <- paste(DNAME, deparse(CRay.L))
    }
    if ( is.null(CRxb.H) ) {
      CRxb.H <- rep(1, length(RTxb.H))
    } else {
      DNAME <- paste(DNAME, deparse(CRxb.H))
    }
    if ( is.null(CRxb.L) ) {
      CRxb.L <- rep(1, length(RTxb.L))
    } else {
      DNAME <- paste(DNAME, deparse(CRxb.L))
    }

    allRT <- c(RTay.H, RTay.L, RTxb.H, RTxb.L)
    allCR <- c(CRay.H, CRay.L, CRxb.H, CRxb.L)

    index <- c(rep(1, length(RTay.H)), rep(2, length(RTay.L)), rep(3, length(RTxb.H)), rep(4, length(RTxb.L))) 
    
    RTmat <- cbind( allRT, allCR, index)
    RT.sort <- sort(RTmat[,1], index.return=TRUE)
    cr.s <- RTmat[RT.sort$ix,2]
    cond.s <- RTmat[RT.sort$ix,3]
    tvec <- RT.sort$x

    Yay.H <- rep(0, length(tvec))
    Yay.L <- rep(0, length(tvec))
    Yxb.H <- rep(0, length(tvec))
    Yxb.L <- rep(0, length(tvec))

    Jay.H <- rep(0, length(tvec))
    Jay.L <- rep(0, length(tvec))
    Jxb.H <- rep(0, length(tvec))
    Jxb.L <- rep(0, length(tvec))

    for (i in 1:length(tvec)) { 
      s <- tvec[i]

      Jay.H[i] <- sum( RTay.H[CRay.H==1]==s )
      Jay.L[i] <- sum( RTay.L[CRay.L==1]==s )

      Jxb.H[i] <- sum( RTxb.H[CRxb.H==1]==s )
      Jxb.L[i] <- sum( RTxb.L[CRxb.L==1]==s )

      Yay.H[i] <- sum(RTay.H >= s)
      Yay.L[i] <- sum(RTay.L >= s)

      Yxb.H[i] <- sum(RTxb.H >= s)
      Yxb.L[i] <- sum(RTxb.L >= s)
    }

    St <- (1-(Jay.H + Jay.L + Jxb.H + Jxb.L) / (Yay.H + Yay.L + Yxb.H + Yxb.L) )
    St[CRay.H==0 | CRay.L==0 | CRxb.H==0 | CRxb.L==0] <- 1
    St <- cumprod(St)
    Lt <- St ^ rho * ( (Yay.H + Yxb.H )*(Yay.L + Yxb.L ) / (Yay.H + Yay.H + Yxb.H + Yxb.L) )

    numer <- 0 
    numer <- numer + sum( (Lt / Yay.H)[cond.s==1 & cr.s==1 & Yay.H > 0 & Yay.L > 0 & Yxb.H > 0 & Yxb.L > 0]  ) 
    numer <- numer - sum( (Lt / Yay.L)[cond.s==2 & cr.s==1 & Yay.H > 0 & Yay.L > 0 & Yxb.H > 0 & Yxb.L > 0]  ) 
    numer <- numer + sum( (Lt / Yxb.H)[cond.s==3 & cr.s==1 & Yay.H > 0 & Yay.L > 0 & Yxb.H > 0 & Yxb.L > 0]  ) 
    numer <- numer - sum( (Lt / Yxb.L)[cond.s==4 & cr.s==1 & Yay.H > 0 & Yay.L > 0 & Yxb.H > 0 & Yxb.L > 0]  ) 
    
    denom <- 0 
    denom <- denom + sum( (Lt / Yay.H)[cond.s==1 & cr.s==1 & Yay.H > 0 & Yay.L > 0 & Yxb.H > 0 & Yxb.L > 0]^2 )
    denom <- denom + sum( (Lt / Yay.L)[cond.s==2 & cr.s==1 & Yay.H > 0 & Yay.L > 0 & Yxb.H > 0 & Yxb.L > 0]^2 )
    denom <- denom + sum( (Lt / Yxb.H)[cond.s==3 & cr.s==1 & Yay.H > 0 & Yay.L > 0 & Yxb.H > 0 & Yxb.L > 0]^2 )
    denom <- denom + sum( (Lt / Yxb.L)[cond.s==4 & cr.s==1 & Yay.H > 0 & Yay.L > 0 & Yxb.H > 0 & Yxb.L > 0]^2 )
    denom <- sqrt(denom)

    STATISTIC <- numer/denom
    #names(STATISTIC) = "CCz"

    pval <- 2*min(pnorm(numer/denom),1-pnorm(numer/denom))
    rval <- list(statistic=STATISTIC, p.value=pval, alternative=ALTERNATIVE,
              method=METHOD) #, data.name=DNAME)
    class(rval) <- "htest"
    return(rval)

}

fPCAresilience<- function(sftData, dimensions, ratio=TRUE, register=c("median","mean","none"), plotPCs=FALSE, acc.cutoff=.70, ...) {
  subjects <- sort(unique(sftData$Subject))
  subjects <- factor(subjects)
  nsubjects <- length(subjects)
  conditions <- sort(unique(sftData$Condition))
  conditions <- factor(conditions)
  nconditions <- length(conditions)
  subj.out <- character()
  cond.out <- character()

  channels <- grep("Channel", names(sftData), value=T)
  nchannels <- length(channels)
  if(nchannels != 2) {
    stop("Invalid number of channels for resilience analysis.")
  }


  tvec <- seq(quantile(sftData$RT,.001), quantile(sftData$RT,.999), 
              length.out=1000)

  midpoint <- floor(length(tvec)/2)
  resAllMat <- numeric()
  varAllMat <- numeric()
  subjVec <- c()
  condVec <- c()

  allRT <- numeric()
  registervals <- numeric()
  good <- logical()

  RTlist <- vector("list", nchannels+1)
  CRlist <- vector("list", nchannels+1)

  ltyvec <- numeric()
  colvec <- numeric()
  condLegend <- levels(conditions)

  # Calculate resilience for each participant in each condition
  for ( cn in 1:nconditions ) {
    if (is.factor(conditions)) {cond <- levels(conditions)[cn]} else {cond <- conditions[cn] }
    condsubjects <- factor(with(sftData, sort(unique(Subject[Condition==cond]))))
    ncondsubjects <- length(condsubjects)
  for ( sn in 1:ncondsubjects ) {
      if (is.factor(condsubjects)) {subj <- levels(condsubjects)[sn]} else {subj <- condsubjects[sn] }

      subjVec <- c(subjVec, subj)
      condVec <- c(condVec, cond)

      ds <- sftData$Subject==subj & sftData$Condition==cond

      # Target Response Times
      usechannel <- ds & apply(sftData[, channels]>0, 1, all)
      RTlist[[1]] <- sftData$RT[usechannel & (sftData$RT < quantile(sftData$RT[usechannel], .975)) ]
      CRlist[[1]] <- sftData$Correct[usechannel & (sftData$RT < quantile(sftData$RT[usechannel], .975))]

      # Single Target Response Times
      for ( ch in 1:nchannels ) {
        usechannel <- ds & sftData[,channels[ch]]>0 & 
                      apply(as.matrix(sftData[,channels[-ch]] < 0), 1, all)
        RTlist[[ch+1]] <- sftData$RT[usechannel & sftData$RT < quantile(sftData$RT[usechannel], .975)]
        CRlist[[ch+1]] <- sftData$Correct[usechannel & sftData$RT < quantile(sftData$RT[usechannel], .975)]
      }
      

      # Check to make sure accuracy on each condition is higher than acc.cutoff
      if(any(lapply(CRlist, mean)<acc.cutoff) | any(lapply(RTlist, length) < 10) ) {
        good <- c(good, FALSE)
        resAllMat <- rbind(resAllMat, rep(NA, length(tvec)))
        varAllMat <- rbind(varAllMat, rep(NA, length(tvec)))
        next
      } else{
        good <- c(good, TRUE)
      }

      # Tracks the amount of offset for each resilience function  (for registering)
      if (register == "median") {
        registervals <- c(registervals, mean(median(RTlist[[1]], median(c(RTlist[2:nconditions],recursive=TRUE)))) )
        shiftn <- midpoint - max( which(tvec < tail(registervals,1)))
      } else if (register == "median") {
        registervals <- c(registervals, mean(c(RTlist,recursive=TRUE)) )
        shiftn <- midpoint - max( which(tvec < tail(registervals,1)))
      } else {
        shiftn <- 0
      }


      resout <- resilience(RTlist, CRlist, ratio=ratio)

      subj.out <- c(subj.out, subj)
      cond.out <- c(cond.out, cond)
      ltyvec <- c(ltyvec, sn)
      colvec <- c(colvec, cn)

      if (ratio) {
        tmin <- max( c(lapply(RTlist, quantile, probs=c(.01)), recursive=TRUE), na.rm=TRUE)
        tmax <- min( c(lapply(RTlist, quantile, probs=c(.99)), recursive=TRUE), na.rm=TRUE)
        ct <- resout$Rt(tvec)
        #ct[tvec < tmin] <- mean(ct[tmin:(tmin+10)])
        #ct[tvec > tmax] <- mean(ct[(tmax+10):tmax])
        ct[tvec < tmin] <- NA
        ct[tvec > tmax] <- NA
        if (register != "none") {
          resAllMat <- rbind(resAllMat, shift(ct, shiftn))
        } else { 
          resAllMat <- rbind(resAllMat, ct)
        } 
      } else {
        if (register != "none") {
          varAllMat <- rbind(varAllMat, shift(resout$Var(tvec), shiftn))
          resAllMat <- rbind(resAllMat, shift(resout$Rt(tvec), shiftn))
        } else { 
          varAllMat <- rbind(varAllMat, resout$Var(tvec))
          resAllMat <- rbind(resAllMat, resout$Rt(tvec))
        } 
      }
    }
  }

  if(register != "none") {
    tvec <- tvec - midpoint
  }
  
  tmin <- min(tvec[!apply(is.na(resAllMat[good,]), 2, all)])
  tmax <- max(tvec[!apply(is.na(resAllMat[good,]), 2, all)])
  tgood <- tvec[tvec >= tmin & tvec <= tmax]
  resGoodMat <- resAllMat[good,tvec >= tmin & tvec <= tmax]
  if (!ratio) { 
    varGoodMat <- varAllMat[good,tvec >= tmin & tvec <= tmax]
  }
  k <- dim(resGoodMat)[1]


  if(plotPCs) {
    dev.new()
    par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0))
    matplot(tgood, t(resGoodMat), type='l', lty=ltyvec, col=colvec, 
            #xlim=xbound,#c(tmin, tmax+50),
            xlim=c(tmin,1300),# ylim=c(-4,4),
            main="Resilience", ylab="C(t)", xlab="Time (Adjusted)")
    if(nconditions <= 5) {
      legend("topright", legend=condLegend, lty=1, col=1:5, cex=.9)
    }
  }

  # Replace NA values in each function with the average resilience across functions.
  resGoodmn <- apply(resGoodMat, 2, mean, na.rm=TRUE)
  for (i in 1:k) {
    resGoodMat[i, is.na(resGoodMat[i,])]  <- resGoodmn[is.na(resGoodMat[i,])]
  }
  if(!ratio) {
    varGoodMat[i, is.na(resGoodMat[i,])]  <- 0
  }
  

  #  subtract mean (across participants and conditions) resilience function 
  resGoodmn <- apply(resGoodMat, 2, mean)
  resGoodMat <- resGoodMat - matrix(resGoodmn, k, length(tgood), byrow=T)

  
  if(plotPCs) {
    dev.new()
    par(mfrow=c(1,2), mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0))
    plot(c(tmin-1000, tmax+1000), c(0,0), type='l', 
           xlim=c(tmin, tmax), #ylim=c(-4,4),
           main="Mean C(t)", ylab="C(t)", xlab="Time (Adjusted)")
    lines(tgood, resGoodmn, lwd=2)
    if (ratio) { abline(0,0, lty=1, col=grey(.4)) }
    matplot(tgood, t(resGoodMat), type='l', lty=ltyvec, col=colvec, 
            #xlim=xbound,#c(tmin, tmax+50),
            xlim=c(tmin,tmax),# ylim=c(-4,4),
            main="C(t)-Mean C(t)", ylab="C(t)", xlab="Time (Adjusted)")
    if(nconditions <= 5) {
      legend("topright", legend=condLegend, lty=1, col=1:5, cex=.9)
    }
  }

  wtvec <- rep(1, length(tgood))
  #if (ratio) {
  #  wtvec <- rep(1, length(tgood))
  #} else {
  #  wtvec <- apply(varGoodMat, 2, sum, na.rm=TRUE) 
  #  wtvec[wtvec<1E-4] <- 1E-4
  #  wtvec <- 1/wtvec
  #  wtvec[is.na(wtvec)] <- 0
  #  wtvec <- wtvec / sum(wtvec)

  #  if (OR) { 
  #    xbound <- c(min(tgood), min(tgood[which(wtvec < 1E-6)]))
  #  } else {
  #    xbound <- c(max(tgood[which(wtvec < 1E-6)]), max(tgood))
  #  }

  #  if(plotPCs) {
  #      dev.new()
  #      par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0))
  #      plot(tgood, wtvec, col='forestgreen', type='l',
  #          xlim=c(tmin, 1300),
  #          main="Weighting Function", xlab="Time (Adjusted)", ylab="")
  #  }
  #}

  
  wtGoodMat <- t(resGoodMat)
  #wtGoodMat <- t(resGoodMat) * matrix(wtvec, nrow=length(tgood), ncol=sum(good))

  if( sum(good)-1 > 4) {
    basis <- create.bspline.basis(rangeval=c(min(tgood),max(tgood)), 
                                  nbasis=sum(good)-1, norder=4)
  } else {
    basis <- create.bspline.basis(rangeval=c(min(tgood),max(tgood)), 
                                  nbasis=sum(good)-1, norder=sum(good)-1)
  }
  resGoodfd <- smooth.basis(tgood, wtGoodMat, basis)
  pcastrGood <- pca.fd(resGoodfd$fd,dimensions)
  if ( dimensions > 1) { 
    pcastrGoodVarmx <- varmx.pca.fd(pcastrGood)
  } else {
    pcastrGoodVarmx <- pcastrGood
  }
 
  if(plotPCs) {
    values <- pcastrGood$values
    dev.new()
    par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0))
    plot(1:5, values[1:5]/sum(values),
         xlim=c(1, 5), ylim=c(0,1), pch=19,
         main="Scree Plot", xlab="Eigenfunction", ylab="Variance Accounted For")
    lines(1:5, values[1:5]/sum(values))
  }
  

  harmmat <- eval.fd(tgood, pcastrGood$harmonics)
  harmmat <- harmmat / (wtvec %*% matrix(1, 1, dimensions))
  facmult <- apply(abs(pcastrGood$scores), 2, mean)

  harmmatV <- eval.fd(tgood, pcastrGoodVarmx$harmonics)
  harmmatV <- harmmatV / (wtvec %*% matrix(1, 1, dimensions))
  facmultV <- apply(abs(pcastrGoodVarmx$scores), 2, mean)

  scoreout <- data.frame(subjVec,condVec)
  for ( i in 1:dimensions) {
    scoreout[[i+2]] <- rep(NA, length(scoreout[[1]]))
    scoreout[[i+2]][good] <- pcastrGood$scores[,i]
  }
  names(scoreout) <- c("Subject","Condition",paste("D",1:dimensions,sep=""))

  scoreoutV <- data.frame(subjVec,condVec)
  for ( i in 1:dimensions) {
    scoreoutV[[i+2]] <- rep(NA, length(scoreoutV[[1]]))
    scoreoutV[[i+2]][good] <- pcastrGoodVarmx$scores[,i]
  }
  names(scoreoutV) <- c("Subject","Condition",paste("D",1:dimensions,sep=""))

  pflist <- vector("list", length=dimensions)
  for (ifac in 1:dimensions) {
    pflist[[ifac]] <- approxfun(tgood,harmmatV[,ifac])
  }

  if(plotPCs) {
    if (ratio) { ylim<-c(0,mean(resGoodmn)+max(facmult)) } else { ylim=c(-1, mean(resGoodmn)+max(facmult)) }
    dev.new()
    par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0), mfrow=c(dimensions,3))
    for ( ifac in 1:dimensions) {
      mainstr <- paste("PC", ifac, "-", floor(100*pcastrGood$varprop[ifac]), "%")
      Wveci <- resGoodmn + facmult[ifac]* harmmat[,ifac]
      plot(tgood, Wveci, type='l', lty=2, main="", xlab="Time (Adjusted)", ylab="",
           ylim=ylim, xlim=c(tmin, tmax))
      lines(tgood, resGoodmn)
      abline(0,0, col=grey(.4))
      mtext(mainstr, side=2, line=1)

      if(ifac==1) {
        mtext("Component Function", side=3, line=.5)
        legend("topright", c("Component", "Mean"), lty=c(2,1))
      }

      plot(tgood, Wveci - resGoodmn, type='l', main="", xlab="Time (Adjusted)", ylab="",
           ylim=ylim, xlim=c(tmin, tmax))
      abline(0,0, col=grey(.4))
      if(ifac==1) {mtext("Component - Mean", side=3, line=.5)}

      plot(scoreout$Subject, scoreout[[ifac+2]], type="n", #ylim=c(-2,2),
           xaxt='n', ylab="", xlab="Subject")
      axis(1,at=1:10, labels=rep("",10), las=0, cex=.1, tck=-.02)
      mtext(side=1, 1:10, at=1:10, line=.05, cex=.7)
      text(scoreout$Subject, scoreout[[ifac+2]], labels=scoreout$Condition,
           col=colvec)
      if(ifac==1) {mtext("Score", side=3, line=.5)}
    }

    dev.new()
    par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0), mfrow=c(dimensions,3))
    for ( ifac in 1:dimensions) {
      mainstr <- paste("PC", ifac, "-", floor(100*pcastrGoodVarmx$varprop[ifac]), "%")

      Wveci <- resGoodmn + facmultV[ifac]* harmmatV[,ifac]

      plot(tgood, Wveci, type='l', lty=2, main="", xlab="Time (Adjusted)", ylab="",
           ylim=ylim, xlim=c(tmin, tmax))
      lines(tgood, resGoodmn)
      abline(0,0, col=grey(.4))
      mtext(mainstr, side=2, line=1)

      if(ifac==1) {
        mtext("Component Function", side=3, line=.5)
        legend("topright", c("Component", "Mean"), lty=c(2,1))
      }

      plot(tgood, Wveci - resGoodmn, type='l', main="", xlab="Time (Adjusted)", ylab="",
           ylim=ylim, xlim=c(tmin, tmax))
      abline(0,0, col=grey(.4))
      if(ifac==1) {mtext("Component - Mean", side=3, line=.5)}

      plot(scoreout$Subject, scoreoutV[[ifac+2]], type="n", #ylim=c(-2,2),
        xaxt='n', ylab="", xlab="Subject")
      axis(1,at=1:10, labels=rep("",10), las=0, cex=.1, tck=-.02)
      mtext(side=1, 1:10, at=1:10, line=.05, cex=.7)
      text(scoreout$Subject, scoreoutV[[ifac+2]], labels=scoreout$Condition,
           col=colvec)
      if(ifac==1) {mtext("Score", side=3, line=.5)}
    }
  }
  
  return(list(Scores=scoreoutV, MeanCT=approxfun(tgood,resGoodmn), PF=pflist, shiftRT=registervals))
}


shift <- function(x, n, wrap=FALSE) {
  # Shift an array (x) by n
  #  positive n shift right; negative n shfit left

  if (abs(n) > length(x) ) {
    if (!wrap ) { return( rep(NA, length(x))) }
    n <- n %% length(x)
  }

  if ( n >= 0 ) {
    s  <- length(x)-n +1
    if (wrap) {
      xout <- c( x[s:length(x)], x[1:(s-1)])
    } else {
      xout <- c(rep(NA,n), x[1:(s-1)])
    }
  } else {
    s <- abs(n)+1
    if (wrap) {
      xout <- c( x[s:length(x)], x[1:(s-1)])
    } else {
      xout <- c( x[s:length(x)], rep(NA, abs(n)))
    }
  }
  return(xout)
}


fPCArdiff<- function(sftData, dimensions, ratio=TRUE, register=c("median","mean","none"), plotPCs=FALSE, acc.cutoff=.70, ...) {
  subjects <- sort(unique(sftData$Subject))
  subjects <- factor(subjects)
  nsubjects <- length(subjects)
  conditions <- sort(unique(sftData$Condition))
  conditions <- factor(conditions)
  nconditions <- length(conditions)
  subj.out <- character()
  cond.out <- character()
  
  channels <- grep("Channel", names(sftData), value=T)
  nchannels <- length(channels)
  if(nchannels != 2) {
    stop("Invalid number of channels for resilience analysis.")
  }
  
  
  tvec <- seq(quantile(sftData$RT,.001), quantile(sftData$RT,.999), 
              length.out=1000)
  
  midpoint <- floor(length(tvec)/2)
  resAllMat <- numeric()
  varAllMat <- numeric()
  subjVec <- c()
  condVec <- c()
  
  allRT <- numeric()
  registervals <- numeric()
  good <- logical()
  
  RTlistL <- vector("list", nchannels+1)
  CRlistL <- vector("list", nchannels+1)
  RTlistH <- vector("list", nchannels+1)
  CRlistH <- vector("list", nchannels+1)
  
  
  ltyvec <- numeric()
  colvec <- numeric()
  condLegend <- levels(conditions)
  
  # Calculate resilience for each participant in each condition
  for ( cn in 1:nconditions ) {
    if (is.factor(conditions)) {cond <- levels(conditions)[cn]} else {cond <- conditions[cn] }
    condsubjects <- factor(with(sftData, sort(unique(Subject[Condition==cond]))))
    ncondsubjects <- length(condsubjects)
    for ( sn in 1:ncondsubjects ) {
      if (is.factor(condsubjects)) {subj <- levels(condsubjects)[sn]} else {subj <- condsubjects[sn] }
      
      subjVec <- c(subjVec, subj)
      condVec <- c(condVec, cond)
      
      ds <- sftData$Subject==subj & sftData$Condition==cond
      
      # Double Target Response Times
      usechannel <- ds & apply(sftData[, channels] > 0, 1, all)
      RTlistH[[1]] <- sftData$RT[usechannel & (sftData$RT < quantile(sftData$RT[usechannel], .975)) ]
      RTlistL[[1]] <- sftData$RT[usechannel & (sftData$RT < quantile(sftData$RT[usechannel], .975)) ]
      
      CRlistH[[1]] <- sftData$Correct[usechannel & (sftData$RT < quantile(sftData$RT[usechannel], .975))]
      CRlistL[[1]] <- sftData$Correct[usechannel & (sftData$RT < quantile(sftData$RT[usechannel], .975))]
      
      # Single Target Response Times
      for ( ch in 1:nchannels ) {
        usechannel <- ds & sftData[,channels[ch]] > 0 & 
          apply(as.matrix(sftData[,channels[-ch]] == -1), 1, all)
        RTlistL[[ch+1]] <- sftData$RT[usechannel & sftData$RT < quantile(sftData$RT[usechannel], .975)]
        CRlistL[[ch+1]] <- sftData$Correct[usechannel & sftData$RT < quantile(sftData$RT[usechannel], .975)]

        usechannel <- ds & sftData[,channels[ch]] > 0 & 
          apply(as.matrix(sftData[,channels[-ch]] == -2), 1, all)
        RTlistH[[ch+1]] <- sftData$RT[usechannel & sftData$RT < quantile(sftData$RT[usechannel], .975)]
        CRlistH[[ch+1]] <- sftData$Correct[usechannel & sftData$RT < quantile(sftData$RT[usechannel], .975)]
      }
      
      
      # Check to make sure accuracy on each condition is higher than acc.cutoff
      if(any(lapply(CRlistL, mean)<acc.cutoff) | any(lapply(RTlistL, length) < 10)  | 
         any(lapply(CRlistH, mean)<acc.cutoff) | any(lapply(RTlistH, length) < 10) ) {
        good <- c(good, FALSE)
        resAllMat <- rbind(resAllMat, rep(NA, length(tvec)))
        varAllMat <- rbind(varAllMat, rep(NA, length(tvec)))
        next
      } else{
        good <- c(good, TRUE)
      }
      
      
      # Tracks the amount of offset for each resilience function  (for registering)
      if (register == "median") {
        registervals <- c(registervals, mean(
            median(RTlistH[[1]], median(c(RTlistH[2:nconditions],recursive=TRUE))),
            median(RTlistL[[1]], median(c(RTlistL[2:nconditions],recursive=TRUE))) ) )
        shiftn <- midpoint - max( which(tvec < tail(registervals,1)))
      } else if (register == "mean") {
        registervals <- c(registervals, mean(c(RTlistH, RTlistL,recursive=TRUE)) )
        shiftn <- midpoint - max( which(tvec < tail(registervals,1)))
      } else {
        shiftn <- 0
      }
      
      resoutL <- resilience(RTlistL, CRlistL, ratio=ratio)
      resoutH <- resilience(RTlistH, CRlistH, ratio=ratio)
      ctL <- resoutL$Rt(tvec)
      ctH <- resoutH$Rt(tvec)
      ct <- ctH - ctL

      
      subj.out <- c(subj.out, subj)
      cond.out <- c(cond.out, cond)
      ltyvec <- c(ltyvec, sn)
      colvec <- c(colvec, cn)
      
      if (ratio) {
        tminL <- max( c(lapply(RTlistL, quantile, probs=c(.01)), recursive=TRUE), na.rm=TRUE)
        tminH <- max( c(lapply(RTlistH, quantile, probs=c(.01)), recursive=TRUE), na.rm=TRUE)
        tmin <- min(c(tminL, tminH))
        tmaxL <- min( c(lapply(RTlistL, quantile, probs=c(.99)), recursive=TRUE), na.rm=TRUE)
        tmaxH <- min( c(lapply(RTlistH, quantile, probs=c(.99)), recursive=TRUE), na.rm=TRUE)
        tmax <- max(c(tmaxL, tmaxH))

        #ct[tvec < tmin] <- mean(ct[tmin:(tmin+10)])
        #ct[tvec > tmax] <- mean(ct[(tmax+10):tmax])
        ct[tvec < tmin] <- NA
        ct[tvec > tmax] <- NA
        if (register != "none") {
          resAllMat <- rbind(resAllMat, shift(ct, shiftn))
        } else { 
          resAllMat <- rbind(resAllMat, ct)
        } 
      } else {
        vt <- resoutH$Var(tvec) + resoutL$Var(tvec)
        if (register != "none") {
          varAllMat <- rbind(varAllMat, shift(vt, shiftn))
          resAllMat <- rbind(resAllMat, shift(ct, shiftn))
        } else { 
          varAllMat <- rbind(varAllMat, vt)
          resAllMat <- rbind(resAllMat, ct)
        } 
      }
    }
  }
  
  if(register != "none") {
    tvec <- tvec - midpoint
  }
  
  tmin <- min(tvec[!apply(is.na(resAllMat[good,]), 2, all)])
  tmax <- max(tvec[!apply(is.na(resAllMat[good,]), 2, all)])
  tgood <- tvec[tvec >= tmin & tvec <= tmax]
  resGoodMat <- resAllMat[good,tvec >= tmin & tvec <= tmax]
  if (!ratio) { 
    varGoodMat <- varAllMat[good,tvec >= tmin & tvec <= tmax]
  }
  k <- dim(resGoodMat)[1]
  
  
  if(plotPCs) {
    dev.new()
    par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0))
    matplot(tgood, t(resGoodMat), type='l', lty=ltyvec, col=colvec, 
            #xlim=xbound,#c(tmin, tmax+50),
            xlim=c(tmin,tmax),# ylim=c(-4,4),
            main="Resilience Difference", ylab="C(t)", xlab="Time (Adjusted)")
    if(nconditions <= 5) {
      legend("topright", legend=condLegend, lty=1, col=1:5, cex=.9)
    }
  }
  
  # Replace NA values in each function with the average resilience across functions.
  resGoodmn <- apply(resGoodMat, 2, mean, na.rm=TRUE)
  for (i in 1:k) {
    resGoodMat[i, is.na(resGoodMat[i,])]  <- resGoodmn[is.na(resGoodMat[i,])]
  }
  if(!ratio) {
    varGoodMat[i, is.na(resGoodMat[i,])]  <- 0
  }
  
  
  #  subtract mean (across participants and conditions) resilience function 
  resGoodmn <- apply(resGoodMat, 2, mean)
  resGoodMat <- resGoodMat - matrix(resGoodmn, k, length(tgood), byrow=T)
  
  
  if(plotPCs) {
    dev.new()
    par(mfrow=c(1,2), mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0))
    plot(c(tmin-1000, tmax+1000), c(0,0), type='l', 
         xlim=c(tmin, tmax), #ylim=c(-4,4),
         main="Mean C(t)", ylab="C(t)", xlab="Time (Adjusted)")
    lines(tgood, resGoodmn, lwd=2)
    if (ratio) { abline(0,0, lty=1, col=grey(.4)) }
    matplot(tgood, t(resGoodMat), type='l', lty=ltyvec, col=colvec, 
            #xlim=xbound,#c(tmin, tmax+50),
            xlim=c(tmin,tmax),# ylim=c(-4,4),
            main="C(t)-Mean C(t)", ylab="C(t)", xlab="Time (Adjusted)")
    if(nconditions <= 5) {
      legend("topright", legend=condLegend, lty=1, col=1:5, cex=.9)
    }
  }
  
  wtvec <- rep(1, length(tgood))
  #if (ratio) {
  #  wtvec <- rep(1, length(tgood))
  #} else {
  #  wtvec <- apply(varGoodMat, 2, sum, na.rm=TRUE) 
  #  wtvec[wtvec<1E-4] <- 1E-4
  #  wtvec <- 1/wtvec
  #  wtvec[is.na(wtvec)] <- 0
  #  wtvec <- wtvec / sum(wtvec)
  
  #  if (OR) { 
  #    xbound <- c(min(tgood), min(tgood[which(wtvec < 1E-6)]))
  #  } else {
  #    xbound <- c(max(tgood[which(wtvec < 1E-6)]), max(tgood))
  #  }
  
  #  if(plotPCs) {
  #      dev.new()
  #      par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0))
  #      plot(tgood, wtvec, col='forestgreen', type='l',
  #          xlim=c(tmin, 1300),
  #          main="Weighting Function", xlab="Time (Adjusted)", ylab="")
  #  }
  #}
  
  
  wtGoodMat <- t(resGoodMat)
  #wtGoodMat <- t(resGoodMat) * matrix(wtvec, nrow=length(tgood), ncol=sum(good))
  
  if( sum(good)-1 > 4) {
    basis <- create.bspline.basis(rangeval=c(min(tgood),max(tgood)), 
                                  nbasis=sum(good)-1, norder=4)
  } else {
    basis <- create.bspline.basis(rangeval=c(min(tgood),max(tgood)), 
                                  nbasis=sum(good)-1, norder=sum(good)-1)
  }
  resGoodfd <- smooth.basis(tgood, wtGoodMat, basis)
  pcastrGood <- pca.fd(resGoodfd$fd,dimensions)
  if ( dimensions > 1) { 
    pcastrGoodVarmx <- varmx.pca.fd(pcastrGood)
  } else {
    pcastrGoodVarmx <- pcastrGood
  }
  
  if(plotPCs) {
    values <- pcastrGood$values
    dev.new()
    par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0))
    plot(1:5, values[1:5]/sum(values),
         xlim=c(1, 5), ylim=c(0,1), pch=19,
         main="Scree Plot", xlab="Eigenfunction", ylab="Variance Accounted For")
    lines(1:5, values[1:5]/sum(values))
  }
  
  
  harmmat <- eval.fd(tgood, pcastrGood$harmonics)
  harmmat <- harmmat / (wtvec %*% matrix(1, 1, dimensions))
  facmult <- apply(abs(pcastrGood$scores), 2, mean)
  
  harmmatV <- eval.fd(tgood, pcastrGoodVarmx$harmonics)
  harmmatV <- harmmatV / (wtvec %*% matrix(1, 1, dimensions))
  facmultV <- apply(abs(pcastrGoodVarmx$scores), 2, mean)
  
  scoreout <- data.frame(subjVec,condVec)
  for ( i in 1:dimensions) {
    scoreout[[i+2]] <- rep(NA, length(scoreout[[1]]))
    scoreout[[i+2]][good] <- pcastrGood$scores[,i]
  }
  names(scoreout) <- c("Subject","Condition",paste("D",1:dimensions,sep=""))
  
  scoreoutV <- data.frame(subjVec,condVec)
  for ( i in 1:dimensions) {
    scoreoutV[[i+2]] <- rep(NA, length(scoreoutV[[1]]))
    scoreoutV[[i+2]][good] <- pcastrGoodVarmx$scores[,i]
  }
  names(scoreoutV) <- c("Subject","Condition",paste("D",1:dimensions,sep=""))
  
  pflist <- vector("list", length=dimensions)
  for (ifac in 1:dimensions) {
    pflist[[ifac]] <- approxfun(tgood,harmmatV[,ifac])
  }
  
  if(plotPCs) {
    if (ratio) { ylim<-c(0,mean(resGoodmn)+max(facmult)) } else { ylim=c(-1, mean(resGoodmn)+max(facmult)) }
    dev.new()
    par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0), mfrow=c(dimensions,3))
    for ( ifac in 1:dimensions) {
      mainstr <- paste("PC", ifac, "-", floor(100*pcastrGood$varprop[ifac]), "%")
      Wveci <- resGoodmn + facmult[ifac]* harmmat[,ifac]
      plot(tgood, Wveci, type='l', lty=2, main="", xlab="Time (Adjusted)", ylab="",
           ylim=ylim, xlim=c(tmin, tmax))
      lines(tgood, resGoodmn)
      abline(0,0, col=grey(.4))
      mtext(mainstr, side=2, line=1)
      
      if(ifac==1) {
        mtext("Component Function", side=3, line=.5)
        legend("topright", c("Component", "Mean"), lty=c(2,1))
      }
      
      plot(tgood, Wveci - resGoodmn, type='l', main="", xlab="Time (Adjusted)", ylab="",
           ylim=ylim, xlim=c(tmin, tmax))
      abline(0,0, col=grey(.4))
      if(ifac==1) {mtext("Component - Mean", side=3, line=.5)}
      
      plot(scoreout$Subject, scoreout[[ifac+2]], type="n", #ylim=c(-2,2),
           xaxt='n', ylab="", xlab="Subject")
      axis(1,at=1:10, labels=rep("",10), las=0, cex=.1, tck=-.02)
      mtext(side=1, 1:10, at=1:10, line=.05, cex=.7)
      text(scoreout$Subject, scoreout[[ifac+2]], labels=scoreout$Condition,
           col=colvec)
      if(ifac==1) {mtext("Score", side=3, line=.5)}
    }
    
    dev.new()
    par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0), mfrow=c(dimensions,3))
    for ( ifac in 1:dimensions) {
      mainstr <- paste("PC", ifac, "-", floor(100*pcastrGoodVarmx$varprop[ifac]), "%")
      
      Wveci <- resGoodmn + facmultV[ifac]* harmmatV[,ifac]
      
      plot(tgood, Wveci, type='l', lty=2, main="", xlab="Time (Adjusted)", ylab="",
           ylim=ylim, xlim=c(tmin, tmax))
      lines(tgood, resGoodmn)
      abline(0,0, col=grey(.4))
      mtext(mainstr, side=2, line=1)
      
      if(ifac==1) {
        mtext("Component Function", side=3, line=.5)
        legend("topright", c("Component", "Mean"), lty=c(2,1))
      }
      
      plot(tgood, Wveci - resGoodmn, type='l', main="", xlab="Time (Adjusted)", ylab="",
           ylim=ylim, xlim=c(tmin, tmax))
      abline(0,0, col=grey(.4))
      if(ifac==1) {mtext("Component - Mean", side=3, line=.5)}
      
      plot(scoreout$Subject, scoreoutV[[ifac+2]], type="n", #ylim=c(-2,2),
           xaxt='n', ylab="", xlab="Subject")
      axis(1,at=1:10, labels=rep("",10), las=0, cex=.1, tck=-.02)
      mtext(side=1, 1:10, at=1:10, line=.05, cex=.7)
      text(scoreout$Subject, scoreoutV[[ifac+2]], labels=scoreout$Condition,
           col=colvec)
      if(ifac==1) {mtext("Score", side=3, line=.5)}
    }
  }
  
  return(list(Scores=scoreoutV, MeanCT=approxfun(tgood,resGoodmn), PF=pflist, shiftRT=registervals))
}
