ucip.test <- function(RT, CR=NULL, OR=NULL, stopping.rule=c("OR","AND","STST"))  {
  METHOD <- "Houpt-Townsend UCIP test"
  DNAME <- deparse(substitute(RT))

  if (is.null(OR)) {
    rule <- match.arg(stopping.rule, c("OR","AND","STST"))
  } else if (OR ==TRUE) {
    rule <- "OR"
  } else {
    rule <- "AND"
  }
  ncond <- length(RT) 
  allRT <- c(RT, recursive=TRUE)
  if ( is.null(CR) ) {
    allCR <- rep(1, length(allRT))
  } else {
    DNAME <- paste(DNAME, "and", deparse(substitute(CR)))
    allCR <- c(CR, recursive=TRUE)
  }
  Nt <- length(allRT)

  index <- numeric()
  for ( i in 1:ncond ) {
    index <- c(index, rep(i, length(RT[[i]])) )
  }

  RTmat <- cbind( allRT, allCR, index)

  RT.sort <- sort(RTmat[,1], index.return=TRUE)
  cr.s <- RTmat[RT.sort$ix,2]
  cond.s <- RTmat[RT.sort$ix,3]
  tvec <- RT.sort$x
  
  if (rule=="OR") {
    ALTERNATIVE <- "response times are different than those predicted by the UCIP-OR model"

    # Y is the number of response times that have not yet occurred
    Yarr <- rep(0, Nt)
    Ymat <- matrix(0, ncond, Nt)
    
    for (j in 1:ncond) { 
      # Set the values of Y to represent the number of RTs that have not 
      #   occurred by time t for condition i
      for (i in 1:Nt ) {Ymat[j,i] <- sum(RTmat[index==j,1] >= tvec[i]) }
    }

    if (ncond > 2) {
      Wv <- Ymat[1,]*( apply(Ymat[2:ncond,], 2, sum)) / apply(Ymat, 2, sum)
    } else {
      Wv <- Ymat[1,]*( Ymat[2:ncond,]) / apply(Ymat, 2, sum)
    }

    # Calculate the numerator ( the difference of the weighted true AV performance
    #   and the weighted predicted UCIP performance
    numer <- sum(Wv[cond.s==1&cr.s==1]/Ymat[1,cond.s==1&cr.s==1])
    for (i in 2:ncond) {
      numer <- numer - sum(Wv[cond.s==i&cr.s==1]/Ymat[i,cond.s==i&cr.s==1])
    }

    # Calculate the denominator (the sum of the estimated variance for each 
    #  cumulative hazard function estimator
    denom <- 0
    for (i in 1:ncond) {
      denom <- denom + sum((Wv[cond.s==i&cr.s==1]/Ymat[i,cond.s==i&cr.s==1])^2)
    }
    denom <- sqrt(denom)
          
  } else {
    if (rule=="AND") {
    ALTERNATIVE <- "response times are different than those predicted by the UCIP-AND model"
    } else if (rule=="STST") {
    ALTERNATIVE <- "response times are different than those predicted by the UCIP-STST model"
    }
    Garr <- rep(0, Nt)
    Gmat <- matrix(0, ncond, Nt)
    for (j in 1:ncond) { 
      for (i in 1:Nt ) {Gmat[j,i] <- sum(RTmat[RTmat[,3]==j,1] <= tvec[i]) }
    }
    if (ncond > 2) {
      Wv <- Gmat[1,]*( apply(Gmat[2:ncond,], 2, sum)) / apply(Gmat, 2, sum)
    } else {
      Wv <- Gmat[1,]*( Gmat[2:ncond,]) / apply(Gmat, 2, sum)
    }
    numer <- -1 * sum(Wv[cond.s==1&cr.s==1]/Gmat[1,cond.s==1&cr.s==1])
    for (i in 2:ncond) {
      numer <- numer + sum(Wv[cond.s==i&cr.s==1]/Gmat[i,cond.s==i&cr.s==1])
    }
    denom <- 0
    for (i in 1:ncond) {
      denom <- denom + sum((Wv[cond.s==i&cr.s==1]/Gmat[i,cond.s==i&cr.s==1])^2)
    }
    denom <- sqrt(denom)
  }
  STATISTIC <- numer/denom
  names(STATISTIC) = "z"

  pval <- 2*min(pnorm(numer/denom),1-pnorm(numer/denom))
  rval <- list(statistic=STATISTIC, p.value=pval, alternative=ALTERNATIVE,
            method=METHOD, data.name=DNAME)
  class(rval) <- "htest"
  return(rval)
}

ucip.id.test <- function(dt.rt, nt.rt, st.rts, dt.cr=NULL, nt.cr=NULL, st.crs=NULL) {
  METHOD <- "Houpt-Townsend UCIP test"

  n_single <- length(st.rts)

  if ( is.null(dt.cr) ) { 
    dt.cr <- rep(1, length(dt.rt))
  } 
  if ( is.null(nt.cr) ) { 
    nt.cr <- rep(1, length(nt.rt))
  } 
  if ( is.null(st.crs) | (length(st.crs) != n_single) ) {
    st.crs <- vector("list", n_single)
    for( i in 1:n_single ) {
      st.crs[[i]] <- rep(1, length(st.rts[[i]]))
    }
  } 

  allRT <- c(dt.rt, nt.rt, c(st.rts, recursive=TRUE))
  allCR <- c(dt.cr, nt.cr, c(st.crs, recursive=TRUE))
  #allRT <- c(RT, recursive=TRUE)
  Nt <- length(allRT)

  index <- c(rep(1, length(dt.rt)), rep(2, length(nt.rt)))
  for ( i in 1:n_single) {
    index <- c(index, rep(i+2, length(st.rts[[i]])))
  }

  RTmat <- cbind( allRT, allCR, index)

  RT.sort <- sort(RTmat[,1], index.return=TRUE)
  cr.s <- RTmat[RT.sort$ix,2]
  cond.s <- RTmat[RT.sort$ix,3]
  tvec <- RT.sort$x
  
  ALTERNATIVE <- "response times are different than those predicted by the adjusted UCIP-AND model"
  Garr <- rep(0, Nt)
  Gmat <- matrix(0, 2+n_single, Nt)
  for (j in 1:(2+n_single)) { 
    for (i in 1:Nt ) {Gmat[j,i] <- sum(RTmat[RTmat[,3]==j,1] <= tvec[i]) }
  }

  Wv <- (Gmat[1,] + Gmat[2,])*(apply(Gmat[3:(2+n_single),], 2, sum)) / apply(Gmat, 2, sum)

  numer <- -1 * (sum(Wv[cond.s==1&cr.s==1]/Gmat[1,cond.s==1&cr.s==1]) + 
		 sum(Wv[cond.s==2&cr.s==1]/Gmat[2,cond.s==2&cr.s==1]))
  for (i in 1:n_single) {
    numer <- numer + sum(Wv[cond.s==(i+2)&cr.s==1]/Gmat[(i+2),cond.s==(i+2)&cr.s==1])
  }
  denom <- 0
  for (i in 1:(2+n_single)) {
    denom <- denom + sum((Wv[cond.s==i&cr.s==1]/Gmat[i,cond.s==i&cr.s==1])^2)
  }
  denom <- sqrt(denom)
  
  STATISTIC <- numer/denom
  names(STATISTIC) = "z"

  pval <- 2*min(pnorm(numer/denom),1-pnorm(numer/denom))
  rval <- list(statistic=STATISTIC, p.value=pval, alternative=ALTERNATIVE,
            method=METHOD)
  class(rval) <- "htest"
  return(rval)
}
