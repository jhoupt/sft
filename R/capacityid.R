capacity.id <- function(dual.target.rt, no.target.rt, single.target.rts, dual.target.cr=NULL, no.target.cr=NULL, single.target.crs=NULL, ratio=TRUE) {
    n_single <- length(single.target.rts)

    if ( is.null(dual.target.cr) ) { 
      dual.target.cr <- rep(1, length(dual.target.rt))
    } 
    if ( is.null(no.target.cr) ) { 
      no.target.cr <- rep(1, length(no.target.rt))
    } 
    if ( is.null(single.target.crs) | (length(single.target.crs) != n_single) ) {
      single.target.crs <- vector("list", n_single)
      for( i in 1:n_single ) {
        single.target.crs[[i]] <- rep(1, length(single.target.rts[[i]]))
      }
    } 
    times <- sort(unique(c(dual.target.rt, no.target.rt, c(single.target.rts, recursive=TRUE)))) 

    rmtest <-  ucip.id.test(dual.target.rt, no.target.rt, single.target.rts, dual.target.cr, no.target.cr, single.target.crs)

    # Find Nelson-Aalen Reverse Cumulative Hazard Estimates
    NAK.dual <- estimateNAK(dual.target.rt, dual.target.cr) 
    NAK.no <- estimateNAK(no.target.rt, no.target.cr)
    NAK.single <- vector("list", n_single)
    for ( i in 1:n_single) { 
        NAK.single[[i]] <- estimateNAK(single.target.rts[[i]], single.target.crs[[i]])
    }

    # Calculate the and capacity coefficient
    denom <- NAK.dual$K(times) + NAK.no$K(times)
    if (!ratio) {
        Var.id <- NAK.dual$Var(times) + NAK.no$Var(times)
    }

    numer <- rep(0, length(times))
    for( i in 1:n_single ) { 
        numer <- numer + NAK.single[[i]]$K(times)
        if(!ratio) {
            Var.id <- Var.id + NAK.single[[i]]$Var(times)
	}
    }

    if (ratio) {
      C.id <- numer / denom
      C.id[is.nan(C.id)] <- NA
      C.id[is.infinite(C.id)] <- NA
      C.id <- approxfun(times, C.id)
      return( list(Ct=C.id, Ctest=rmtest) )
    } else {
      C.id <- denom - numer
      C.id <- approxfun(c(times,Inf), c(C.id,0))
      Var.id <- approxfun(c(times,Inf), c(Var.id,0))
      return( list(Ct=C.id, Var=Var.id, Ctest=rmtest) )
    }
}
