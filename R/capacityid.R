capacity.id <- function(dt.rt, nt.rt, st.rts, dt.cr=NULL, nt.cr=NULL, st.crs=NULL, ratio=TRUE) {
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
    times <- sort(unique(c(dt.rt, nt.rt, c(st.rts, recursive=TRUE)))) 

    rmtest <-  ucip.id.test(dt.rt, nt.rt, st.rts, dt.cr, nt.cr, st.crs)

    # Find Nelson-Aalen Reverse Cumulative Hazard Estimates
    NAK.dual <- estimateNAK(dt.rt, dt.cr) 
    NAK.no <- estimateNAK(nt.rt, nt.cr)
    NAK.single <- vector("list", n_single)
    for ( i in 1:n_single) { 
        NAK.single[[i]] <- estimateNAK(st.rts[[i]], st.crs[[i]])
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
