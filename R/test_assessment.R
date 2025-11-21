
rate.c <- .015
rate.i <- rate.c * .5#2/3

c1c.12 <- rexp(1000, rate.c)
c1i.12 <- rexp(1000, rate.i)
RT.pa <- pmin(c1c.12, c1i.12)
CR.pa <- c1c.12 < c1i.12

c1c.12.l <- rexp(1000, rate.c * .7)
c1i.12.l <- rexp(1000, rate.i * .9)
RT.px.limited <- pmin(c1c.12.l, c1i.12.l)
CR.px.limited <- c1c.12.l < c1i.12.l

c1c.12.u <- rexp(1000, rate.c)
c1i.12.u <- rexp(1000, rate.i)
RT.px.unlimited <- pmin(c1c.12.u, c1i.12.u)
CR.px.unlimited <- c1c.12.u < c1i.12.u

c1c.12.s <- rexp(1000, rate.c * 1.3)
c1i.12.s <- rexp(1000, rate.i * 1.1)
RT.px.super <- pmin(c1c.12.s, c1i.12.s)
CR.px.super <- c1c.12.s < c1i.12.s

tvec <- sort(unique(c(RT.pa, RT.px.limited, RT.px.unlimited, RT.px.super)))

a.correct.fast.limited <- assessment(RT = list(RT.px.limited, RT.pa),
                                     CR = list(CR.px.limited, CR.pa),
				     stopping.rule="STST",
				     correct = TRUE, fast = TRUE)

a.correct.fast.unlimited <- assessment(RT = list(RT.px.unlimited, RT.pa),
                                     CR = list(CR.px.unlimited, CR.pa),
				     stopping.rule="STST",
				     correct = TRUE, fast = TRUE)

a.correct.fast.super <- assessment(RT = list(RT.px.super, RT.pa),
                                     CR = list(CR.px.super, CR.pa),
				     stopping.rule="STST",
				     correct = TRUE, fast = TRUE)
par(mfrow=c(2,2))
plot(tvec, a.correct.fast.limited(tvec), type='l', col='red', lty=4, 
     xlim=c(0,quantile(tvec, .98)),
     ylim=c(0,2), ylab="A(t)", xlab="Time", main="Correct and Fast")
lines(tvec, a.correct.fast.unlimited(tvec), col='blue', lty=5)
lines(tvec, a.correct.fast.super(tvec), col='darkgreen', lty=1)
legend("topright", c("Super", "Unlimited", "Limited"), 
       lty=c(1,5, 4), col=c("darkgreen", "blue", "red"))
abline(h=1, col='grey')



a.correct.slow.limited <- assessment(RT = list(RT.px.limited, RT.pa),
                                     CR = list(CR.px.limited, CR.pa),
				     stopping.rule="STST",
				     correct = TRUE, fast = FALSE)

a.correct.slow.unlimited <- assessment(RT = list(RT.px.unlimited, RT.pa),
                                     CR = list(CR.px.unlimited, CR.pa),
				     stopping.rule="STST",
				     correct = TRUE, fast = FALSE)

a.correct.slow.super <- assessment(RT = list(RT.px.super, RT.pa),
                                     CR = list(CR.px.super, CR.pa),
				     stopping.rule="STST",
				     correct = TRUE, fast = FALSE)
plot(tvec, a.correct.slow.limited(tvec), type='l', col='red', lty=4, 
     xlim=c(0,quantile(tvec, .98)),
     ylim=c(0,2), ylab="A(t)", xlab="Time", main="Correct and Slow")
lines(tvec, a.correct.slow.unlimited(tvec), col='blue', lty=5)
lines(tvec, a.correct.slow.super(tvec), col='darkgreen', lty=1)
abline(h=1, col='grey')


a.incorrect.fast.limited <- assessment(RT = list(RT.px.limited, RT.pa),
                                     CR = list(CR.px.limited, CR.pa),
				     stopping.rule="STST",
				     correct = FALSE, fast = TRUE)

a.incorrect.fast.unlimited <- assessment(RT = list(RT.px.unlimited, RT.pa),
                                     CR = list(CR.px.unlimited, CR.pa),
				     stopping.rule="STST",
				     correct = FALSE, fast = TRUE)

a.incorrect.fast.super <- assessment(RT = list(RT.px.super, RT.pa),
                                     CR = list(CR.px.super, CR.pa),
				     stopping.rule="STST",
				     correct = FALSE, fast = TRUE)
plot(tvec, a.incorrect.fast.limited(tvec), type='l', col='red', lty=4, 
     xlim=c(0,quantile(tvec, .98)),
     ylim=c(0,2), ylab="A(t)", xlab="Time", main="Incorrect and Fast")
lines(tvec, a.incorrect.fast.unlimited(tvec), col='blue', lty=5)
lines(tvec, a.incorrect.fast.super(tvec), col='darkgreen', lty=1)
abline(h=1, col='grey')



a.incorrect.slow.limited <- assessment(RT = list(RT.px.limited, RT.pa),
                                     CR = list(CR.px.limited, CR.pa),
				     stopping.rule="STST",
				     correct = FALSE, fast = FALSE)

a.incorrect.slow.unlimited <- assessment(RT = list(RT.px.unlimited, RT.pa),
                                     CR = list(CR.px.unlimited, CR.pa),
				     stopping.rule="STST",
				     correct = FALSE, fast = FALSE)

a.incorrect.slow.super <- assessment(RT = list(RT.px.super, RT.pa),
                                     CR = list(CR.px.super, CR.pa),
				     stopping.rule="STST",
				     correct = FALSE, fast = FALSE)
plot(tvec, a.incorrect.slow.limited(tvec), type='l', col='red', lty=4, 
     xlim=c(0,quantile(tvec, .98)),
     ylim=c(0,2), ylab="A(t)", xlab="Time", main="Incorrect and Slow")
lines(tvec, a.incorrect.slow.unlimited(tvec), col='blue', lty=5)
lines(tvec, a.incorrect.slow.super(tvec), col='darkgreen', lty=1)
abline(h=1, col='grey')



