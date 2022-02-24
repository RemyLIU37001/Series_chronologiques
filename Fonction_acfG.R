

##########################################################
# Graphique d'afc avec seuil pour un test global
##########################################################


acfG <- function(x, lag.max=NULL) {
  n <- length(x)
  if (is.null(lag.max)) 
    lag.max <- floor(10*log10(n))
  lag.max <- as.integer(min(lag.max, n-1L))
  rho <- acf(x, lag.max=lag.max, plot=F)$acf
  if ( all.equal(rho[1],1)==TRUE )
    rho <- rho[-1]
  nh <- length(rho) 
  out <- sum( abs(rho) > 1.96/sqrt(n) ) 
  alpha.global <- 0.05
  alpha <- 1 - (1-alpha.global)^(1/nh)
  seuil.global <- qnorm(1-alpha/2, sd=1) / sqrt(n)
  acf(x, lag.max=lag.max)
  abline(h=seuil.global, lty=3, col="red")
  abline(h=-seuil.global, lty=3, col="red")
  #
  testG <- binom.test(out, n=nh, p=0.05)
  return(testG)
}