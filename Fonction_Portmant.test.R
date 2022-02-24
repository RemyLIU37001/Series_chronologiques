

##########################################################
# Test d'indépendance de Portmanteau
##########################################################



Portmant.test <- function(vec)
{
  n <- length(vec)
  autocor <- acf(vec, plot = FALSE, na.action = na.pass)
  H <- floor(n/4)
  H <- min( length(autocor$acf) -1 , H)   
  stat <- n * sum(autocor$acf[2:(H+1)] ^2 )
  pval <- 1 - pchisq(stat, df=H)
  res <- c(H,stat,pval)
  names(res) <- c("d.f.","stat"," p-value")
  return(res)
}


PortmantLB.test <- function(vec)
{
  n <- length(vec)
  autocor <- acf(vec, plot = FALSE, na.action = na.pass)
  H <- floor(n/4)
  H <- min( length(autocor$acf) -1 , H) 
  h <- 1:H  
  stat <- n * (n+2) * sum(autocor$acf[2:(H+1)] ^2 / (n-h) )
  pval <- 1 - pchisq(stat, df=H)
  res <- c(H,stat,pval)
  names(res) <- c("d.f.","stat"," p-value")
  return(res)
}
