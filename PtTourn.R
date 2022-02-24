
#########################################################
# Test d'ind�pendance des points tournants
#########################################################



PtTourn.test <- function (vec)
{
   n <- length(vec)
   nT <- 0
   for (i in 2:(n-1)){
       if (((vec[i] > vec[i-1]) & (vec[i] > vec[i+1])) | ( (vec[i] < vec[i-1])
        & (vec[i] < vec[i+1]))) 
          {nT <- nT+1}
           
}
Tobs <- ( nT - 2*(n-2)/3 ) / sqrt( (16*n-29)/90 )
p <- 2 * (1- pnorm(abs(Tobs)))
res <- c(n,nT,Tobs,p)
names(res) <- c("n","nT","stat"," p-valeur")
return(res)
}