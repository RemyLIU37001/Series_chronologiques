rm(list = ls())
library(ggplot2)
library(forecast)
library(tseries)
library(TSA)

#Soit Xt le processus stationnaire solution de l’équation:

#                               Xt = 0.4Xt−1 + εt − θεt−1,
#  ε est un bruit blanc Gaussien de variance 1 et θ = ( moyenne de vos mois de naissance)/13 (vous
#  pouvez arrondir la valeur de θ à 2 chiffres après la virgule).
#  Donc pour nous Theta = (12+4)/13= 1.23

##1. Simuler sous R une réalisation de longeur n = 200 du processus stationnaire solution de cette
##   équation. Nommez-le simul.
simul<-arima.sim(list(order=c(1,0,1), ma=-1.23, ar=0.4), n=200)
plot(simul)

##2. Charger les données depuis le fichier ’data.txt’. Le fichier contient 200 valeurs de la variable x.
##   Créez une nouvelle série chronologique comme la somme de x + simul.

###2.1 Chargement les données 
x<- read.csv("data.txt",header = TRUE,sep = " ")
### On utilise le commande "read.csv" pour charger le fichier "data.txt", et le séparateur entre les données est 
### espace.
x<-ts(x)
plot(x)
###2.2 Création la nouvelle série chronologique
serie<-x+ simul

##3. Ajustez le modèle SARIMA approprié aux données en suivant les étapes décrites dans les diapositives
##   221 et 222 du cours.

###3.1 Tracer la série pour proposer une structure.
plot(serie)
#### Nous avons un modèle additif. 

###3.2 Tracer l'acf pour détecter une tendance ou une composante saisonnière.
acf(serie)
#### Dans la série que nous avons eu, nous pouvons trouver qu'il y a une tendance par l'acf, et la période c'est 10.

###3.3 On élimine d'abord la composante saisonnière.

####3.3.1 On calcule la période à l'aide des périodogrammes
periodog <- function(x, graph = T){
  freqs <- (2 * pi * (0:floor(length(x)/2)))/length(x)
  periodo <- (Mod(fft(x - mean(x)))^2/(2 * pi * length(x)))[1:(floor(length(x)/2) +
                                                                 1)]
  if (graph == T) {
    plot(freqs, periodo, type = "l", main = "Periodogram")
  }
  else {
    return(data.frame(periodo = periodo, freqs = freqs))
  }}
periodog(serie)
spectrum(serie)
#### Par ça nous pouvons savoir que la période de la saisonnalité c'est 10(par 
#### Periodogram c'est presque 0.1, mais nous pouvons savoir plus clair avec le spectrum(0.1) ) 

####3.3.2 On valide la période r à l'aide du test de saisonnalité
##### On excute les commandes 
Saison.test <- function(vec, d) {
  n <- length(vec)
  ns <- n%/%d
  n <- ns * d
  vec <- vec[1:n]
  rangs.vec <- matrix(rep(0, n), nrow = d, ncol = ns)
  for (j in 1:ns) {
    saisonj <- vec[(d * (j - 1) + 1):(d * j)]
    rangs.vec[order(saisonj), j] <- 1:d
  }
  rangs <- apply(rangs.vec, 1, sum)
  stat <- 12 * sum((rangs - ns * (d + 1)/2)^2)/(ns * d * (d + 1))
  pval <- 1 - pchisq(stat, df = d - 1)
  res <- c(d, d - 1, stat, pval)
  names(res) <- c("P´eriode", "d.f.", "Tobs", "p-valeur")
  return(res)
}

Saison.test(serie,d=10)

## H0) : pas de saisonnalité de période r dans la série
## H1) : présence d'une saisonnalité de période r

### Le p-valeur est 0 qui est inférieur à 0.05, donc on rejette H0, il y a de la saisonnalité de période 12 dans 
### la série.

####3.3.3 On différencie la série pour faire disparaître la composante saisonnière
serie1 <- diff(serie , lag=10)
plot(serie1)

###3.4 On trace de nouveau l'acf pour voir s'il reste une tendance. Si oui, on l'élimine aussi.
acf(serie1)
#### Dans la série1, nous trouvons que la tendance c'est diminue mais,il y en a aussi. donc on va continuer
#### à l'éliminer.

serie2<-( diff(serie1, order=10) )
plot(serie2)

#### Par rapport à la graphique de serie1, l'acf de serie2 c'est autour de 0, donc maintenant nous conclusons
#### qu'il n'y a pas de la tendance dans la série.

####3.4.1 On confime l'existence de la tendance à l'aide des tests des montées (PtMont.test) 
####      ou des discordances (PtDisc.test).

##### le test des montées :
PtMont.test <- function (vec)
{
  n <- length(vec)
  nM <- 0
  for (i in 1:(n-1)){
    if (vec[i] < vec[i+1]) {
      nM <- nM+1         
    }}
  Tobs <- ( nM - (n-1)/2 ) / sqrt( (n+1)/12 )
  p <- 2 * (1- pnorm(abs(Tobs)))
  res <- c(n,nM,Tobs,p)
  names(res) <- c("n","nM","stat"," p-valeur")
  return(res)
}


##### le test des des discordances


PtDisc.test <- function (vec)
{
  n <- length(vec)
  res <- cor.test(1:n, vec, method="kendall")
  nD=floor((1-res$estimate)/4*n*(n-1))
  Tobs <- ( 1 - 4*nD/(n*(n-1)) ) / sqrt( 2*(2*n+5)/(9*n*(n-1)) )
  p <- 2 * (1- pnorm(abs(Tobs)))
  res <- c(n,nD,Tobs,p)
  names(res) <- c("n","nD","stat"," p-valeur")
  return(res)
}

#### Le résultat du test des montées
PtMont.test(serie2)#p-valeur = 0.6152282

#### Le résultat du test des discordances 
PtDisc.test(serie2)#p-valeur = 0.95047987

#### H0) : pas de tendance dans la série
#### H1) : présence d’une tendance

#### Dans les 2 résultat des 2 test, nous trouvons que les p-valeur tout sont supérieur à 0.05, donc on ne rejette 
#### pas H0, donc maintenant, il n'y a pas de tendence dans la série.
acf(serie2)
#### Nous ne pouvons pas trouver la tendance par l'acf de la dernière série.

## 4. On détermine le modèle de la série différenciée.

###4.1 On trace l'acf (acfMA), la pacf et l'eacf pour déterminer les ordres p et q maximaux.

#### Les commandes pour acfMA
acfMA <- function(x) {
  n <- length(x)
  rho <- acf(x, plot = F)$acf[-1]
  nlags <- length(acf(x, plot = F)$lag)
  wh <- c(1, 1 + 2 * cumsum(rho^2))
  seuilsSup <- 1.96 * sqrt(wh/n)
  seuilsSup <- as.ts(seuilsSup)
  seuilsSup <- ts(seuilsSup, start = 1/frequency(x),
                  end = nlags/frequency(x), frequency = frequency(x))
  acf(x)
  lines(seuilsSup, lty = 2, col = "red",
        lwd = 2)
  lines(-seuilsSup, lty = 2, col = "red",
        lwd = 2)
}
acfMA(serie2)
pacf(serie2)
eacf(serie2)
## P= 0 q = 1 





### 4.2 On estime les coefficients du modèle SARIMA(p,d,q)(0,1,0)[r] (Arima).

#### Ici on a utilisé la fonction 'auto.arima()' de la package 'forecast' pour choisir le mieux modèles
modelebien <- auto.arima(serie2);modelebien
#### On trouve que dans notre meilleur modèle c'est ARIMA(4,0,1) qui est avec 0 moyenne.
#### Le résultat nous montre que le meilleur modèles est ARIMA(4,0,1)(1,0,2)[12]

### 4.3 On valide le modèle avec la commande tsdiag(modele),(les résidus doivent former un bruit blanc).
tsdiag(modelebien)
#### on vérifie que les résidus sont bruits blanc par la distributions des résidus.
#### et puis dans l'acf des résidus, on trouve que après lag=0,il n'y a pas de autorégression d'entre eux, c'est à
#### dire que notre modèles est simulé mieux.

acf(resid(modelebien))
### 4.4 On vérifie qu'on ne peut pas simplier encore le modèle.

### 4.5 On effectue des prévisions ( prev(modele) ).

prevmodele <- forecast(modelebien, h = 2)
prevmodele

### 4.6 On vérifie Si les innovations ( modele$res ) sont gaussiennes.
shapiro.test(modelebien$residuals)#p-value = 0.5287
####H0: Les données de l'échantillon ne sont pas significativement différentes de la distribution normale.
####H1: les données de l'échantillon sont très différentes de la distribution normale.

#### P—valeur est 0.5287, qui est supérieur à 0.05, donc on ne rejette pas H0, donc les innovations suivent la 
#### loi gaussinne.

# Test de normalite jarque.bera.test
jarque.bera.test(serie2)#p-value = 0.987

####H0 : les données suivent une loi normale.
####H1 : les données ne suivent pas une loi normale.
#### P—valeur est 0.987, qui est supérieur à 0.05, donc on ne rejette pas H0, donc les innovations suivent la 
#### loi gaussinne.

## 5. Estimer les prévisions à l’ordre 1 et 2, i.e., pour X201,X202(avec l'intervalle confiance)
prediction <-forecast(modelebien,h=2,level=c(99.5))
prediction
plot(prediction)

### On prédit que pour X201, le point Forecast = 2.3435082
###               pour X202, le point Forecast = -0.3511171.(L'intervall confiance est 5%)

