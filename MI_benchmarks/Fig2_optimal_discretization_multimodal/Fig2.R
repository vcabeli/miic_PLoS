## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----imports, message=FALSE----------------------------------------------
library(MASS)
library(copula)
library(distr)
library(miic)


## ----marginals-----------------------------------------------------------
mixnorm1 <- UnivarMixingDistribution(Norm(0,1), Norm(4,2), Norm(6,0.7))
dmixnorm1 <- d(mixnorm1); qmixnorm1 <- q(mixnorm1); pmixnorm1 <- p(mixnorm1)


mixnorm2 <- UnivarMixingDistribution(Norm(-3,2), Norm(6,0.5), Norm(9,0.5))
dmixnorm2 <- d(mixnorm2); qmixnorm2 <- q(mixnorm2); pmixnorm2 <- p(mixnorm2)


## ----copula--------------------------------------------------------------
mv.NE <- mvdc(gumbelCopula(5), margins=c("mixnorm1","mixnorm2"), paramMargins=list(list(),list()))


## ----plot, fig.width = 6, fig.height = 6---------------------------------
N = 500
data = rMvdc(N, mv.NE)
print(discretizeMutual(data[,1], data[,2])$info)

