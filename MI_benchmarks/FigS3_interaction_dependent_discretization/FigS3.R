## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- fig.height=6, fig.width=6------------------------------------------
library(miic)

N = 1000

X = runif(N, min = 0, max=0.3)
noise = c(rnorm(N, sd=0.4))
Y1 = log(X) + noise
miic::discretizeMutual(X, Y1)$infok
Y2 = X**5 + noise/1900
miic::discretizeMutual(X, Y2)$infok

