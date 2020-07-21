## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----setup python--------------------------------------------------------
library(reticulate)

# Import minepy from python
minepy <- import("minepy")
mine = minepy$MINE(alpha=0.6, c=15, est="mic_approx")


## ------------------------------------------------------------------------
library(MASS)
library(miic)
library(dplyr)
library(ggplot2)

# Run parameters
N=1000
nreps=10000

# Distribution parameters
mu = c(0, 0) # Means
s1 = 1 # Standard deviations
s2 = 1

miic_res_nbins = data.frame(rho=numeric(), ncutsx=numeric(), ncutsy=numeric())
mine_res_nbins = data.frame(rho=numeric(), ncutsx=numeric(), ncutsy=numeric())

for(i in 1:nreps){
  if((100*i/nreps)%%10==0){print(paste0(100*i/nreps,"%"))}

  rho = runif(1, -1, 1)
  sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2) # Covariance matrix
  # Bivariate normal with N samples
  bvn <- mvrnorm(N, mu = mu, Sigma = sigma )

  # miic
  res = miic::discretizeMutual(bvn[,1], bvn[,2], plot=F)
  miic_res_nbins[nrow(miic_res_nbins)+1,] = list(rho, length(res$cutpoints1), length(res$cutpoints2))
  
  # MINE
  mine$compute_score(bvn[,1], bvn[,2])
  mine_matrix = mine$get_score()
  mic = mine$mic()

  cuts_X = which(unlist(lapply(lapply(mine_matrix, function(x) x == mic), any)))[1]
  cuts_Y = which(mine_matrix[[cuts_X]] == mic)[1]
  
  mine_res_nbins[nrow(mine_res_nbins)+1,] = list(rho, cuts_X+1, cuts_Y+1)
}


## ----plot miic, fig.width = 6, fig.height = 6----------------------------
miic_res_nbins %>%
  group_by(ncutsx, ncutsy) %>%
  summarise(avgrho = mean(abs(rho)), avgRI = mean(-log(1-rho**2)/2), count = n()) %>%
  ggplot(aes(x=ncutsx-1, y=ncutsy-1)) +
    geom_point(aes(fill=avgRI, size=count), shape=22) +
    scale_fill_gradient(low = 'white', high = '#0013a3', limits=c(0,3)) +
    scale_x_continuous(breaks = 1:50, minor_breaks = c()) +
    scale_y_continuous(breaks = 1:50, minor_breaks = c()) +
    coord_fixed(ratio = 1, xlim = c(1,max(c(miic_res_nbins$ncutsx, miic_res_nbins$ncutsy))-1),
                ylim = c(1,max(c(miic_res_nbins$ncutsx, miic_res_nbins$ncutsy))-1)) +
    scale_size_continuous(range = c(2,10), limits = c(0,nreps)) +
    theme_minimal() +
    ylab("# bins on Y") +
    xlab("# bins on X")



## ----plot mine, fig.width = 6, fig.height = 6----------------------------
mine_res_nbins %>%
  group_by(ncutsx, ncutsy) %>%
  summarise(avgrho = mean(abs(rho)), avgRI = mean(-log(1-rho**2)/2), count = n()) %>%
  ggplot(aes(x=ncutsx, y=ncutsy)) +
    geom_point(aes(fill=avgRI, size=count), shape=22) +
    scale_fill_gradient(low = 'white', high = '#0013a3', limits=c(0,3)) +
    scale_x_continuous(breaks = 1:50, minor_breaks = c()) +
    scale_y_continuous(breaks = 1:50, minor_breaks = c()) +
    #coord_fixed(ratio = 1, xlim = c(1,max(c(mine_res_nbins$ncutsx, mine_res_nbins$ncutsy))-1),
    #            ylim = c(1,max(c(mine_res_nbins$ncutsx, mine_res_nbins$ncutsy))-1)) +
    scale_size_continuous(range = c(2,10), limits = c(0,nreps)) +
    theme_minimal() +
    ylab("# bins on Y") +
    xlab("# bins on X")

