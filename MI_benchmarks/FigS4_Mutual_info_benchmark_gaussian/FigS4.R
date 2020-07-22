## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----imports-------------------------------------------------------------
library(MASS)
library(ggplot2)
library(dplyr)
library(miic)
library(doParallel) # Multithreading
library(doSNOW) # Multithreading

################################################################################
# KNN estimation methods with jidt
library("rJava")
.jinit()
.jaddClassPath("/home/vcabeli/Downloads/jidt/infodynamics.jar")

jidt_compute_mi = function(X,Y, k=3){
  miCalc<-.jnew("infodynamics/measures/continuous/kraskov/MutualInfoCalculatorMultiVariateKraskov1")
  .jcall(miCalc,"V","setProperty", "k", as.character(k))
  .jcall(miCalc,"V","initialise")
  .jcall(miCalc, "V", "setObservations", X, Y)

  return(.jcall(miCalc,"D","computeAverageLocalOfObservations"))
}
################################################################################


## ----MI function---------------------------------------------------------
# Estimate MI with all methods for a given N, rho
gaussian_MI_estimation = function(N, rho, rep){

  results = data.frame(N=numeric(), rho=numeric(), method=character(), MI=numeric(), rep=numeric(), stringsAsFactors = F)

  # Distribution parameters
  mu = c(0, 0)
  s1 = 1
  s2 = 1
  # Covariance matrix
  sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2)
  # Bivariate normal distribution
  bvn <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS

  res = miic::discretizeMutual(bvn[,1], bvn[,2], plot=F, cplx="nml", maxbins=N)
  results[nrow(results)+1,] = list(N, rho, "miic", res$info, rep)
  results[nrow(results)+1,] = list(N, rho, "knn3", jidt_compute_mi(bvn[,1], bvn[,2], k = 3), rep)
  results[nrow(results)+1,] = list(N, rho, "knn10", jidt_compute_mi(bvn[,1], bvn[,2], k = 10), rep)
  results[nrow(results)+1,] = list(N, rho, "knn20", jidt_compute_mi(bvn[,1], bvn[,2], k = 20), rep)
  results
}


## ----run benchmarks------------------------------------------------------
#############
# Run all settings
Nmin=100; Nmax=10000; steps=7; nreps=100
Nvals = ceiling(10**(seq(log10(Nmin), log10(Nmax), length.out = steps)))
rhos = c(10**(seq(-2,-0.4,length.out=10)), seq(0.5, 0.9, by=0.13))

# Set up parallel backend
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
registerDoSNOW(cl)

# Set up progress bar
pb <- txtProgressBar(max = nreps*length(Nvals)*length(rhos), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Launch parallel computation
results <- foreach(rep=1:nreps, .combine=rbind, .inorder = F,
                   .packages = c("miic", "MASS", "rJava")) %:%
  foreach(rho=rhos, .combine=rbind, .inorder = F) %:%
    foreach(N=Nvals, .combine=rbind, .inorder = F, .options.snow=opts) %dopar% {
      .jinit()
      .jaddClassPath("/home/vcabeli/Downloads/jidt/infodynamics.jar")
      gaussian_MI_estimation(N, rho, rep)
    }
# Close backend
close(pb)
stopCluster(cl)


## ----plot, fig.width=4, fig.height=9-------------------------------------
##########
# Plot results
estimation = results %>% group_by(method,rho,N) %>% summarize(mean_MI = mean(MI), sd = sd(MI)) %>%
    ggplot(aes(x=rho, y=mean_MI, color=method)) +
    geom_line() +
    geom_errorbar(aes(ymin=mean_MI-sd, ymax=mean_MI+sd), alpha=0.5) +
    ylab("Estimated I(X;Y)") +
    scale_color_discrete(drop=TRUE, limits = levels(results$method)) +
    theme_bw() +
    facet_grid(N ~ .) +
    theme(legend.position = "top")
plot(estimation)

error = results %>% group_by(method,rho,N) %>% mutate(error = MI-(-log(1-rho**2)/2)) %>% summarize(MSE = mean(error**2)) %>%
    ggplot(aes(x=rho, y=MSE, color=method)) +
    geom_line() +
    scale_color_discrete(drop=TRUE, limits = levels(results$method)) +
    scale_y_continuous(trans='log10') +
    theme_bw() +
    facet_grid(N ~ .)+#, scales = "free_y") +
    theme(legend.position = "top")
plot(error)

variance = results %>% group_by(method,rho,N) %>% summarize(variance = sd(MI)) %>%
    ggplot(aes(x=rho, y=variance, color=method)) +
    geom_line() +
    ylab("Estimator variance") +
    scale_color_discrete(drop=TRUE, limits = levels(results$method)) +
    scale_y_continuous(trans='log10') +
    theme_bw() +
    facet_grid(N ~ .)+#, scales="free_y") +
    theme(legend.position = "top")
plot(variance)

