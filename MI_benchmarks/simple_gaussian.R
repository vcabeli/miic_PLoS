library(MASS)
library(ggplot2)
library(dplyr)
library(miic)
library(doParallel)
library(doSNOW)

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


#Estimate MI with all methods for a given N, rho
gaussian_MI_estimation = function(N, rho, rep){
  results = data.frame(N=numeric(), rho=numeric(), method=character(), MI=numeric(), rep=numeric(), stringsAsFactors = F)
    
  mu1 = 1
  mu2 = 10
  mu = c(mu1, mu2)
  s1 = 2
  s2 = 5
  sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2) # Covariance matrix


  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS

  res = miic::discretizeMutual(bvn1[,1], bvn1[,2], plot=F, cplx="nml", maxbins=N)
  results[nrow(results)+1,] = list(N, rho, "miic", res$info, rep)
  results[nrow(results)+1,] = list(N, rho, "knn3", jidt_compute_mi(bvn1[,1], bvn1[,2], k = 3), rep)
  results[nrow(results)+1,] = list(N, rho, "knn10", jidt_compute_mi(bvn1[,1], bvn1[,2], k = 10), rep)
  results[nrow(results)+1,] = list(N, rho, "knn20", jidt_compute_mi(bvn1[,1], bvn1[,2], k = 20), rep)
  #results[nrow(results)+1,] = list(rho, "analytical", -log(1-rho**2)/2)
  results
}



#############
# Run all settings
Nmin=100; Nmax=10000; steps=7; reps=20
#Nvals = ceiling(10**(seq(log10(Nmin), log10(Nmax), length.out = steps)))
Nvals = c(100, 200, 1000, 10000)
#rhos = seq(-0.95, 0.95, by=0.05)
rhos = c(10**(seq(-2,-0.4,length.out=10)), seq(0.5, 0.9, by=0.13))

# Set up parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
registerDoSNOW(cl)

# Set up progress bar
pb <- txtProgressBar(max = reps*length(Nvals)*length(rhos), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Launch parallel computation
results <- foreach(rep=1:reps, .combine=rbind, .inorder = F,
                   .packages = c("miic", "MASS", "rJava")) %:% 
  foreach(rho=rhos, .combine=rbind, .inorder = F) %:%
    foreach(N=Nvals, .combine=rbind, .inorder = F, .options.snow=opts) %dopar% {
      .jinit()
      .jaddClassPath("/home/vcabeli/Downloads/jidt/infodynamics.jar")
      gaussian_MI_estimation(N, rho, rep)
}
close(pb)
stopCluster(cl)



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
ggsave("simple_gaussian_estimation_bis.pdf", device = "pdf", plot = estimation,
       width = 3, height=8, units = "in")

error = results %>% group_by(method,rho,N) %>% mutate(error = MI-(-log(1-rho**2)/2)) %>% summarize(MSE = mean(error**2)) %>%
    ggplot(aes(x=rho, y=MSE, color=method)) +
    geom_line() +
    scale_color_discrete(drop=TRUE, limits = levels(results$method)) +
    scale_y_continuous(trans='log10') +
    theme_bw() + 
    facet_grid(N ~ .)+#, scales = "free_y") +
    theme(legend.position = "top")
plot(error)
ggsave("simple_gaussian_error_bis.pdf", device = "pdf", plot = error,
       width = 3, height=8, units = "in")

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
ggsave("simple_gaussian_variance_bis.pdf", device = "pdf", plot = variance,
       width = 3, height=8, units = "in")




source("~/work/projects/miic_dyn_programming_MI/simple_gaussian_plots.R")

simple_gaussian_plot(results)

results_2 = results %>% group_by(rho, method,N) %>% mutate(error = (-log(1-rho**2)/2-MI)**2)
simple_gaussian_plot_error(results_2)

results_2 = results %>% group_by(rho, method,N) %>% mutate(error = -log(1-rho**2)/2-MI) %>% summarize(mean=mean(MI), MISD=sd(MI))
simple_gaussian_plot_SD(results_2)

results_2 = results %>% filter(method %in% c("analytical", "knn3", "knn10", "knn20", "miic"))
#results_2 = results %>% filter(method %in% c("analytical", "miic", "miick", "ef_max"))
results_2$method = factor(results_2$method)

results_2 %>% group_by(rho, method) %>%
  summarize(mean=median(MI), sd=sd(MI)) %>%
  ggplot() +
    geom_line(aes(x=rho, y=mean, group=method, color=method)) +
    geom_ribbon(aes(x=rho, ymin=mean-sd, ymax=mean+sd, group=method, fill=method), alpha = 0.1) +
    scale_color_discrete(drop=TRUE, limits = levels(results_2$method)) +
    ylab("estimated MI")

results_2 %>% group_by(N, rho, method) %>%
  summarize(mean=mean(MI), sd=sd(MI)) %>%
  ggplot()+
    geom_line(aes(x=rho, y=sd, group=method, color=method)) +
    ylab("SD") +
    scale_color_discrete(drop=TRUE, limits = levels(results_2$method)) +
    ylab("estimation standard deviation")

results_3 = results_2 %>% filter(method != 'analytical')





g= results_3 %>% filter(rho!=0) %>% group_by(method,rho) %>% mutate(error = MI/(-log(1-rho**2)/2)) %>% summarize(mean_error = mean(error), sd = sd(error)) %>%
    ggplot(aes(x=rho, y=mean_error, color=method)) +
    geom_hline(yintercept = 1, lty='dashed', alpha=0.5) +
    geom_line() +
    geom_errorbar(aes(ymin=mean_error-sd, ymax=mean_error+sd), alpha=0.5) +
    ylab("estimated MI / exact MI") +
    scale_color_discrete(drop=TRUE, limits = levels(results_2$method)) #+
plot(g)





hdfx = function(x) dnorm(x, mean = mu1, sd = s1)

pdfy = function(y) dnorm(y, mean = mu2, sd = s2)
pdfxcondy = function(x,y) dnorm(x, mu1 + (s1/s2)*rho*(y-mu2), (1-rho**2)*s1**2)
pdfycondx = function(x,y) dnorm(y, mu2 + (s2/s1)*rho*(x-mu1), (1-rho**2)*s2**2)
pdfxy = function(x,y) 1/(2*pi*s1*s2*sqrt(1-rho**2)) * exp(
                      -(1/(2*(1-rho**2)) * ( (x-mu1)**2/s1**2 + (y-mu2)**2/s2**2 - 2*rho*(x-mu1)*(y-mu2)/(s1*s2) )))


HX  = -integrate(function(x) pdfx(x)*log(pdfx(x)), lower = -30, upper = 30)$value
HY  = -integrate(function(y) pdfy(y)*log(pdfy(y)), lower = -30, upper = 30)$value

InnerIntegral = Vectorize(function(y) { integrate(function(x) pdfxy(x, y)*log(pdfxy(x,y)), -30, y)$value})
HXY = -integrate(InnerIntegral , -30, 30)$value
InnerIntegral = Vectorize(function(x) { integrate(function(y) pdfxy(x, y)*log(pdfxy(x,y)), -30, x)$value})
HXY = HXY -integrate(InnerIntegral , -30, 30)$value

Cm = matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2) # Covariance matrix
HXY2 = 1/2*log((2*pi*exp(1))**2*det(Cm))

# Compute information
HX + HY - HXY
HX + HY - HXY2

-log(1-rho**2)/2



pdfxy2 = function(x) pdfxy(x[1], x[2])
cubature::adaptIntegrate(function(x) -pdfxy2(x) * log(pdfxy2(x)), lowerLimit = c(-30,-30), upperLimit = c(30,30))$integral


pdfmix = function(x) 0.5*pdfx(x)+0.5*pdfy(x)
plot(curve(pdfmix, from = -10, to = 30, n = 200))

