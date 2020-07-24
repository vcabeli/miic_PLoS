FigS6 Conditional mutual information estimation for multivariate
Gaussian distributions
================



This script aims to benchmark the miic conditional mutual information
estimation for Gaussian multivariate distributions. Four-dimensional
normal distributions
![P(X,Y,Z\_1,Z\_2)](https://latex.codecogs.com/png.latex?P%28X%2CY%2CZ_1%2CZ_2%29
"P(X,Y,Z_1,Z_2)") are sampled for
![N=100](https://latex.codecogs.com/png.latex?N%3D100 "N=100") to
![5,000](https://latex.codecogs.com/png.latex?5%2C000 "5,000") samples
![100](https://latex.codecogs.com/png.latex?100 "100") times for each
correlation coefficient
![\\rho=\\rho\_{XY}](https://latex.codecogs.com/png.latex?%5Crho%3D%5Crho_%7BXY%7D
"\\rho=\\rho_{XY}") between
![0.05](https://latex.codecogs.com/png.latex?0.05 "0.05") and
![0.95](https://latex.codecogs.com/png.latex?0.95 "0.95"). The other
pairwise correlation coefficients are fixed as
![\\rho\_{XZ\_1}=\\rho\_{XZ\_2}=\\rho\_{YZ\_1}=\\rho\_{YZ\_2}=\\lambda=0.7](https://latex.codecogs.com/png.latex?%5Crho_%7BXZ_1%7D%3D%5Crho_%7BXZ_2%7D%3D%5Crho_%7BYZ_1%7D%3D%5Crho_%7BYZ_2%7D%3D%5Clambda%3D0.7
"\\rho_{XZ_1}=\\rho_{XZ_2}=\\rho_{YZ_1}=\\rho_{YZ_2}=\\lambda=0.7") and
![\\rho\_{Z\_1Z\_2}=0.9](https://latex.codecogs.com/png.latex?%5Crho_%7BZ_1Z_2%7D%3D0.9
"\\rho_{Z_1Z_2}=0.9"). The conditional mutual information
![I(X;Y|Z\_1,Z\_2)](https://latex.codecogs.com/png.latex?I%28X%3BY%7CZ_1%2CZ_2%29
"I(X;Y|Z_1,Z_2)") was then estimated using the proposed optimum
partitioning scheme as well as with kNN conditional information
estimates as in Fig. S4.

To run this script you will need to download the
[jidt](https://github.com/jlizier/jidt) library (Lizier 2014) which
contains many mutual information estimators, and install `rJava` to be
able to call Java function from R. If you want to use multithreading,
you will also need to install `doParallel` and `doSNOW`.

``` r
################################################################################
# KNN estimation methods with jidt
library("rJava")
.jinit()
.jaddClassPath("/home/vcabeli/Downloads/jidt/infodynamics.jar")


jidt_compute_cmi = function(X,Y,Z, k=3){
  cmiCalc<-.jnew("infodynamics/measures/continuous/kraskov/ConditionalMutualInfoCalculatorMultiVariateKraskov1")
  .jcall(cmiCalc,"V","setProperty", "k", as.character(k))
  .jcall(cmiCalc,"V","initialise", as.integer(1), as.integer(1), ncol(Z))
  .jcall(cmiCalc, "V", "setObservations", .jarray(matrix(X, ncol=1), dispatch = T),
                                          .jarray(matrix(Y, ncol=1), dispatch = T), 
                                          .jarray(Z, dispatch = T))
  
  return(.jcall(cmiCalc,"D","computeAverageLocalOfObservations"))
}
```

![\\rho](https://latex.codecogs.com/png.latex?%5Crho "\\rho") values
closed to zero, mimick \`\`V-structures’’ as they correspond to pairwise
independence but conditional dependence; by constrast ![\\rho = 2
\\lambda^2 / (1+\\rho\_{Z\_1Z\_2})
\\simeq 0.5158](https://latex.codecogs.com/png.latex?%5Crho%20%3D%202%20%5Clambda%5E2%20%2F%20%281%2B%5Crho_%7BZ_1Z_2%7D%29%20%5Csimeq%200.5158
"\\rho = 2 \\lambda^2 / (1+\\rho_{Z_1Z_2}) \\simeq 0.5158") corresponds
to conditional independence, while ![\\rho
\> 0.5158](https://latex.codecogs.com/png.latex?%5Crho%20%3E%200.5158
"\\rho \> 0.5158") implies that
![X](https://latex.codecogs.com/png.latex?X "X") and
![Y](https://latex.codecogs.com/png.latex?Y "Y") share more information
than the indirect flow through
![Z\_1](https://latex.codecogs.com/png.latex?Z_1 "Z_1") and
![Z\_2](https://latex.codecogs.com/png.latex?Z_2 "Z_2").

![Benchmark setup](FigS6_files/FigS6.png)

``` r
# Distribution parameters
mu=c(0,0,0,0)
s1=1; s2=1; s3=1; s4=1
rho2=0.7
rho3=0.7
rho4=0.9

# Estimate CMI with all methods for a given N, rho
gaussian_CMI_estimation = function(N, rho, rep, noise_Z){
  results = data.frame(N=numeric(), rho=numeric(), method=character(), MI=numeric(), rep=numeric(), stringsAsFactors = F)
    
  sigma=matrix(c(
    s1^2,       rho*s1*s2,  rho2*s1*s3, rho2*s1*s4,
    rho*s1*s2,  s2^2,       rho3*s2*s3, rho3*s2*s4,
    rho2*s1*s3, rho3*s2*s3, s3^2,       rho4*s3*s4,
    rho2*s1*s4, rho3*s2*s4, rho4*s3*s4, s4^2),
    nrow=4,ncol=4,byrow=TRUE);
  
  data=mvrnorm(n = N, mu, sigma, tol = 1e-2, empirical = TRUE, EISPACK = FALSE);

  X = rank(data[,1])
  Y = rank(data[,2])
  Z = data[,c(3,4)]
  
  if(noise_Z > 0){
    Z = cbind(Z, matrix(rnorm(N*noise_Z), ncol = noise_Z))
    Z = Z[,sample(1:(2+noise_Z))] # shuffle column order
  }
  Z = apply(Z, 2, rank)
  
  res = miic::discretizeMutual(X, Y, matrix_u = Z, plot=F)
  results[nrow(results)+1,] = list(N, rho, "miic", res$info, rep)
  results[nrow(results)+1,] = list(N, rho, "miic_k", res$infok, rep)
  results[nrow(results)+1,] = list(N, rho, "knn3",  jidt_compute_cmi(X, Y, Z, k = 3),  rep)
  results[nrow(results)+1,] = list(N, rho, "knn10", jidt_compute_cmi(X, Y, Z, k = 10), rep)
  results[nrow(results)+1,] = list(N, rho, "knn20", jidt_compute_cmi(X, Y, Z, k = 20), rep)
  results
}
```

``` r
#############
# Run all settings
Nvals = c(100, 200, 1000, 10000)
cond_indep_rho = 0.5158
rhos = c(cond_indep_rho - 10**seq(-2.8,log10(cond_indep_rho-0.3),length.out=6),
         cond_indep_rho + 10**seq(-2.8,log10(0.7-cond_indep_rho),length.out=6))
rhos = sort(c(rhos, 0.05, 0.2, 0.8, 0.95))
noise_Z = 0
reps = 100

# Set up parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
registerDoSNOW(cl)

# Set up progress bar
pb <- txtProgressBar(max = reps*length(Nvals)*length(rhos), style = 3)
```

    ##   |                                                                         |                                                                 |   0%

``` r
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Launch parallel computation
results <- foreach(rep=1:reps, .combine=rbind, .inorder = F,
                   .packages = c("miic", "MASS", "rJava")) %:% 
  foreach(rho=rhos, .combine=rbind, .inorder = F) %:%
    foreach(N=Nvals, .combine=rbind, .inorder = F, .options.snow=opts) %dopar% {
      source("~/work/projects/conditional_independence_testing/jidt_functions.R")
      gaussian_CMI_estimation(N, rho, rep, noise_Z)
}
```

    ##   |                                                                         |                                                                 |   1%  |                                                                         |=                                                                |   1%  |                                                                         |=                                                                |   2%  |                                                                         |==                                                               |   2%  |                                                                         |==                                                               |   3%  |                                                                         |==                                                               |   4%  |                                                                         |===                                                              |   4%  |                                                                         |===                                                              |   5%  |                                                                         |====                                                             |   5%  |                                                                         |====                                                             |   6%  |                                                                         |====                                                             |   7%  |                                                                         |=====                                                            |   7%  |                                                                         |=====                                                            |   8%  |                                                                         |======                                                           |   8%  |                                                                         |======                                                           |   9%  |                                                                         |======                                                           |  10%  |                                                                         |=======                                                          |  10%  |                                                                         |=======                                                          |  11%  |                                                                         |=======                                                          |  12%  |                                                                         |========                                                         |  12%  |                                                                         |========                                                         |  13%  |                                                                         |=========                                                        |  13%  |                                                                         |=========                                                        |  14%  |                                                                         |=========                                                        |  15%  |                                                                         |==========                                                       |  15%  |                                                                         |==========                                                       |  16%  |                                                                         |===========                                                      |  16%  |                                                                         |===========                                                      |  17%  |                                                                         |===========                                                      |  18%  |                                                                         |============                                                     |  18%  |                                                                         |============                                                     |  19%  |                                                                         |=============                                                    |  19%  |                                                                         |=============                                                    |  20%  |                                                                         |=============                                                    |  21%  |                                                                         |==============                                                   |  21%  |                                                                         |==============                                                   |  22%  |                                                                         |===============                                                  |  22%  |                                                                         |===============                                                  |  23%  |                                                                         |===============                                                  |  24%  |                                                                         |================                                                 |  24%  |                                                                         |================                                                 |  25%  |                                                                         |=================                                                |  25%  |                                                                         |=================                                                |  26%  |                                                                         |=================                                                |  27%  |                                                                         |==================                                               |  27%  |                                                                         |==================                                               |  28%  |                                                                         |===================                                              |  28%  |                                                                         |===================                                              |  29%  |                                                                         |===================                                              |  30%  |                                                                         |====================                                             |  30%  |                                                                         |====================                                             |  31%  |                                                                         |====================                                             |  32%  |                                                                         |=====================                                            |  32%  |                                                                         |=====================                                            |  33%  |                                                                         |======================                                           |  33%  |                                                                         |======================                                           |  34%  |                                                                         |======================                                           |  35%  |                                                                         |=======================                                          |  35%  |                                                                         |=======================                                          |  36%  |                                                                         |========================                                         |  36%  |                                                                         |========================                                         |  37%  |                                                                         |========================                                         |  38%  |                                                                         |=========================                                        |  38%  |                                                                         |=========================                                        |  39%  |                                                                         |==========================                                       |  39%  |                                                                         |==========================                                       |  40%  |                                                                         |==========================                                       |  41%  |                                                                         |===========================                                      |  41%  |                                                                         |===========================                                      |  42%  |                                                                         |============================                                     |  42%  |                                                                         |============================                                     |  43%  |                                                                         |============================                                     |  44%  |                                                                         |=============================                                    |  44%  |                                                                         |=============================                                    |  45%  |                                                                         |==============================                                   |  45%  |                                                                         |==============================                                   |  46%  |                                                                         |==============================                                   |  47%  |                                                                         |===============================                                  |  47%  |                                                                         |===============================                                  |  48%  |                                                                         |================================                                 |  48%  |                                                                         |================================                                 |  49%  |                                                                         |================================                                 |  50%  |                                                                         |=================================                                |  50%  |                                                                         |=================================                                |  51%  |                                                                         |=================================                                |  52%  |                                                                         |==================================                               |  52%  |                                                                         |==================================                               |  53%  |                                                                         |===================================                              |  53%  |                                                                         |===================================                              |  54%  |                                                                         |===================================                              |  55%  |                                                                         |====================================                             |  55%  |                                                                         |====================================                             |  56%  |                                                                         |=====================================                            |  56%  |                                                                         |=====================================                            |  57%  |                                                                         |=====================================                            |  58%  |                                                                         |======================================                           |  58%  |                                                                         |======================================                           |  59%  |                                                                         |=======================================                          |  59%  |                                                                         |=======================================                          |  60%  |                                                                         |=======================================                          |  61%  |                                                                         |========================================                         |  61%  |                                                                         |========================================                         |  62%  |                                                                         |=========================================                        |  62%  |                                                                         |=========================================                        |  63%  |                                                                         |=========================================                        |  64%  |                                                                         |==========================================                       |  64%  |                                                                         |==========================================                       |  65%  |                                                                         |===========================================                      |  65%  |                                                                         |===========================================                      |  66%  |                                                                         |===========================================                      |  67%  |                                                                         |============================================                     |  67%  |                                                                         |============================================                     |  68%  |                                                                         |=============================================                    |  68%  |                                                                         |=============================================                    |  69%  |                                                                         |=============================================                    |  70%  |                                                                         |==============================================                   |  70%  |                                                                         |==============================================                   |  71%  |                                                                         |==============================================                   |  72%  |                                                                         |===============================================                  |  72%  |                                                                         |===============================================                  |  73%  |                                                                         |================================================                 |  73%  |                                                                         |================================================                 |  74%  |                                                                         |================================================                 |  75%  |                                                                         |=================================================                |  75%  |                                                                         |=================================================                |  76%  |                                                                         |==================================================               |  76%  |                                                                         |==================================================               |  77%  |                                                                         |==================================================               |  78%  |                                                                         |===================================================              |  78%  |                                                                         |===================================================              |  79%  |                                                                         |====================================================             |  79%  |                                                                         |====================================================             |  80%  |                                                                         |====================================================             |  81%  |                                                                         |=====================================================            |  81%  |                                                                         |=====================================================            |  82%  |                                                                         |======================================================           |  82%  |                                                                         |======================================================           |  83%  |                                                                         |======================================================           |  84%  |                                                                         |=======================================================          |  84%  |                                                                         |=======================================================          |  85%  |                                                                         |========================================================         |  85%  |                                                                         |========================================================         |  86%  |                                                                         |========================================================         |  87%  |                                                                         |=========================================================        |  87%  |                                                                         |=========================================================        |  88%  |                                                                         |==========================================================       |  88%  |                                                                         |==========================================================       |  89%  |                                                                         |==========================================================       |  90%  |                                                                         |===========================================================      |  90%  |                                                                         |===========================================================      |  91%  |                                                                         |===========================================================      |  92%  |                                                                         |============================================================     |  92%  |                                                                         |============================================================     |  93%  |                                                                         |=============================================================    |  93%  |                                                                         |=============================================================    |  94%  |                                                                         |=============================================================    |  95%  |                                                                         |==============================================================   |  95%  |                                                                         |==============================================================   |  96%  |                                                                         |===============================================================  |  96%  |                                                                         |===============================================================  |  97%  |                                                                         |===============================================================  |  98%  |                                                                         |================================================================ |  98%  |                                                                         |================================================================ |  99%  |                                                                         |=================================================================|  99%  |                                                                         |=================================================================| 100%

``` r
close(pb)
```

``` r
stopCluster(cl)
```

The analytical value of the conditional mutual information is derived as
follows; given the ![4
\\times 4](https://latex.codecogs.com/png.latex?4%20%5Ctimes%204
"4 \\times 4") covariance matrix
![{\\boldsymbol{\\Sigma}}](https://latex.codecogs.com/png.latex?%7B%5Cboldsymbol%7B%5CSigma%7D%7D
"{\\boldsymbol{\\Sigma}}") and its four ![2
\\times 2](https://latex.codecogs.com/png.latex?2%20%5Ctimes%202
"2 \\times 2") partitions
![{\\boldsymbol{\\Sigma}}\_{ij}](https://latex.codecogs.com/png.latex?%7B%5Cboldsymbol%7B%5CSigma%7D%7D_%7Bij%7D
"{\\boldsymbol{\\Sigma}}_{ij}"), we first compute the conditional
covariance matrix ![\\overline{\\boldsymbol\\Sigma} =
\\boldsymbol\\Sigma\_{11} - \\boldsymbol\\Sigma\_{12}
\\boldsymbol\\Sigma\_{22}^{-1}
\\boldsymbol\\Sigma\_{21}](https://latex.codecogs.com/png.latex?%5Coverline%7B%5Cboldsymbol%5CSigma%7D%20%3D%20%5Cboldsymbol%5CSigma_%7B11%7D%20-%20%5Cboldsymbol%5CSigma_%7B12%7D%20%5Cboldsymbol%5CSigma_%7B22%7D%5E%7B-1%7D%20%5Cboldsymbol%5CSigma_%7B21%7D
"\\overline{\\boldsymbol\\Sigma} = \\boldsymbol\\Sigma_{11} - \\boldsymbol\\Sigma_{12} \\boldsymbol\\Sigma_{22}^{-1} \\boldsymbol\\Sigma_{21}")
where
![\\boldsymbol\\Sigma\_{22}^{-1}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%5CSigma_%7B22%7D%5E%7B-1%7D
"\\boldsymbol\\Sigma_{22}^{-1}") is the generalized inverse of
![\\boldsymbol\\Sigma\_{22}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%5CSigma_%7B22%7D
"\\boldsymbol\\Sigma_{22}"). The partial correlation between
![X](https://latex.codecogs.com/png.latex?X "X") and
![Y](https://latex.codecogs.com/png.latex?Y "Y") is obtained as
![\\rho\_{XY\\cdot Z\_1Z\_2} = {\\overline{\\boldsymbol\\Sigma}\_{12}
/{\\sqrt{\\overline{\\boldsymbol\\Sigma}\_{11} \*
\\overline{\\boldsymbol\\Sigma}\_{22}}}
}](https://latex.codecogs.com/png.latex?%5Crho_%7BXY%5Ccdot%20Z_1Z_2%7D%20%3D%20%7B%5Coverline%7B%5Cboldsymbol%5CSigma%7D_%7B12%7D%20%2F%7B%5Csqrt%7B%5Coverline%7B%5Cboldsymbol%5CSigma%7D_%7B11%7D%20%2A%20%5Coverline%7B%5Cboldsymbol%5CSigma%7D_%7B22%7D%7D%7D%20%7D
"\\rho_{XY\\cdot Z_1Z_2} = {\\overline{\\boldsymbol\\Sigma}_{12} /{\\sqrt{\\overline{\\boldsymbol\\Sigma}_{11} * \\overline{\\boldsymbol\\Sigma}_{22}}} }"),
and the analytical conditional mutual information for a multivariate
normal distribution is given by
![I(X;Y|Z\_1,Z\_2)=-\\log(1-\\rho\_{XY\\cdot
Z\_1Z\_2}^2)/2](https://latex.codecogs.com/png.latex?I%28X%3BY%7CZ_1%2CZ_2%29%3D-%5Clog%281-%5Crho_%7BXY%5Ccdot%20Z_1Z_2%7D%5E2%29%2F2
"I(X;Y|Z_1,Z_2)=-\\log(1-\\rho_{XY\\cdot Z_1Z_2}^2)/2").

``` r
###########################################
# Plots
get_true_CMI = function(rho){
  rho = unique(rho)
  #mu=c(0,0,0,0)
  #s1=1; s2=1; s3=1; s4=1
  #rho2=0.7
  #rho3=0.7
  #rho4=0.9
  sigma=matrix(c(
    s1^2,       rho*s1*s2,  rho2*s1*s3, rho2*s1*s4,
    rho*s1*s2,  s2^2,       rho3*s2*s3, rho3*s2*s4,
    rho2*s1*s3, rho3*s2*s3, s3^2,       rho4*s3*s4,
    rho2*s1*s4, rho3*s2*s4, rho4*s3*s4, s4^2),
    nrow=4,ncol=4,byrow=TRUE);
  
  sigma11=sigma[1:2,1:2]
  sigma22=sigma[3:4,3:4]
  sigma12=sigma[1:2,3:4]
  sigma21=sigma[3:4,1:2]
  
  sigmac=sigma11-sigma12%*%ginv(sigma22)%*%sigma21
  #sigmac=0.5*(sigmac+t(sigmac))
  
  rr=sigmac[1,2]/sqrt(sigmac[1,1])/sqrt(sigmac[2,2])
  MItrue=-0.5*log(1-rr^2)
  
  MItrue
}


estimation = results %>% filter(method %in% c("miic", "knn3", "knn10", "knn20")) %>%
                         group_by(method,rho, N) %>%
                         summarize(mean_MI = mean(MI), sd = sd(MI)) %>%
    ggplot(aes(x=rho, y=mean_MI, color=method)) +
    geom_line() +
    geom_errorbar(aes(ymin=mean_MI-sd, ymax=mean_MI+sd), alpha=0.5, width=0) +
    ylab("Estimated I(X;Y|Z1,Z2)") +
    scale_color_discrete(drop=TRUE, limits = levels(results$method)) +
    theme_bw() + 
    facet_grid(N ~ .) +
    theme(legend.position = "top")
plot(estimation)
```

![](FigS6_files/figure-gfm/plot-1.png)<!-- -->

``` r
error = results %>% filter(method %in% c("miic", "knn3", "knn10", "knn20")) %>%
                    group_by(method,rho, N) %>%
                    mutate(error = MI-get_true_CMI(rho)) %>% summarize(MSE = mean(error**2)) %>%
    ggplot(aes(x=rho, y=MSE, color=method)) +
    geom_line() +
    scale_color_discrete(drop=TRUE, limits = levels(results$method)) +
    scale_y_continuous(trans='log10') +
    theme_bw() + 
    facet_grid(N ~ .) +
    theme(legend.position = "top")
plot(error)
```

![](FigS6_files/figure-gfm/plot-2.png)<!-- -->

``` r
variance = results %>% filter(method %in% c("miic", "knn3", "knn10", "knn20")) %>%
                       group_by(method,rho,N) %>%
                       summarize(variance = sd(MI)) %>%
    ggplot(aes(x=rho, y=variance, color=method)) +
    geom_line() +
    ylab("Estimator variance") +
    scale_color_discrete(drop=TRUE, limits = levels(results$method)) +
    #scale_y_continuous(trans='pseudo_log') +
    theme_bw() + 
    facet_grid(N ~ .) +
    theme(legend.position = "top")
plot(variance)
```

![](FigS6_files/figure-gfm/plot-3.png)<!-- -->

While conditional mutual information estimation is much harder than
pairwise MI estimation, we can see once again that `miic` gives decent
estimation for all ![\\rho](https://latex.codecogs.com/png.latex?%5Crho
"\\rho") and especially at conditional independence.

# References

<div id="refs" class="references">

<div id="ref-lizier_jidt_2014">

Lizier, Joseph T. 2014. “JIDT: An Information-Theoretic Toolkit for
Studying the Dynamics of Complex Systems.” *Frontiers in Robotics and
AI* 1: 11.

</div>

</div>
