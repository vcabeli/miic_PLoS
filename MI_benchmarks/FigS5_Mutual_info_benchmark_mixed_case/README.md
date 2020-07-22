FigS5 Mutual information estimation of mixed variable
================

Here we propose benchmarks to evaluate our MI estimation for mixed cases
with experiments taken from (Gao et al. 2017) for increasing sample size
![N](https://latex.codecogs.com/png.latex?N "N").

We will use `reticulate` to call python functions from R to use Gao et
al’s method for estimating mutual information between mixed variables
(Gao et al. 2017) (“mixed KSG Gao”) which you can download here :
<https://github.com/wgao9/mixed_KSG/blob/master/mixed.py>

``` r
library(reticulate)
################################################################################
# Mixed KSG from Gao et al. 2017
source_python("~/work/projects/miic_dyn_programming_MI/mixed_KSG/mixed.py")
################################################################################
```

We will use `rJava` to use the JIDT toolbox (Lizier 2014)
(<https://github.com/jlizier/jidt>) for other KSG estimators, namely the
original Kraskov (Kraskov, Stögbauer, and Grassberger 2004) with noise
(“noisy KSG”), Ross’s approach also based on KSG (Ross 2014) (“mixed
KSG Ross”), as well as a method based of pdf estimation with kernels
(“kernel”). In all KNN approaches, we keep the parameter
![k=5](https://latex.codecogs.com/png.latex?k%3D5 "k=5") by
default.

``` r
################################################################################
# KNN estimation methods with JIDT
library("rJava")
.jinit()
.jaddClassPath("/home/vcabeli/Downloads/jidt/infodynamics.jar")

jidt_compute_mi_ksg = function(X,Y, k=5, noise=T){
  miCalc<-.jnew("infodynamics/measures/continuous/kraskov/MutualInfoCalculatorMultiVariateKraskov1")
  .jcall(miCalc,"V","setProperty", "k", as.character(k))
  .jcall(miCalc,"V","setProperty", "addNoise", ifelse(noise, "true", "false"))
  .jcall(miCalc,"V","initialise")
  .jcall(miCalc, "V", "setObservations", X, Y)
  
  return(.jcall(miCalc,"D","computeAverageLocalOfObservations"))
}

jidt_compute_mi_ksg_mixed = function(X,Y, k=5, noise=T){
  miCalc<-.jnew("infodynamics/measures/mixed/kraskov/MutualInfoCalculatorMultiVariateWithDiscreteKraskov")
  .jcall(miCalc,"V","setProperty", "k", as.character(k))
  .jcall(miCalc,"V","setProperty", "addNoise", ifelse(noise, "true", "false"))
  .jcall(miCalc,"V","initialise", as.integer(1), length(unique(Y)))
  .jcall(miCalc,"V","setObservations",
         .jarray(lapply(X, .jarray),"[D", dispatch = TRUE),
         .jarray(Y, "I", dispatch=TRUE))
  
  return(.jcall(miCalc,"D","computeAverageLocalOfObservations"))
}

jidt_compute_mi_kernel = function(X,Y){
  miCalc<-.jnew("infodynamics/measures/continuous/kernel/MutualInfoCalculatorMultiVariateKernel")
  .jcall(miCalc,"V","initialise", as.integer(1), as.integer(1))
  .jcall(miCalc, "V", "setObservations", X, Y)
  
  return(.jcall(miCalc,"D","computeAverageLocalOfObservations"))
}
################################################################################
```

Defining a plot function and the number of runs:

``` r
plot_res = function(res, ground_truth){
  plot = res %>% group_by(N,method) %>%
    summarise(MSE = mean((value - ground_truth)**2), sd = sd((value - ground_truth)**2)) %>%
    ggplot(aes(x=N, y=MSE, color=method)) +
    #geom_errorbar(aes(ymin=MSE-sd, ymax=MSE+sd)) +
    geom_line() +
    geom_point() +
    scale_y_log10(lim=c(-6,1)) +
    theme_bw()+
    theme(legend.position = "top")
  
  return(plot)
}

nreps = 20
```

## Model 1

Model 1 has ![X](https://latex.codecogs.com/png.latex?X "X") and
![Y](https://latex.codecogs.com/png.latex?Y "Y") partly continuous and
partly discrete.

``` r
###
# Model 1
res_1 = data.frame(N=numeric(), method=character(), value=numeric(), stringsAsFactors = F)
rho = 0.9
prop = 0.5
sigma = matrix(c(1, rho, rho, 1), 2)

ground_truth = -prop*0.5*log(det(sigma)) + prop*(log(2)+rho*log(rho)+(1-rho)*log(1-rho))-prop*log(prop)-prop*log(prop)

for(N in seq(100,3100, by=500)){
  print(N)
  for(rep in 1:nreps){
    cont_part = mvrnorm(N, mu = c(0,0), Sigma = sigma)
    Xcont = cont_part[,1]
    Ycont = cont_part[,2]
    Xdisc = rbinom(N,1,0.5)
    Ydisc = (Xdisc + rbinom(N,1,1-rho))%%2
    Xdisc = 2*Xdisc-matrix(1,1,N)
    Ydisc = 2*Ydisc-matrix(1,1,N)
    samples = sample(1:N*2, N, replace = F)
    X = c(Xcont,Xdisc)[samples]
    Y = c(Ycont,Ydisc)[samples]
    
    res_1[nrow(res_1)+1,] = list(N, "miic", miic::discretizeMutual(X,Y,maxbins = N, plot = F)$info)
    res_1[nrow(res_1)+1,] = list(N, "mixed KSG Gao", Mixed_KSG(matrix(X), matrix(Y)))
    res_1[nrow(res_1)+1,] = list(N, "noisy KSG", jidt_compute_mi_ksg(X,Y))
    res_1[nrow(res_1)+1,] = list(N, "kernel", jidt_compute_mi_kernel(X,Y))
    res_1[nrow(res_1)+1,] = list(N, "equal freq (N^1/3)",  mutinformation(discretize(X, disc="equalwidth", nbins = N**(1/3)), 
                                                                      discretize(Y, disc="equalwidth", nbins = N**(1/3))))
  }
}
```

    ## [1] 100
    ## [1] 600
    ## [1] 1100
    ## [1] 1600
    ## [1] 2100
    ## [1] 2600
    ## [1] 3100

``` r
example = miic::discretizeMutual(X,Y)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

![](FigS5_files/figure-gfm/model%201-1.png)<!-- -->

``` r
plot_res(res_1, ground_truth)
```

    ## Warning in trans$transform(limits): NaNs produced

![](FigS5_files/figure-gfm/model%201-2.png)<!-- -->

## Model 2

Model 2 has a truly continuous
![Y](https://latex.codecogs.com/png.latex?Y "Y") made of a mixture of
normal distributions which means are based on the discrete values of
![X](https://latex.codecogs.com/png.latex?X "X").

``` r
###
# Model 2
res_2 = data.frame(N=numeric(), method=character(), value=numeric(), stringsAsFactors = F)
nlevels = 5
ground_truth = log(nlevels) - (nlevels-1)*log(2)/nlevels

for(N in seq(100,3100, by=500)){
  print(N)
  for(rep in 1:nreps){
    X = sample(1:nlevels, N, replace = T)-1
    Y = sapply(X, function(x) runif(1, x, x+2))
    res_2[nrow(res_2)+1,] = list(N, "miic", miic::discretizeMutual(X,Y,is_discrete = c(T,F), maxbins = N, plot = F)$info)
    res_2[nrow(res_2)+1,] = list(N, "mixed KSG Gao", Mixed_KSG(matrix(X), matrix(Y)))
    res_2[nrow(res_2)+1,] = list(N, "mixed KSG Ross", jidt_compute_mi_ksg_mixed(Y,as.integer(X)))
    res_2[nrow(res_2)+1,] = list(N, "noisy KSG", jidt_compute_mi_ksg(as.double(X),Y))
    res_2[nrow(res_2)+1,] = list(N, "kernel", jidt_compute_mi_kernel(as.double(X),Y))
    res_2[nrow(res_2)+1,] = list(N, "equal freq (N^1/3)",  mutinformation(X,
                                                                          discretize(Y, disc="equalfreq", nbins = N**(1/3))))
  }
}
```

    ## [1] 100
    ## [1] 600
    ## [1] 1100
    ## [1] 1600
    ## [1] 2100
    ## [1] 2600
    ## [1] 3100

``` r
example = miic::discretizeMutual(X,Y)
```

![](FigS5_files/figure-gfm/model%202-1.png)<!-- -->

``` r
plot_res(res_2, ground_truth)
```

    ## Warning in trans$transform(limits): NaNs produced

![](FigS5_files/figure-gfm/model%202-2.png)<!-- -->

## Model 4

Model 4 has a truly continuous
![X](https://latex.codecogs.com/png.latex?X "X") and a discrete
![Y](https://latex.codecogs.com/png.latex?Y "Y") based on a poisson law
with ![X](https://latex.codecogs.com/png.latex?X "X").

``` r
###
# Model 4
res_4 = data.frame(N=numeric(), method=character(), value=numeric(), stringsAsFactors = F)
p = 0

ground_truth = (1-p)*0.3012

for(N in seq(100,3100, by=500)){
  print(N)
  for(rep in 1:nreps){
    min_bin_size = 0
    while(min_bin_size <= 5){
      X=rexp(N,1)
      Z=rpois(N,X)
      Y0=rbinom(N,1,1-p)
      Y=Y0*Z
      Y = factor(Y)
      levels(Y) = 0:nlevels(Y)
      Y = as.numeric(Y)-1
      
      min_bin_size = min(table(Y))
    }
    
    res_4[nrow(res_4)+1,] = list(N, "miic", miic::discretizeMutual(X,Y,is_discrete = c(F,T), maxbins = N, plot = F)$info)
    res_4[nrow(res_4)+1,] = list(N, "mixed KSG Gao", Mixed_KSG(matrix(X), matrix(Y)))
    res_4[nrow(res_4)+1,] = list(N, "mixed KSG Ross", jidt_compute_mi_ksg_mixed(X,as.integer(Y)))
    res_4[nrow(res_4)+1,] = list(N, "noisy KSG", jidt_compute_mi_ksg(X,as.double(Y)))
    res_4[nrow(res_4)+1,] = list(N, "kernel", jidt_compute_mi_kernel(X,as.double(Y)))
    res_4[nrow(res_4)+1,] = list(N, "equal freq (N^1/3)", mutinformation(discretize(X, disc="equalwidth", nbins = N**(1/3)),
                                                                         Y))
  }
}
```

    ## [1] 100
    ## [1] 600
    ## [1] 1100
    ## [1] 1600
    ## [1] 2100
    ## [1] 2600
    ## [1] 3100

``` r
example = miic::discretizeMutual(X,Y)
```

![](FigS5_files/figure-gfm/model%204-1.png)<!-- -->

``` r
plot_res(res_4, ground_truth)
```

    ## Warning in trans$transform(limits): NaNs produced

![](FigS5_files/figure-gfm/model%204-2.png)<!-- -->

## Model 5

Model 5 is the same as model 4 with slightly less information.

``` r
###
# Model 5
res_5 = data.frame(N=numeric(), method=character(), value=numeric(), stringsAsFactors = F)
p = 0.15

ground_truth = (1-p)*0.3012

for(N in seq(100,3100, by=500)){
  print(N)
  for(rep in 1:nreps){
    min_bin_size = 0
    while(min_bin_size <= 5){
      X=rexp(N,1)
      Z=rpois(N,X)
      Y0=rbinom(N,1,1-p)
      Y=Y0*Z
      Y = factor(Y)
      levels(Y) = 0:nlevels(Y)
      Y = as.numeric(Y)-1
      
      min_bin_size = min(table(Y))
    }
    
    res_5[nrow(res_5)+1,] = list(N, "miic", miic::discretizeMutual(X,Y,is_discrete = c(F,T), maxbins = N, plot = F)$info)
    res_5[nrow(res_5)+1,] = list(N, "mixed KSG Gao", Mixed_KSG(matrix(X), matrix(Y)))
    res_5[nrow(res_5)+1,] = list(N, "mixed KSG Ross", jidt_compute_mi_ksg_mixed(X,as.integer(Y)))
    res_5[nrow(res_5)+1,] = list(N, "noisy KSG", jidt_compute_mi_ksg(X,as.double(Y)))
    res_5[nrow(res_5)+1,] = list(N, "kernel", jidt_compute_mi_kernel(X,as.double(Y)))
    res_5[nrow(res_5)+1,] = list(N, "equal freq (N^1/3)",  mutinformation(discretize(X, disc="equalwidth", nbins = N**(1/3)),
                                                                          Y))
  }
}
```

    ## [1] 100
    ## [1] 600
    ## [1] 1100
    ## [1] 1600
    ## [1] 2100
    ## [1] 2600
    ## [1] 3100

``` r
example = miic::discretizeMutual(X,Y)
```

![](FigS5_files/figure-gfm/model%205-1.png)<!-- -->

``` r
plot_res(res_5, ground_truth)
```

    ## Warning in trans$transform(limits): NaNs produced

![](FigS5_files/figure-gfm/model%205-2.png)<!-- -->

# References

<div id="refs" class="references">

<div id="ref-gao_estimating_2017">

Gao, Weihao, Sreeram Kannan, Sewoong Oh, and Pramod Viswanath. 2017.
“Estimating Mutual Information for Discrete-Continuous Mixtures.”
*arXiv:1709.06212 \[Cs, Math\]*, September.
<http://arxiv.org/abs/1709.06212>.

</div>

<div id="ref-kraskov_estimating_2004">

Kraskov, Alexander, Harald Stögbauer, and Peter Grassberger. 2004.
“Estimating Mutual Information.” *Physical Review E* 69 (6): 066138.

</div>

<div id="ref-lizier_jidt_2014">

Lizier, Joseph T. 2014. “JIDT: An Information-Theoretic Toolkit for
Studying the Dynamics of Complex Systems.” *Frontiers in Robotics and
AI* 1: 11.

</div>

<div id="ref-ross_mutual_2014">

Ross, Brian C. 2014. “Mutual Information Between Discrete and Continuous
Data Sets.” *PLOS ONE* 9 (2): e87357.
<https://doi.org/10.1371/journal.pone.0087357>.

</div>

</div>
