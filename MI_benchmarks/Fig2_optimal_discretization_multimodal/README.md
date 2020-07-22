Fig2 Optimum bivariate discretization for mutual information estimation
================

The proposed information-maximizing discretization scheme is illustrated
for a joint distribution defined as a Gumbel bivariate copula with
parameter ![\\theta
= 5](https://latex.codecogs.com/png.latex?%5Ctheta%20%3D%205
"\\theta = 5") and marginal distributions chosen as Gaussian mixtures
with three equiprobable peaks and respective means and variances,
![\\mu\_X =
\\{0,4,6\\}](https://latex.codecogs.com/png.latex?%5Cmu_X%20%3D%20%5C%7B0%2C4%2C6%5C%7D
"\\mu_X = \\{0,4,6\\}"), ![\\sigma\_X =
\\{1,2,0.7\\}](https://latex.codecogs.com/png.latex?%5Csigma_X%20%3D%20%5C%7B1%2C2%2C0.7%5C%7D
"\\sigma_X = \\{1,2,0.7\\}") and ![\\mu\_Y =
\\{-3,6,9\\}](https://latex.codecogs.com/png.latex?%5Cmu_Y%20%3D%20%5C%7B-3%2C6%2C9%5C%7D
"\\mu_Y = \\{-3,6,9\\}") , ![\\sigma\_Y =
\\{2,0.5,0.5\\}](https://latex.codecogs.com/png.latex?%5Csigma_Y%20%3D%20%5C%7B2%2C0.5%2C0.5%5C%7D
"\\sigma_Y = \\{2,0.5,0.5\\}").

``` r
library(MASS)
library(copula)
library(distr)
library(miic)
```

Defining marginals as multimodal distributions :

``` r
mixnorm1 <- UnivarMixingDistribution(Norm(0,1), Norm(4,2), Norm(6,0.7))
dmixnorm1 <- d(mixnorm1); qmixnorm1 <- q(mixnorm1); pmixnorm1 <- p(mixnorm1)


mixnorm2 <- UnivarMixingDistribution(Norm(-3,2), Norm(6,0.5), Norm(9,0.5))
dmixnorm2 <- d(mixnorm2); qmixnorm2 <- q(mixnorm2); pmixnorm2 <- p(mixnorm2)
```

Defining the 2d distribution with a [Gumbel
Copula](https://en.wikipedia.org/wiki/Copula_\(probability_theory\)#Archimedean_copulas)
:

``` r
mv.NE <- mvdc(gumbelCopula(5), margins=c("mixnorm1","mixnorm2"), paramMargins=list(list(),list()))
```

Generating samples and running the dynamic discretization :

``` r
for(N in c(500,1000,10000)){
  data = rMvdc(N, mv.NE)
  print(discretizeMutual(data[,1], data[,2])$info)
}
```

![](Fig2_files/figure-gfm/plot-1.png)<!-- -->

    ## [1] 0.98682

![](Fig2_files/figure-gfm/plot-2.png)<!-- -->

    ## [1] 1.10994

![](Fig2_files/figure-gfm/plot-3.png)<!-- -->

    ## [1] 1.14539
