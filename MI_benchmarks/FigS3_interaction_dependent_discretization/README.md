FigS3 Interaction-dependent optimum discretization.
================

Optimum bivariate partitions obtained from
![N=1,000](https://latex.codecogs.com/png.latex?N%3D1%2C000 "N=1,000")
samples of two different joint distributions
![P(X,Y)](https://latex.codecogs.com/png.latex?P%28X%2CY%29 "P(X,Y)")
sharing the same sampling of ![X](https://latex.codecogs.com/png.latex?X
"X") taken from a uniform distribution on
![\[0, 0.3\]](https://latex.codecogs.com/png.latex?%5B0%2C%200.3%5D
"[0, 0.3]"), but with different dependences for
![Y](https://latex.codecogs.com/png.latex?Y "Y"). {}
![Y](https://latex.codecogs.com/png.latex?Y "Y") is defined as
![\\log(X) +
\\epsilon\_1](https://latex.codecogs.com/png.latex?%5Clog%28X%29%20%2B%20%5Cepsilon_1
"\\log(X) + \\epsilon_1"), and {}
![Y](https://latex.codecogs.com/png.latex?Y "Y") is defined as ![X^5 +
\\epsilon\_2](https://latex.codecogs.com/png.latex?X%5E5%20%2B%20%5Cepsilon_2
"X^5 + \\epsilon_2"), where
![\\epsilon\_1](https://latex.codecogs.com/png.latex?%5Cepsilon_1
"\\epsilon_1") and
![\\epsilon\_2](https://latex.codecogs.com/png.latex?%5Cepsilon_2
"\\epsilon_2") are Gaussian noise terms chosen so that the mutual
informations of both examples are comparable,
{![I(X;Y)\\simeq 0.75](https://latex.codecogs.com/png.latex?I%28X%3BY%29%5Csimeq%200.75
"I(X;Y)\\simeq 0.75")}.

This example shows that the optimum partition for
![X](https://latex.codecogs.com/png.latex?X "X") depends on its specific
relation with ![Y](https://latex.codecogs.com/png.latex?Y "Y") and needs
to be discretized with finer partitions in {} at low
![X](https://latex.codecogs.com/png.latex?X "X") values for which
![Y\\simeq \\log
X](https://latex.codecogs.com/png.latex?Y%5Csimeq%20%5Clog%20X
"Y\\simeq \\log X") varies the most and in {} at higher
![X](https://latex.codecogs.com/png.latex?X "X") values for ![Y\\simeq
X^5](https://latex.codecogs.com/png.latex?Y%5Csimeq%20X%5E5
"Y\\simeq X^5").

``` r
library(miic)

N = 1000

X = runif(N, min = 0, max=0.3)
noise = c(rnorm(N, sd=0.4))
Y1 = log(X) + noise
miic::discretizeMutual(X, Y1)$infok
```

![](FigS3_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

    ## [1] 0.63026

``` r
Y2 = X**5 + noise/1900
miic::discretizeMutual(X, Y2)$infok
```

![](FigS3_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

    ## [1] 0.66306
