FigS2 Adaptive information-maximizing partitions depending on
interaction strength
================

To assess the range in bin numbers depending on the strength of
interaction between variables, we generated N = 1,000 independent
samples for 10,000 Gaussian bivariate distributions with a uniformly
distributed correlation coefficient ρ in \[−1, 1\]. The real mutual
information (RI) of Gaussian bivariate distributions can be computed
directly as ![RI(X; Y ) = - log(1 - \\rho^2
)/2](https://latex.codecogs.com/png.latex?RI%28X%3B%20Y%20%29%20%3D%20-%20log%281%20-%20%5Crho%5E2%20%29%2F2
"RI(X; Y ) = - log(1 - \\rho^2 )/2") (Cover and Thomas 2012). For each
pair (X, Y ), we estimated the mutual information with the proposed
optimum bivariate discretization as well as the Maximal Information
Coefficient (Reshef et al. 2011) using the `minepy` package. We counted
the number of bins on each dimension for the optimal discretization
found by both methods.

To run this script you will need a working python installation with the
[minepy package](https://minepy.readthedocs.io/en/latest/) (Albanese et
al. 2013).

``` r
library(reticulate)

# Import minepy from python
minepy <- import("minepy")
mine = minepy$MINE(alpha=0.6, c=15, est="mic_approx")
```

Thanks to `reticulate`, We can call both the `miic` R package and the
python `minepy` on the same data.

``` r
library(MASS)
library(miic)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:MASS':
    ## 
    ##     select

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
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
```

    ## [1] "10%"
    ## [1] "20%"
    ## [1] "30%"
    ## [1] "40%"
    ## [1] "50%"
    ## [1] "60%"
    ## [1] "70%"
    ## [1] "80%"
    ## [1] "90%"
    ## [1] "100%"

We can now plot the number of bins found on X and Y by the optimal
discretization on our normal bivariate distributions. Since the
marginals are the same and the joint distribution is symmetrical, we
expect roughly the same number of bins on X and Y. For the `miic`
discretization, we also expect more bins for stronger correlations and
less for weaker correlations.

### miic results

``` r
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
```

![](FigS2_files/figure-gfm/plot%20miic-1.png)<!-- -->

The information-maximizing partition proposed in the present paper
behaves as expected: the number of bins on each variable is roughly
similar and scales monotonically with the strength of the interaction
between variables. This implies that additional bins are only introduced
when their associated complexity cost is justified by a larger gain in
mutual information. Conversely, when the information between X and Y
approaches zero, both variables are partitioned into fewer and fewer
bins until a single bin is selected for each variable, when they are
inferred to be independent, given the available data.

### MINE results

``` r
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
```

![](FigS2_files/figure-gfm/plot%20mine-1.png)<!-- -->

The partition chosen to estimate the Maximal Information Coefficient
(Reshef et al. 2011) is very different, regardless of the interaction
strength, as it systematically corresponds to an unbalanced distribution
of bins between the two variables, with one variable usually partitioned
into the maximum number of bins (set by default to
![\\text{floor}(N^{0.6}/2)
= 31](https://latex.codecogs.com/png.latex?%5Ctext%7Bfloor%7D%28N%5E%7B0.6%7D%2F2%29%20%3D%2031
"\\text{floor}(N^{0.6}/2) = 31")) while the other is discretized into
two levels only. This result is not unexpected, however, as the Maximal
Information Coefficient is defined by maximizing the mutual information
of the discretized variables over the grid,
![I(\[X\]\_{\\Delta\_x};\[Y\]\_{\\Delta\_y})](https://latex.codecogs.com/png.latex?I%28%5BX%5D_%7B%5CDelta_x%7D%3B%5BY%5D_%7B%5CDelta_y%7D%29
"I([X]_{\\Delta_x};[Y]_{\\Delta_y})"), normalized by the minimum of
![\\text{log}
\\Delta\_x](https://latex.codecogs.com/png.latex?%5Ctext%7Blog%7D%20%5CDelta_x
"\\text{log} \\Delta_x") and ![\\text{log}
\\Delta\_y](https://latex.codecogs.com/png.latex?%5Ctext%7Blog%7D%20%5CDelta_y
"\\text{log} \\Delta_y"). Indeed, maximizing the normalized mutual
information is done by partitioning as few samples as possible into the
maximum number of bins in one dimension (as sketched in Fig. 1), while
simultaneously minimizing the number of bins, and thus
![\\text{log}\\Delta\_i](https://latex.codecogs.com/png.latex?%5Ctext%7Blog%7D%5CDelta_i
"\\text{log}\\Delta_i") , in the other dimension. See further discussion
in (Kinney and Atwal 2014).

# References

<div id="refs" class="references">

<div id="ref-albanese_minerva_2013">

Albanese, Davide, Michele Filosi, Roberto Visintainer, Samantha
Riccadonna, Giuseppe Jurman, and Cesare Furlanello. 2013. “Minerva and
Minepy: A c Engine for the MINE Suite and Its R, Python and MATLAB
Wrappers.” *Bioinformatics (Oxford, England)* 29 (3): 407–8.
<https://doi.org/10.1093/bioinformatics/bts707>.

</div>

<div id="ref-cover_elements_2012">

Cover, Thomas M., and Joy A. Thomas. 2012. *Elements of Information
Theory*. John Wiley & Sons.

</div>

<div id="ref-kinney_equitability_2014">

Kinney, Justin B., and Gurinder S. Atwal. 2014. “Equitability, Mutual
Information, and the Maximal Information Coefficient.” *Proceedings of
the National Academy of Sciences* 111 (9): 3354–9.
<https://doi.org/10.1073/pnas.1309933111>.

</div>

<div id="ref-reshef_detecting_2011">

Reshef, David N., Yakir A. Reshef, Hilary K. Finucane, Sharon R.
Grossman, Gilean McVean, Peter J. Turnbaugh, Eric S. Lander, Michael
Mitzenmacher, and Pardis C. Sabeti. 2011. “Detecting Novel Associations
in Large Datasets.” *Science (New York, N.y.)* 334 (6062): 1518–24.
<https://doi.org/10.1126/science.1205438>.

</div>

</div>
