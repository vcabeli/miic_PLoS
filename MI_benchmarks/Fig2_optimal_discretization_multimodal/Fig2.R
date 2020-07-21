library(MASS)
library(copula)
library(distr)
library(miic)

# The proposed information-maximizing discretization scheme is illustrated for a joint distribution
# defined as a Gumbel bivariate copula with parameter θ=5 and marginal distributions chosen as 
# Gaussian mixtures with three equiprobable peaks and respective means and variances, 
# μX={0,4,6}, σX={1,2,0.7} and μY={−3,6,9},σY={2,0.5,0.5}.

mixnorm1 <- UnivarMixingDistribution(Norm(0,1), Norm(4,2), Norm(6,0.7))
rmixnorm1 <- r(mixnorm1)
dmixnorm1 <- d(mixnorm1)
qmixnorm1 <- q(mixnorm1)
pmixnorm1 <- p(mixnorm1)


mixnorm2 <- UnivarMixingDistribution(Norm(-3,2), Norm(6,0.5), Norm(9,0.5))
rmixnorm2 <- r(mixnorm2)
dmixnorm2 <- d(mixnorm2)
qmixnorm2 <- q(mixnorm2)
pmixnorm2 <- p(mixnorm2)


mv.NE <- mvdc(gumbelCopula(5), margins=c("mixnorm1","mixnorm2"), paramMargins=list(list(),list()))


for(N in c(500, 1000, 5000, 10000)){
    data = rMvdc(N, mv.NE)
    print(discretizeMutual(data[,1], data[,2])$info)
}
