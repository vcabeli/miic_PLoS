# Learning clinical networks from medical records based on information estimates in mixed-type data

This is the repository containing the code and figures for the 2020 paper published in PLoS 
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007866

![NML-optimal 2d discretization](https://journals.plos.org/ploscompbiol/article/figure/image?size=large&id=info:doi/10.1371/journal.pcbi.1007866.g002)

## Installation

You need the package `miic` and its dependencies to run the code.
Use the package installer directly in R for the latest stable CRAN version
```R
install.packages("miic", dep = T)
```

or install the package from the tarball (https://github.com/miicTeam/miic_R_package) using the console command line:
```bash
git clone https://github.com/miicTeam/miic_R_package.git
cd miic_R_package
R CMD INSTALL .
```


## License
[MIT](https://choosealicense.com/licenses/mit/)
