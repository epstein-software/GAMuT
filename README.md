![R](Rlogo.png)
## GAMuT:  Gene Association with Multiple Traits
### Introduction
In [Broadaway et al. (2016)](http://www.cell.com/ajhg/abstract/S0002-9297(16)00052-5),
we created a test based on kernel distance-covariance methodology called GAMuT for gene-based association testing of rare variants with multiple phenotypes in a sample of unrelated subjects. The analyzed phenotypes can be either continuous or categorical in nature.  The method also allows for adjustment of covariates using adjusted residuals. 

GAMuT is implemented as a set of R routines, which can be executed within an R environment. 

### The R environment
R is a widely-used, free and open source software environment for statistical computing and graphics. The most recent version of R can be downloaded from the 
[Comprehensive R Archive Network (CRAN)](http://cran.r-project.org/)
CRAN provides precompiled binary versions of R for Windows, MacOS, and select Linux distributions that are likely sufficient for many users' needs.  Users can also install R from source code;  however this may require a significant amount of effort.  For specific details on how to compile, install, and manage R and R-packages, refer to the manual [R Installation and Administration](http://cran.r-project.org/doc/manuals/r-release/R-admin.html).


### R packages required for analysis
GAMuT requires the installation of the following R library:

[CompQuadForm](https://cran.r-project.org/web/packages/CompQuadForm/index.html)

The easiest method to install these packages is with the following command entered in an R shell:

    install.packages("CompQuadForm")

One can also [install R packages from the command line]
(http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).


### Running GAMuT
For the example, we provide sample files consisting of 1000 unrelated variants
possessing 6 phenotypes and who are sequenced for a set of rare variants in a gene of interest.
We show in our example R code how to implement the GAMuT test to perform an association test
between the 6 phenotypes and the rare variants within the gene of interest. 

In the [GAMuT example analysis page](http://genetics.emory.edu/labs/epstein/software/gamut/GAMuT-example-analysis.html),
we provide a step-by-step instruction on how the GAMuT code operates here.


### Questions and technical support
For questions or concerns with the GAMuT R routines, please contact
[Richard Duncan](mailto:rduncan@emory.edu) and 
[Michael Epstein](mailto:mpepste@emory.edu)

We appreciate any feedback you have with our site and instructions.
