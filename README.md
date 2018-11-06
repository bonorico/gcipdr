# gcipdr
Gaussian Copula (based) Individual Person Data Reproduction

# Installation instructions in R (follow in order):

## step 1: open your R and run command below
```
library(devtools)
install_github("bonorico/gcipdr")
```
# step 2: Important !

gcipdr depends on CRAN archived package 'JohnsonDistribution'. You must install this dependency from source by doing the following:
```
 url <- "https://cran.r-project.org/src/contrib/Archive/JohnsonDistribution/JohnsonDistribution_0.24.tar.gz"
 pkgFile <- "JohnsonDistribution_0.24.tar.gz"
 download.file(url = url, destfile = pkgFile)

 install.packages(pkgs=pkgFile, type="source", repos=NULL)

 unlink(pkgFile)
```

# ~+~+~+~+~+~+~+~  ~+~+~+~+~+~+~+~  ~+~+~+~+~+~+~+~  ~+~+~+~+~+~+~+~ ~+~+~+~+~+~+~+~  ~+~+~+~+~+~+~+~ ~+~+~+~+~+~+~+~   

# Package Description (from DESCRIPTION file)

The proposed Gaussian Copula technique reconstructs (unavailable) original IPD from simple IPD summaries only (IPD empirical marginal moments and correlation matrix). Ultimate focus is recovery of original IPD inferences from the reconstructed IPD. This methods can be useful there where original raw data (IPD), and inferences associated to it, are not available but primarily sought: it is assumed only simple distributional IPD summaries can be recovered as analysis' input data only. Applications are envisioned in fields like Statistical Disclosure Control, Distributed Computing (within infrastructures like DataShield -- see https://cran.datashield.org/), and, provided the required above IPD summaries can be retrieved, Research Synthesis (like meta-analysis) and Reproduction.

# Further Details

This package is a beta version and extensive editing, documenting, and
amelioration shall be expected in the future.


# Installation Notes:

Executing command 'install_github()' above should cause to fresh (re)install LinkingTo packages listed in the DESCRIPTION file (Rcpp, RcppArmadillo). In my case the command also caused to fresh (re)install other dependencies ('cubature' and 'ggplot2') but not all of them. The full list of dependencies ('moment', 'parallel', 'cubature', 'mvtnorm', 'ggplot2') is issued in a warning at installation time. Install any missing dependency manually. 

If you want to avoid automatically fresh re-installing any of the above dependencies, use the following ALTERNATIVE INSTALLATION PLAN (Linux commands -- should work similarly with other OS):

1) Download all content of the gcipdr repository on your local machine (in a folder called 'gcipdr')

2) From Terminal: cd into 'your_local_directory/gcipdr' and run: 
                                         
                                                                  R CMD build gcipdr
                                                                   
3) From R:

                                     pkg <- "your_local_directory/gcipdr_0.0.tar.gz"
                            install.packages(pkgs=pkg, type="source", repos=NULL)
                          
4) From R: all remnant dependencies must be installed manually, if yet not installed on your machine.



Installation works fine on my local machine: Linux 3.13.0-147-generic #196-Ubuntu x86_64 (with 4 Intel CPUs) GNU/Linux. 
At installation you could receive the following 

WARNING: support for OpenMP requires C++11/C++14; add -std=c++11 or -std=c++14 to compiler flags

You should be able to just ignore it.

# Help files

After installation type '?DataRebuild' or '?Simulate.many.datasets' in the R command line to obtain more information on gcipdr modules.
