## Test environments
* local windows install, R 4.0.3
* ubuntu-20.04 release (on github actions)
* ubuntu-20.04 devel (on github actions)
* macOS-latest release (on github actions)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Katherine Brumberg <kbrum@wharton.upenn.edu>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Gurobi (12:38)
  Unweighted (3:15)
  gurobi (13:43)
  
  "Unweighted" is how we describe the resulting strata of our method.
  
  "Gurobi" is the name of a commerical optimization software, with related R package "gurobi."
  
* checking package dependencies ... NOTE
Package suggested but not available for checking: 'gurobi'

As "gurobi" is a commercial software, the R package needs to be installed through a download from their website and is not available on CRAN. This package is not required for my package as I provide an alternative using "Rglpk" from CRAN, but it is recommended. I have checked for the gurobi package in the code before using it with requireNamespace(), and have mentioned gurobi and its installation in the description of the package.

## Downstream dependencies
There are currently no downstream dependencies for this package.
