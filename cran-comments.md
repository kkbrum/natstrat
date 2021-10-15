## Updated package
This is an update to an existing package. In this major version I have added 
additional functionality, which is almost entirely backwards compatible.

## Test environments
* local windows install, R 4.0.3
* local windows install, R 4.1.1
* ubuntu-20.04 release (on github actions)
* ubuntu-20.04 devel (on github actions)
* macOS-latest release (on github actions)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs on release or devel.

Only in my local windows install, using version 4.0.3 of R,
there is 1 WARNING:
     Warning: package 'ggplot2' was built under R version 4.0.5
As my package does not rely heavily on ggplot2, there should not be any issues
with users using my package with earlier versions of R. There do not seem
to have been any ggplot2 changes that affect my package.

There were the following NOTEs:

* Suggests or Enhances not in mainstream repositories:
  gurobi
* Package suggested but not available for checking: 'gurobi'
  
  'gurobi' is a commercial optimization software, and is not required for the use of my package, although it is recommended. Instructions for installation are included in the description.
  







## Downstream dependencies
There are currently no downstream dependencies for this package.
