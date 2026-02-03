## Updated package (formerly archived)
This is a patch to update my email address. All other changes are just to get checks/tests running smoothly.

## Test environments
* local macOS install, R 4.2.2
- macOS builder
- win-builder (devel and release)
- github actions:
    - {os: windows-latest, r: 'release'}
    - {os: ubuntu-latest,   r: 'devel', http-user-agent: 'release'}
    - {os: ubuntu-latest,   r: 'release'}
    - {os: ubuntu-latest,   r: 'oldrel-1'}

## R CMD check results
There were no ERRORs or WARNINGs.

There were the following NOTEs:

Maintainer: 'Katherine Brumberg <kbrum@umich.edu>'

New submission

Package was archived on CRAN

Possibly misspelled words in DESCRIPTION:
  Brumberg (14:66, 14:118)
  al (14:78, 14:130)
  et (14:75, 14:127)

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2025-06-12 as email to the maintainer is
    undeliverable.

Suggests or Enhances not in mainstream repositories:
  gurobi
  
  'gurobi' is a commercial optimization software, and is not required for the use of my package, although it is recommended. Instructions for installation are included in the description.


## Downstream dependencies
There are currently no downstream dependencies for this package.
