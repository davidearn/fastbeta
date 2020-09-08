# fastbeta
Fast estimation of time-varying infectious disease transmission rates

## Installation

### Local

To install from a local repository:

* Check that you are in the `devel` branch with `git branch`.
* Check that you are in the root directory of the package with `ls`.
* Run `make`.

### Remote

To install from the remote repository, run this script in R:

```r
if (!require(remotes)) install.packages("remotes")
remotes::install_github("davidearn/fastbeta", ref = "devel", build_vignettes = TRUE)
```

## Documentation

Package documentation can be browsed like so:

```r
library(fastbeta)

## List of exported functions
ls("package:fastbeta")

## Vignette
vignette("fastbeta-vignette")

## Help pages
#options(help_type = "html") # if running R from command line
?"fastbeta-package"   # package
?function_name        # exported function "function_name"
?"class_name-methods" # S3 methods for class "class_name"
```
