SepTest
================

## Overview

**SepTest** provides statistical tools for testing *first-order
separability* in spatio-temporal point processes, that is, whether a
spatio-temporal intensity function can be expressed as the product of
purely spatial and purely temporal components.

The package implements several hypothesis testing procedures, including
exact and asymptotic methods for both Poisson and non-Poisson processes.
Available methods include global envelope tests, chi-squared-type
statistics, and a novel Hilbert–Schmidt Independence Criterion
(HSIC)–based test using block or pure permutation schemes.

The package covers the simulation studies and applications presented in
Ghorbani et al. (2021, 2025).

## Installing the package

You can install the development version of SepTest from GitHub using:

``` r
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
remotes::install_github("mghorbani01/SepTest")
```

    ## Using GitHub PAT from the git credential store.

    ## Downloading GitHub repo mghorbani01/SepTest@HEAD

    ## 
    ## ── R CMD build ─────────────────────────────────────────────────────────────────
    ##          checking for file 'C:\Users\mdgi0001\AppData\Local\Temp\RtmpeeZKKf\remotes10083d61057\mghorbani01-SepTest-0c57ff3/DESCRIPTION' ...  ✔  checking for file 'C:\Users\mdgi0001\AppData\Local\Temp\RtmpeeZKKf\remotes10083d61057\mghorbani01-SepTest-0c57ff3/DESCRIPTION'
    ##       ─  preparing 'SepTest': (2.4s)
    ##    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
    ##       ─  checking for LF line-endings in source and make files and shell scripts (343ms)
    ##   ─  checking for empty or unneeded directories
    ##       ─  building 'SepTest_0.0.1.tar.gz' (751ms)
    ##      
    ## 

    ## Installing package into 'C:/Users/mdgi0001/AppData/Local/R/win-library/4.5'
    ## (as 'lib' is unspecified)

## Examples

``` r
library(SepTest)

if (requireNamespace("dHSIC", quietly = TRUE)) {

  set.seed(123)

  # Simulated spatio-temporal data
  X <- cbind(
    runif(100),          # spatial x-coordinate
    runif(100),          # spatial y-coordinate
    runif(100, 0, 10)    # temporal component
  )

  # Pure permutation HSIC test
  result <- dHS.test(
    sim.procedure = "pure_per",
    X = X,
    nsim = 199,
    bandwidth = 0.05
  )
  print(result$p.value)

  # Block permutation HSIC test
  result_block <- dHS.test(
    sim.procedure = "block_per",
    X = X,
    nblocks = 5,
    nperm = 100,
    bandwidth = 0.05
  )
  print(result_block$p.value.bw)
}
```

    ## [1] 0.595
    ## [1] 0.4950495

## Documentation

For full documentation of each function, use:

``` r
help(package = "SepTest")
```

## Bug Reports

Please report bugs or feature requests at:

<https://github.com/mghorbani01/SepTest/issues>

When reporting issues, please include a minimal reproducible example and
your sessionInfo() output.

## License

GPL-3
