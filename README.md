# Escort

## Installation

You can install `Escort` via GitHub:

``` r
remotes::install_github("xiaorudong/Escort")
```

## Libraries

After successful installation, `Escort` and related packages need be loaded into the working space:

``` r
library(Escort)
library(cluster)
library(mclust)
library(RMTstat)
```

## Running Escort

To load the package into R:

```{r}
library(Escort)
```

To use Escort as a Shiny App:
```{r}
Escort::shinyEscort()
```
