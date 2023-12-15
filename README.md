# Escort
[![DOI](https://zenodo.org/badge/718898361.svg)](https://zenodo.org/doi/10.5281/zenodo.10392544)

`Escort` is an R package designed to assist users in detecting trajectories, evaluating and ranking various analysis choices. It follows a three-step approach, guiding users through the analysis process by offering goodness-of-fit evaluations for embeddings that represent a range of analysis choices. This encompasses features like feature selection, dimension reduction, and trajectory inference methods, along with their respective hyperparameters.


## Installation

You can install `Escort` via GitHub:

``` r
if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
    
remotes::install_github("xiaorudong/Escort")
```

## Dependencies

The following packages are required for installing Escort. If installation fails, you may need to manually install the dependencies using the function 'install.packages' for CRAN packages or 'BiocManager::install' for Bioconductor packages.


* CRAN packages: 
* Bioconductor packages: 


For any other installation issues/questions please leave a message on our [Github Issues](https://github.com/xiaorudong/Escort/issues).

## Running Escort

After successful installation, `Escort` and related packages need be loaded into the working space:

``` r
library(Escort)
library(cluster)
library(mclust)
library(RMTstat)
```

To use Escort as a Shiny App:

```r
library(shiny)
library(shinydashboard)

Escort::shinyEscort()
```

Many more details and examples on how to use the Shiny App or functions directly in R are located in the vignette:

* [Vignette for R package]()
* [Vignette for Shiny App](vignettes/shiny_vignette.md)


# Contact/Maintainer

* Xiaoru Dong (xdong1@ufl.edu)  
Department of Biostatistics, University of Florida

* [Rhonda Bacher, PhD](https://www.rhondabacher.com) (rbacher@ufl.edu)  
Department of Biostatistics, University of Florida


