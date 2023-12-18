# Escort
[![DOI](https://zenodo.org/badge/718898361.svg)](https://zenodo.org/doi/10.5281/zenodo.10392544)

<img src="https://github.com/xiaorudong/Escort/blob/main/vignettes/escort_hex_v2.png"  width = "120">

`Escort` is a framework for evaluating a single-cell RNA-seq dataset’s suitability for trajectory inference and for quantifying trajectory properties influenced by analysis decisions. Escort is an R package designed to guide users through the trajectory inference process by offering goodness-of-fit evaluations for embeddings that represent a range of analysis decisions such as feature selection, dimension reduction, and trajectory inference method-specific hyperparameters.


## Installation

You can install `Escort` via GitHub:

``` r
if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
    
remotes::install_github("xiaorudong/Escort")
```

## Dependencies

The following packages are required for installing Escort. If installation fails, you may need to manually install the dependencies using the function 'install.packages' for CRAN packages or 'BiocManager::install' for Bioconductor packages.

* Bitbucket: 
``` r
    remotes::install_bitbucket("scLCA/single_cell_lca")
```
* CRAN packages: alphahull, fastICA, Rtsne, umap, rstatix, DT, FNN, grDevices, Hmisc, jmuOutlier, shotGroups, parallelDist, robustbase, scales, SCORPIUS, sfsmisc
shiny, shinycssloaders, shinydashboard, shinyjs, shinyWidgets, dplyr
* Bioconductor packages: SingleCellExperiment, SC3, slingshot, clusterProfiler, org.Hs.eg.db, org.Mm.eg.db, edgeR, limma

    
For any other installation issues/questions please leave a message on our [Github Issues](https://github.com/xiaorudong/Escort/issues).

## Running Escort

After successful installation, `Escort` and related packages need be loaded into the working space:

``` r
library(Escort)
```

To use Escort as a Shiny App:

```r
library(shiny)
library(shinydashboard)

Escort::shinyEscort()
```

Many more details and examples on how to use the Shiny App or functions directly in R are located in the vignette:

* [Vignette for R package](https://www.rhondabacher.com/vignettes/Escort.html)
* [Vignette for Shiny App](https://www.rhondabacher.com/vignettes/shiny_vignette.html)


# Contact/Maintainer

* [Xiaoru Dong](https://xiaorudong.github.io) (xdong1@ufl.edu)  
Department of Biostatistics, University of Florida

* [Rhonda Bacher, PhD](https://www.rhondabacher.com) (rbacher@ufl.edu)  
Department of Biostatistics, University of Florida


