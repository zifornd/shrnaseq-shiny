# shrnaseq-shiny

This repository houses the a Shiny app to visualise and explore the output of the shRNAseq workflow of shRNA-seq and CRISPR-Cas9 genetic screen analysis. Before lauching this app you will need to run the workflow on your data. Details on how to run the shrnaseq workflow on your shRNA-seq or CRISPR-Cas9 genetic screening data can be found [here](https://github.com/zifornd/shrnaseq/).

## Usage

The easiest way to run this application is to install the shiny package in RStudio and use the `runGitHub` function as below.

```R
library(shiny)
shiny::runGitHub('shrnaseq-shiny', 'zifornd')
```

Alternatively you can clone the git repository, 

```bash 
git clone https://github.com/zifornd/shrnaseq
```
then use the `runApp()` function in RStudio as below.
```R
library(shiny)
setwd("~/shrnaseq-shiny")
runApp()
```

Once the application is loaded, upload the `shiny.rds` file found in `results` directory of your successfully run shrnaseq workflow. 

## Citations

For more information about Shiny from RStudio, see [here](https://shiny.rstudio.com/). 

### CRAN
- rlang
- shiny
- DT
- ggplot2
- RColorBrewer
- heatmaply
- plotly
- reshape
- tidyverse
- bslib
- scales
- ggrepel

### Bioconductor
- limma
- edgeR