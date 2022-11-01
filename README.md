# shrnaseq-shiny

Shiny application for the output of the [shRNAseq](https://github.com/zifornd/shrnaseq) snakemake workflow of shRNA-seq and CRISPR-Cas9 genetic screen analysis using edgeR.

The easiest way to run this application is to install the shiny package in RStudio and use the `runGitHub` function as below.

```R
library(shiny)
shiny::runGitHub('shrnaseq-shiny', 'zifornd')
```

Alternatively you can clone the git repository, then use the `runApp()` function in RStudio as below.

```bash 
git clone https://github.com/zifornd/shrnaseq
```

```R
library(shiny)
setwd("~/shrnaseq-shiny")
runApp()
```

Once the application is loaded, upload the `shiny.rds` file found in `results` directory of your successfully run snakemake workflow.
