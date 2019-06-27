# proteomicsliteratureminer
proteomicsliteratureminer is a new tool that aims to help researchers reduce time spent on literature research post analysis and streamline the decision about which proteins or genes are the most interesting and most promising for follow-up experiments.

<!-- badges: start -->
<!-- badges: end -->

The goal of proteomicsliteratureminer is to streamline the process of literature retrieval and provides result categorisation to assist researchers select appropraite leads for further research.

## Installation

You can install the released version of proteomicsliteratureminer from [SIH-GIT](https://github.com/Sydney-Informatics-Hub/proteomicsliteratureminer) with:

``` r
install.packages(devtools) # only if devtools is not installed
install_github("Sydney-Informatics-Hub/proteomicsliteratureminer")
```

## Example

This is a basic example which shows you how to create an Excel file which has the Pubmed results. potentialmarker is a R dataframe that is part of the R package, for description of its contents, run the following R command
``` r
library(proteomicsliteratureminer)
?potentialmarker
```

``` r
library(proteomicsliteratureminer)
pubmedMiner_entry(potentialmarker, output.file = "potential_marker_pubmed_results.xlsx")
```
## Citing
When using Lit-helper please cite: Steffen P, Wu J, Hariharan S, Molloy MP, Schluter H, Lit-helper A bioinformatics tool for prioritizing biological leads from omics data using literature mining.
