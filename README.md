# vissE - Visualising Set Enrichment Analysis Results <img src="https://github.com/Bioconductor/BiocStickers/blob/devel/vissE/vissE.png" alt="logo" align="right" height="140" width="120"/>

[![R-CMD-check](https://github.com/DavisLaboratory/vissE/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/DavisLaboratory/vissE/actions)
[![codecov](https://codecov.io/gh/DavisLaboratory/vissE/branch/main/graph/badge.svg?token=8JHZB1GN26)](https://codecov.io/gh/DavisLaboratory/vissE)
[![BioC status](https://bioconductor.org/shields/years-in-bioc/vissE.svg)](https://bioconductor.org/packages/vissE/)

This package implements the vissE algorithm to summarise results of gene-set analyses. Usually, the results of a gene-set enrichment analysis (e.g using limma::fry, singscore or GSEA) consist of a long list of gene-sets. Biologists then have to search through these lists to determines emerging themes to explain the altered biological processes. This task can be labour intensive therefore we need solutions to summarise large sets of results from such analyses.

This package provides an approach to provide summaries of results from gene-set enrichment analyses. It exploits the relatedness between gene-sets and the inherent hierarchical structure that may exist in pathway databases and gene ontologies to cluster results. For each cluster of gene-sets vissE identifies, it performs text-mining to automate characterisation of biological functions and processes represented by the cluster.

An additional power of vissE is to perform a novel type of gene-set enrichment analysis based on the network of similarity between gene-sets. Given a list of genes (e.g. from a DE analysis), vissE can characterise said list by first identifying all other gene-sets that are similar to it, following up with clustering the resulting gene-sets and finally performing text-mining to reveal emerging themes.

In addition to these analyses, it provides visualisations to assist the users in understanding the results of their experiment.

[**Check out the full tutorial**](https://davislaboratory.github.io/vissE/articles/vissE.html)

## Installation

vissE can be installed from Bioconductor directly as follows:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("vissE")
```
