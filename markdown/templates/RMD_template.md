---
title: "Template"
author: "Sung Rye Park"
date: "`r format(Sys.Date())`"
output:  
  rmdformats::robobook: 
    code_folding: hide 
    number_sections: TRUE
    toc_depth: 6
    toc_float: true
---

<style> 
#TOC { 
  top: 1%; 
  opacity: 0.5; 
} 
#TOC:hover { 
  opacity: 1; 
} 
</style>

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=TRUE, warning=FALSE, message=FALSE, results = 'asis')
options(warn = F)

library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(DT)
```

**Read me**

This is a template for general markdown report.

The details regarding ____ are currently confidential.

<br>
<hr>
<br>

# Introduction of the project 


# Import data

## Information of data 


# UMAP

# Clustering

# Cell type annotation

# DEG

# GSEA

# AUC

# Additional 1

# Additional 2 

