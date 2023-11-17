#ChIP-seq work#

# Initial rmarkdown setting 

---
title: "Analysis"
author: "Sung Rye Park"
date: "`r format(Sys.Date())`"
output:
  html_document:
    toc: yes
    code_folding: hide 
    number_sections: FALSE
    toc_depth: 6
    toc_float: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(eval=T, echo=F, warning=FALSE, message=FALSE, results = 'asis')
options(warn = F)

library(cowplot)
library(dplyr)
library(ggplot2)
library(DT)
library(tidyr)
```

<br>
<hr>
<br>

