## initial setting
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=TRUE, fig.align = "center", message=F, warning=F, out.width='80%', fig.asp=.75, fig.cap='Figure')
```
## add line
---
## load packages
```{r, echo=F}
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(knitr)
library(kableExtra)
```
---
## data import
```{r}
dir='~/Desktop/DF/DFCI_Paweletz/2023_Takeda/'
obj.srt = readRDS(paste0(dir, 'rds/MT.all.raw.singlets.23.05.16.rds'))
```
---
## print table with DT::datatable()
```{r, echo=F}
obj.srt@meta.data$orig.ident %>% table() %>% data.frame()  %>% DT::datatable(class = 'cell-border stripe', colnames = c('sample', 'number of cells'), rownames = F)
obj.srt@meta.data$predicted_doublet %>% table() %>% data.frame() %>% DT::datatable(colnames = c('doublet', 'number of cells'), rownames = F)
```
---
## print table with print()

```{r, echo=F}
averaged_data <- obj.srt@meta.data %>%
  group_by(orig.ident, sample) %>%
  summarise(average_of_UMI = mean(nCount_RNA)) %>% tidyr::spread(sample, average_of_UMI)
print(averaged_data)
```
---
## print table with kable 
```{r results='asis', echo=F}
# Create a data frame with the sample information
sample_info <- data.frame(
  original_data = c("30212", "30268", "30269", "30271", "30273"),
  sample = c("MT30212", "MT30268", "MT30269", "MT30271", "MT30273"),
  condition = c("C,P, T,COM", "C,P, T,COM", "C,P, T,COM", "C,P, T,COM", "C,T,COM"),
  response = c("R", "R", "NR", "R", "NR"),
  number_of_cells = c(3411, 95597, 28651, 24425, 11871),
  COM = c(867, 23687, 4553, 3646, 5173),
  CTL = c(883, 24149, 7118, 6091, 3929),
  PEM = c(1005, 20001, 10693, 7853, 0),
  TAK = c(656, 27760, 6287, 6835, 2769)
)

# Display the table using kable with center alignment
kable(sample_info, format = "html", caption = "Sample Information Table") %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE, latex_options = "hold_position=false")
```
---
## Include image from the storage
```{r,out.width = '80%', echo=F}
dir <- '~/Desktop/DF/DFCI_Paweletz/2023_Takeda/'
knitr::include_graphics(paste0(dir,'data/all/fig/Fig1.all_samples_number_of_cells.png'))
```
---
