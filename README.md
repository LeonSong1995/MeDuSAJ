# MeDuSAJ

## Installation
```R
install.packages("devtools")

##Please install the "Seurat" first. (https://satijalab.org/seurat/)
install.packages("Seurat")
library(Seurat)

##R version need > 3.5.0
devtools::install_github("LeonSong1995/MeDuSAJ", build_vignettes=F)

## Estimate cell state (type) abundance
### required annotation:
sce_use$cellType
sce_use$sampleID
### run MeDuSAJ
ct.es = CTdcv(bulk = bulk,ncpu = 2,sce = sce_use,data_type = 'count')


```
