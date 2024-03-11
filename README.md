# MeDuSAJ

## Installation
```R
install.packages("devtools")

##Please install the "Seurat" first. (https://satijalab.org/seurat/)
install.packages("Seurat")
library(Seurat)

##R version need > 3.5.0
devtools::install_github("LeonSong1995/MeDuSAJ", build_vignettes=F)
library(MLM) # This is the MeDuSAJ package

## Estimate cell state (type) abundance
### required annotation:
sce_use$cellType # your label of cell types or cell states for each cell.
sce_use$sampleID # sample Id for each cell.

## Find marker genes
Idents(sce_use) = sce_use$cellType
mk = FindAllMarkers(sce_merge,only.pos = T,verbose = T,min.pct = 0.1) # you can use the FindAllMarkers functions in Seurat to find marker genes for each cell state (cell type). 
mk$cluster = as.vector(mk$cluster)

mk_gene = sapply(unique(mk$cluster), function(clu){
  mk_temp = mk[mk$cluster==clu,]
  mk_temp[mk_temp$p_val_adj < 0.05,'gene'][1:50]  
})
mk_gene = unique(c(mk_gene))

### run MeDuSAJ
ct.es = CTdcv(bulk = bulk,ncpu = 2,sce = sce_use,data_type = 'count',gene = mk_gene)
```
