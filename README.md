# MeDuSAJ
MeDuSAJ supports cell-state deconvolution for annotated `cell states (cell types)`. MeDuSAJ is more robust for estimating cell state (cell type) abundance for rare cells, albeit with a slightly increased computational burden. 


## Installation
in shell command:
```shell
git clone https://github.com/LeonSong1995/MeDuSAJ.git
cd MeDuSAJ
R CMD install `pwd`
```
in R:
```R
install.packages('YourPath/MeDuSAJ',type='source',repos=NULL)
```

## How to use
```R
## Load the example data
example.sce = MeDuSAJ::example.sce
example.bulk = MeDuSAJ::example.bulk

## Estimate cell state (type) abundance
## You can bin your continuous cell trajectory into an uneven discrete cell state. 
### required annotation:
example.sce$cellLabel # your label of cell types or cell states for each cell.
example.sce$sampleID # sample Id for each cell.


## Find marker genes
Idents(example.sce) = example.sce$cellLabel
### You can use the FindAllMarkers functions in Seurat to find marker genes for each cell state (cell type). 
mk = FindAllMarkers(example.sce,only.pos = T,verbose = T,min.pct = 0.1) 
mk$cluster = as.vector(mk$cluster)

mk_gene = sapply(unique(mk$cluster), function(clu){
  mk_temp = mk[mk$cluster==clu,]
  mk_temp[mk_temp$p_val_adj < 0.05,'gene'][1:50]  
})
mk_gene = unique(c(mk_gene))


### Run MeDuSAJ
ct.es = CTdcv(bulk = example.bulk,ncpu = 2,sce = example.sce,data_type = 'count',gene = mk_gene)
```
