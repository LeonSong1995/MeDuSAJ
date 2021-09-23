
#' @keywords internal
cluster = function(XY,aj){
  currXY=XY
  nbins = aj
  breaks = seq(min(currXY)-10^-7,max(currXY)+10^-7, (max(currXY)-min(currXY)+2*10^-7)/ceiling(nbins))
  cellLocationOnGrid = rep(NA,length(currXY))
  for(currBreakIndex in 2:length(breaks)){
    cellLocationOnGrid[which(currXY>=breaks[currBreakIndex-1] & currXY<breaks[currBreakIndex])]=currBreakIndex-1 
  }
  names(cellLocationOnGrid) = rownames(currXY)
  return(cellLocationOnGrid)
}


#' @keywords internal
geneAsso = function(space,exprsData,mode,maxgene,cov=NULL,family,k){  
  

  CellCluster = sapply(names(table(mode)),function(id){names(which(mode ==id))})
  
  #rename cells
  cellNamesForAnova = unlist(sapply(seq(1:length(CellCluster)),function(cell){
    paste0(cell,rep('c',length(CellCluster[[cell]])))
  }))

  #Cluster-level profile
  groupRef = t(aggregate(t(exprsData),by=list(mode),FUN=mean)[,-1])
  colnames(groupRef) = paste0(unique(mode),'c')
  genes_to_take = rownames(groupRef)
  cellNum = ncol(exprsData)
  
  ####cov##########
  U = rep(1,cellNum)
  if(!is.null(cov)){
    U = cov
  }

  #ANOVA
  dat = cbind(rep(0,length(genes_to_take)),exprsData[genes_to_take,])
  group = c("",cellNamesForAnova)
  dmat <- stats::model.matrix(~ group)
  fit <- limma::eBayes(limma::lmFit(dat, dmat)[,-1])
  fitF = fit$F
  names(fitF) = genes_to_take
 
  res = data.frame(colnames(groupRef)[apply(groupRef,1,which.max)],fitF)
  colnames(res)= c("group", "F_score")
  rownames(res)= rownames(groupRef)
  
  
  
  allGenes = sapply(unique(cellNamesForAnova),function(bin){
    index = which(as.character(res$group)==as.character(bin))
    temp = res[index,]
    temp = temp[order(as.numeric(temp$F_score),decreasing = T),]
    rownames(temp)
  })
  


  #select genes
  geneNumberCluster = ceiling(maxgene/ncol(groupRef))
  
  topGenes = unlist(sapply(1:length(allGenes),function(i){
    temp = allGenes[[i]]
    na.omit(temp[1:min(geneNumberCluster,length(temp))])
  }))
  topGenes = unlist(c(topGenes))

  #remove the over-weighted genes
  # expM = apply(exprsData[topGenes,],1,mean)
  # topGenes=topGenes[-which(expM>(mean(expM)+2*sd(expM)))]
  
  #Choose gene 
  # k = c()
  # bestKappa = Inf
  # bestG = 0
  # bestGenes = c()
  # for (i in  maxgene:10){
  #   selectedGenes = allGenes[1:i]
  #   currgroupRef = groupRef[match(selectedGenes, row.names(groupRef)),]
  #   newKappa = kappa(currgroupRef)
  #   k = c(k,newKappa)
  #  print(newKappa)
  #   
  #   if (newKappa <= bestKappa){
  #     bestKappa = newKappa
  #     bestG = i
  #     bestGenes = unique(selectedGenes)
  #   }
  # }
  # 
  
  return(list('bestGenes'=topGenes,'gall'=allGenes,'F'=fitF))
}

#' @keywords internal
VariableGenes = function(a, ratio) {
  count_nonZeros = length(which(a > min(a)))
  if (count_nonZeros/length(a) > ratio) {TRUE} else{FALSE}
}

#' @keywords internal
geneSelect = function(exprsData,cellType,space,bulk,maxgene,aj,cov=NULL,family='Gaussian',k=10){
 
   #Choose cell-expressed genes
  eligibleGenes = apply(exprsData,1,function(gene){VariableGenes(as.numeric(as.matrix(gene)),0.1)})
  exprsData = exprsData[eligibleGenes,]
  
  commgene = intersect(rownames(exprsData),rownames(bulk))
  exprsData = exprsData[commgene,]
  print(dim(exprsData))
  
  mode = cluster(space, aj)
  gw = geneAsso(space=space,exprsData=exprsData,mode=mode,maxgene=maxgene,cov=cov,family,k=k)
  g = gw$bestGenes
  F_va = gw$F
  
  return(list('g'=g,'gall'=gw$gall,'Fv'=F_va))
}


