# Base function
# Author: Liyang Song <songliyang@westlake.edu.cn>
# Advisor: Jian Yang, Xiwei Sun
# Copyright: Liyang Song

#############################################################################################

#' @keywords internal
basis = function(exprsData,MetaData,bulk,ct.cell.size,data_type, Filter, BatchCorrect, SF, gene){

	cell_name = colnames(exprsData)
	countmat = exprsData
	ct.id = MetaData[cell_name,]$cellType
	sample.id = MetaData[cell_name,]$sampleID
	ct_sample.id = paste(ct.id,sample.id, sep = '%')

	if (data_type == "count"){
	mean.mat = sapply(unique(ct_sample.id), function(id){
	  y = as.matrix(countmat[, ct_sample.id %in% id])
	  apply(y,1,sum, na.rm = TRUE)/sum(y)
	})
	}else{
	mean.mat = sapply(unique(ct_sample.id), function(id){
	  y = as.matrix(countmat[, ct_sample.id %in% id])
	  apply(y,1,mean)
	})
	}
	# GeneIM = Impute(exprsData,MetaData,gene,SF)
	

	mean.id = do.call('rbind',strsplit(unique(ct_sample.id), split = '%'))

	sum.mat2 = sapply(unique(sample.id), function(sid){
	sapply(unique(ct.id), function(id){
	  y = as.matrix(countmat[, ct.id %in% id & sample.id %in% sid])
	  sum(y)/ncol(y)
	})
	})
	rownames(sum.mat2) = unique(ct.id)
	colnames(sum.mat2) = unique(sample.id)


	# library size factor calculated from the samples:
	if (is.null(ct.cell.size)){
	sum.mat = rowMeans(sum.mat2, na.rm = T)
	} else {
	if (is.null(names(ct.cell.size))){
	  stop("Cell size factor vector requires cell type names...")
	} else {
	  sum.mat = ct.cell.size
	}}


	base = sapply(unique(mean.id[,1]), function(id){
		y = as.matrix(mean.mat[,mean.id[,1] %in% id])
		y = apply(y,1,mean, na.rm = TRUE)
		y = (y / (sum(y)+1e-200))*SF
		return(y)
	})
	# base = sapply(GeneIM,rowMeans)
	# rownames(base) = gene
	
	# Adding an extremely small positive number (1e-200) to avoid infinity.
	data_cellType = sapply(unique(ct.id), function(id){
		y = countmat[, ct.id %in% id]
		y = sweep(y,2,colSums(y)+1e-200,'/')*SF
		if (Filter==TRUE){
			y = cellFilter(y[gene, ])
		}
		return(y[gene, ])
	 })
	
	# data_cellType = sapply(GeneIM,function(ct){
		# A = cellFilter(ct)
		# rownames(A) = gene
		# return(A)
	# })

	if(is.null(ncol(bulk))){
		bulk = ((bulk/(sum(bulk)+1e-200))*SF)[gene]
		if(BatchCorrect==TRUE){warning('Sample Size too samll (<20 samples), cannot remove BatchEffect!')}
	} else{
		if (BatchCorrect==TRUE){
			if(ncol(bulk)<20){
				warning('Sample Size too samll (<20 samples), cannot remove BatchEffect!')
				bulk = (sweep(bulk,2,colSums(bulk)+1e-200,'/')*SF)[gene,]
			}else{
				message('Removing BatchEffect with ComBat...')
				bulk = (BatchCorrection(exprsData,MetaData,bulk,base,SF,gene))[gene,]
			}
		}else{
			bulk = (sweep(bulk,2,colSums(bulk)+1e-200,'/')*SF)[gene,]
		}
	}


	return(list(base = base[gene,],
				data_cellType = data_cellType,
				bulk = bulk,
				cellSize = sum.mat))
}


#' @keywords internal
cellFilter<-function(exprsData){
  d <- as.matrix(dist(log(t(exprsData+1))))
  diag(d) <- 1e+10
  d_near<- apply(d,2,min)
  d_q <- quantile(d_near,0.25) + 1.5*(quantile(d_near,0.75)-quantile(d_near,0.25))
  return(exprsData[,which(d_near<d_q)])
}


#' @keywords internal
BatchCorrection = function(exprsData,MetaData,bulk,basis,SF,gene){

	exprsData = sweep(exprsData,2,colSums(exprsData)+1e-200,'/')*SF
	bulk = sweep(bulk,2,colSums(bulk)+1e-200,'/')*SF
	basis = sweep(basis,2,colSums(basis)+1e-200,'/')*SF

	cell_name = colnames(exprsData)
	ct.id = MetaData[cell_name,]$cellType
	sample.id = MetaData[cell_name,]$sampleID

	cell_number = table(MetaData$sampleID)

	#Estimating cell type proportion by Non-negative least squares (nnls)...
	PesudoFraction = t(apply(bulk[gene,],2,function(b){
				nnls::nnls(basis[gene,],b)$x
				}))
	colnames(PesudoFraction) = colnames(basis)

	cell_select = 1e+3 * round(sweep(PesudoFraction,1,rowSums(PesudoFraction),'/'),3)

	#Constructing single-cell data-based artificial bulk RNA-Seq data by the nnls estimated cell-type proportions...
	PesudoSample = apply(cell_select,1,function(psid){
						Pesudo = sapply(unique(ct.id),function(id){
									s = psid[id]
									exprsData[gene,sample(cell_name[ct.id %in% id],s,replace=T)]
								})
						rowMeans(matrix(unlist(Pesudo),nrow=length(gene)))
					})
	rownames(PesudoSample) = gene

	Merge = log2(cbind(bulk[gene,],PesudoSample[gene,])+1)
	Batch = c(rep('B',ncol(bulk)),rep('P',ncol(PesudoSample)))
	Merge_new = sva::ComBat(dat=Merge[gene,], batch=Batch, par.prior=TRUE, prior.plots=F)
	Merge_new = as.matrix(2^Merge_new - 1)

	BulkTransfer = pmax(Merge_new[,1:ncol(bulk)],0)
	return(BulkTransfer = BulkTransfer)
}
