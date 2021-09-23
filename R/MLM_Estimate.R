# Deep Deconvolution
# Author: Liyang Song <songliyang@westlake.edu.cn>
# Adviser: Jian Yang, Xiwei Sun
# Copyright: Jian Yang
#############################################################################################################
#' Cell-type level deconvolution
#'
#' @param bulk A matrix containing bulk RNA-Seq data. Each row corresponds to a certain gene and each column to a certain sample.
#' @param sce A 'Seurat' object containing the single-cell RNA-Seq data. Meta data of the 'Seurat' object must includes 'cellType' and 'sampleID'.
#' @param gene A character vector of the gene names to use as signature for the deconvolution. We summarized signature genes for the 64 human cell-types. please see 'sg'.
#' @param data_type A character of the type of the single-cell RNA-Seq data, including 'count', 'tpm', 'rpkm','fpkm', and 'cpm'.
#' @param select.ct A character vector of the names of the target cell-types. The default value is NULL. With default value, all cell-types in the single-cell data will be used.
#' @param RanSplit A character (or character vector) to split the random components. The default value is NULL. With default value, all cells, excepting those in target the cell-type, will be fitted as one random component.
#' @param ct.cell.size A character vector of the cell-size (total mRNA amount) of the selected cell-types. The default value is NULL.
#' @param BatchCorrect A Boolean variable to determine whether to run 'ComBat' to correct batch effects between single-cell RNA-Seq and bulk RNA-Seq or not. The default value is FALSE.
#' @param Filter A Boolean variable to determine whether to filter the outlier in single-cell data or not. The default value is FALSE.
#' @param SF Scaling factor. The default value is 1e+3.
#' @param ncpu The number of CPU cores to be used.
#' @param iter_max The maximum iterations of REML.
#'
#' @return A list with elements:
#'   *ct.pro: matrix of cell-type proportions estimated by mixed linear model (sample x cell-type);
#'   *ct.pro.p: matrix of p value (χ²(df=1)) for the cell-type proportions estimated by the mixed linear model (sample x cell-type);
#'   *cellSize: vector of cell sizes with labeled cell-type names.
#'
#' @export
#'
#' @examples
#' library(MLM)
#' ct.es = CTdcv(bulk = example.bulk,sce = example.sce,gene = example.gene,data_type = 'count')
#'
CTdcv = function(bulk,sce,gene=NULL,data_type, select.ct = NULL, RanSplit=NULL, ct.cell.size = NULL,BatchCorrect=F,Filter=T,SF=1e+3,ncpu=NULL,iter_max=1000){

	print("Thanks for using MLM to perform bulk deconvolution analysis.")

	#Checking the correct format of the reference single-cell data input
	if (!(data_type %in% c('count','tpm','rpkm','cpm','fpkm'))){
		stop('Please input correct data_type: count/tpm/rpkm/cpm/fpkm.')
	}
	if(!("Seurat" %in% class(sce))){
		stop('Please input Seurat Object.')
	}
	if(!is.null(RanSplit)){
		if (is.na(match(RanSplit,colnames(sce@meta.data)))){
		stop('Do not know how to split randomp components, please check your MetaData: Seurat_Obj@meta.data.')
		}
	}
	if(!'sampleID' %in% colnames(sce@meta.data)){
		stop('Please input sampleID, check your MetaData: Seurat_Obj@meta.data.')
	}
	if(!'cellType' %in% colnames(sce@meta.data)){
		stop('Please input cellType, check your MetaData: Seurat_Obj@meta.data.')
	}
	if(length(unique(sce$cellType))==1){
		stop('MLM can not work with only one cell type inputted.')
	}
	if(!is.null(select.ct)){
		if(is.na(table(sce$cellType %in% select.ct)['TRUE'])){
			stop('No cell types selected. Please check the select.ct!')
		}else{
			sce = sce
		}
	}else{
		select.ct = unique(sce$cellType)
	}

	MetaData = sce@meta.data
	exprsData = as.matrix(sce@assays$RNA@counts)
	bulk = bulk[rowSums(bulk)>0,]

	#finding signature gene
	if(is.null(gene)){
		print('Finding signature genes with "FindMarkers".')
		gene = SignatureGenerator(sce[,sce$cellType %in% select.ct])
		print(length(gene))
	}

	#Checking the signature genes input
	commonGene = intersect(rownames(exprsData),rownames(bulk))
	gene = intersect(commonGene,gene)
	if(length(gene)<10){
		stop('Too few signature genes (signature genes < 10).')
	}

	bulk = as.matrix(bulk)[commonGene,]
	exprsData = exprsData[commonGene,]


	MetaData$cellType = as.vector(MetaData$cellType)
	MetaData$sampleID = as.vector(MetaData$sampleID)

	#Preparing for the basic running information (fixed/random components, cell size...)
	print("Data preparing.")
	Info = basis(bulk = bulk, exprsData = exprsData, MetaData=MetaData, ct.cell.size = ct.cell.size,
				data_type=data_type,gene=gene,BatchCorrect = BatchCorrect,Filter=Filter,SF=SF)


	base = Info$base[,select.ct]
	bulk = Info$bulk
	data_cellType = Info$data_cellType
	cellSize = Info$cellSize
	type_n = length(data_cellType)


	#Preparing for the random components
	rancmp = sapply(seq(1,type_n), function(i){
				Z_new = data_cellType
				Z_new[i] = NULL
				Rest =  matrix(unlist(Z_new),nrow = length(gene))
				cellName = unlist(lapply(Z_new,colnames))
				colnames(Rest) = cellName
				Rest})

	rancmp = sapply(seq(1,type_n)[colnames(base) %in% select.ct],function(i){rancmp[[i]]})

	ct_name = colnames(base)
	bk_name = colnames(bulk)


	#Run cell-type level deconvolution parallelly
	if(is.null(ncpu)){
		ncpu = max(1, parallel::detectCores() - 1)
	}
	print(paste(paste("Running cell-type level deconvolution with",ncpu,sep=" "),'cores.',sep=" "))
	cl = parallel::makeCluster(ncpu)
	parallel::clusterExport(cl=cl, varlist=c("RunReml", "base", "bulk","rancmp","MetaData","RanSplit","iter_max", "bk_name","reml"),
			  envir=environment())
	doSNOW::registerDoSNOW(cl)
	pb = utils::txtProgressBar(min = 1, max = ncol(bulk), style = 3)
	progress = function(n) setTxtProgressBar(pb, n)
	opts = list(progress = progress)
	`%dopar2%` = foreach::`%dopar%`
	runNumber = NULL
	estimate = foreach::foreach(runNumber = 1:ncol(bulk), .options.snow = opts) %dopar2% {

		r = RunReml(runNumber,base=base,bulk=bulk,rancmp=rancmp,MetaData=MetaData, RanSplit=RanSplit, iter_max,bk_name)
		b = r['b',]
		p = r['p',]
		setTxtProgressBar(pb, runNumber)
		cbind(b,p)
	}
	parallel::stopCluster(cl)
	close(pb)
	b = t(sapply(1:length(estimate), function(o){estimate[[o]][,1]}))
	p = t(sapply(1:length(estimate), function(o){estimate[[o]][,2]}))

	colnames(b) = colnames(p) = ct_name
	rownames(b) = rownames(p) = bk_name


	#Adjusting for the cell size when count matrix inputted
	if(data_type=='count'){
		b = sweep(b,2,cellSize[colnames(b)],'/')
		b = sweep(b,1,rowSums(b),'/')
	}else{
		warning("The estimated cell type proportions are not comparable among different cell types!")
	}

	#Adjusting the negative value
	if(min(b)<0){b = b + abs(min(b))}

	return(list('ct.pro'=b,'ct.pro.p'=p,'cellSize'= cellSize))
}


#' @keywords internal
RunReml = function(sid,base, bulk, rancmp, MetaData, RanSplit, iter_max, bk_name){
  x = as.matrix(base)
  y = bulk[,sid]
  type_n = ncol(x)

  result = sapply(seq(1,type_n),function(id){
    fixcmp = as.matrix(x[,id])
    if (!is.null(RanSplit)){
      # Splitting random components based on the inputted RanSplit
      sepInfo = MetaData[intersect(colnames(rancmp[[id]]),rownames(MetaData)),RanSplit]
      Z = sapply(unique(sepInfo),function(sid){
        rancmp[[id]][,sepInfo %in% sid]
      })
    }else{
      Z = list(rancmp[[id]])
    }
	start = rep(var(y)/(length(Z)+1),(length(Z)+1))
    mlmfit = reml(start,X = fixcmp,y = y,Z = Z,maxiter = iter_max)
	b = mlmfit[[1]]
	p = 1-pchisq((b*b)/diag(mlmfit[[2]]),df=1)
	varcmp = mlmfit[[3]]

    return(c('b'=b,'p'=p))
  })

  colnames(result) = colnames(x)
  return(result)
}



####################################################################################################
#' Single-cell level deconvolution
#' @param bulk A matrix containing bulk RNA-Seq data. Each row corresponds to a certain gene and each column to a certain sample.
#' @param sce A 'Seurat' object containing the single-cell RNA-Seq data. Meta data of the 'Seurat' object must includes 'cellType' and 'space'.
#' @param select.ct A character (or character vector) of the names of the target cell-types. The default value is NULL. With default value, all cell-types in the single-cell data will be used.
#' @param ncpu The number of CPU cores to be used.
#' @param smoothing A Boolean variable to determine whether to smooth the BLUP along cell-space or not. The default value is TRUE.
#' @param gene A character vector of the gene names to use as signature for the deconvolution.The default value is NULL. With default, MLM will select genes to differentiate cells in the target cell-type with ANOVA.
#' @param SF Scaling factor. The default value is 1e+3.
#'
#' #Gene selecting
#' @param op A numeric variable to determine the overlapping between cell-clusters. The default value is 0.2.
#' @param maxgene A numeric variable to determine the maximum number of genes to be selected. The default value is 1e+3.
#' @param nbin A numeric variable to determine the number of cell-clusters. The default value is 0.2.
#'
#' @return A list with cell_density-cell_name matrix:
#'
#' @export
#'
#' @examples
#' library(MLM)
#' SCdcv(bulk = example.bulk,sce = example.sce,select.ct = 'alpha')
SCdcv = function(bulk,sce,select.ct,ncpu=NULL,smoothing=TRUE,gene=NULL,aj=0.2,maxgene=1000,SF=1e+3){

  print("Thanks for using MLM to perform bulk deconvolution analysis.")

	#Checking the correct format of the reference single-cell data input
	if(!("Seurat" %in% class(sce))){
		stop('Please input Seurat Object.')
	}
	if(!'cellType' %in% colnames(sce@meta.data)){
		stop('Please input cellType, check your MetaData: Seurat_Obj@meta.data.')
	}
	if(!'space' %in% colnames(sce@meta.data)){
		stop('Please input space, check your MetaData: Seurat_Obj@meta.data.')
	}
	if(!is.null(select.ct)){
		if(is.na(table(sce$cellType %in% select.ct)['TRUE'])){
			stop('No cell types selected. Please check the select.ct!')
		}
	}else{
		select.ct = unique(sce$cellType)
	}

	#Preparing the reference
	space = as.matrix(sce$space)
	cellType = as.matrix(sce$cellType)
	exprsData = as.matrix(sce@assays$RNA@counts)
	bulk = bulk

	#Scaling
	exprsDataS = as.matrix(sweep(exprsData,2,colSums(exprsData)+1e-200,'/')*SF)
	bulkS = as.matrix(sweep(bulk,2,colSums(bulk)+1e-200,'/')*SF)

	result = list()
	i = 1
	for (ct in select.ct){
		currCT =names(cellType[which(cellType==ct),])
		currSpace = as.matrix(sort(space[currCT,]))
		currExpS = as.matrix(exprsDataS[,rownames(currSpace)])
		# sd = apply(currExpS,1,sd)
		# currExpS = sweep(currExpS,1,sd,'/')
		currType = as.matrix(cellType[rownames(currSpace),])


		if(is.null(gene)){
			#Choose genes
			print('selecting genes')
			gm= geneSelect(exprsData=currExpS,cellType=currType,space=currSpace, bulk=bulk,aj = aj,maxgene)
			gene = gm[[1]]
			mode = cluster(currSpace,0.5)
		}

		gene = intersect(intersect(gene,rownames(bulk)),rownames(exprsData))
		print(paste(length(gene),'genes were selected to do deep deconvolution.'))

		#Generate fix components
		if (length(unique(cellType))==1){
			fixcmp = as.matrix(rep(1,length(gene)))
		}else{
			fixcmp = sapply(unique(cellType)[-which(unique(cellType)==ct)], function(ty){
				temp = as.matrix(exprsDataS[gene,which(cellType==ty)])
				if(ncol(temp)>1){
					temp = rowMeans(temp)
				}
				temp
			})
			fixcmp = as.matrix(cbind(fixcmp, rep(1,length(gene))))
		}
		rownames(fixcmp) = gene
		# rancmp = sapply(unique(mode),function(cell){
		# 		list(as.matrix(currExpS[gene,which(mode==cell)]))
		# 	})
		rancmp = list(as.matrix(currExpS[gene,]))

		currbulk = bulkS[gene,]

		#Run single-cell level deconvolution parallelly
		if(is.null(ncpu)){
			ncpu = max(1, parallel::detectCores() - 1)
		}
		print(paste(paste("Running single-cell level deconvolution with",ncpu,sep=" "),'cores.',sep=" "))
		cl = parallel::makeCluster(ncpu)
		parallel::clusterExport(cl=cl, varlist=c("reml2","reml", "fixcmp", "rancmp",'currbulk'),
				  envir=environment())
		doSNOW::registerDoSNOW(cl)
		pb = utils::txtProgressBar(min = 1, max = ncol(currbulk), style = 3)
		progress = function(n) setTxtProgressBar(pb, n)
		opts = list(progress = progress)
		`%dopar2%` = foreach::`%dopar%`
		runNumber = NULL
		estimate = foreach::foreach(runNumber = 1:ncol(currbulk), .options.snow = opts) %dopar2% {
			setTxtProgressBar(pb, runNumber)
			y <- currbulk[,runNumber]
			# e = c()
			# sd2 = c()
			# start = rep(var(y)/2,2)
			# for (j in 1:ncol(currExpS)){
			# 	indexTemp = which(mode == mode[j])
			#   fixcmpTemp = cbind(currExpS[gene,j],fixcmp)
			#   rancmp = as.matrix(currExpS[gene, -indexTemp])
			#   r =reml(start = start,X = fixcmpTemp,Z = list(rancmp),y = y,maxiter = 1e+3)
			#   e = c(e,r[[1]][1])
			#   sd2 = c(sd2,diag(r[[2]])[1])
			#   start = r[[3]]
			# }
			# P_value=1-pchisq((e*e)/sd2,df=1)

			Vi = reml2(X = fixcmp,y = y,Z = rancmp,maxiter=1e+4,start=rep((var(y)/2),length(rancmp)+1))
			e = Vi[[1]][1,]
			z = e/(Vi[[1]][2,])
			P_value=1-pchisq(z^2,1)
			rbind(e,P_value)
		}
		parallel::stopCluster(cl)
		close(pb)
		deepcp = t(sapply(estimate,function(e){e[1,]}))
		p = t(sapply(estimate,function(e){e[2,]}))
		colnames(deepcp) = colnames(p) = colnames(currExpS)
		rownames(deepcp) = rownames(p) = colnames(bulk)

		#Smoothing
		chosenNeigList = lapply(1:nrow(currType),function(cellIndex){
			cellDist = abs(currSpace[cellIndex,]- currSpace)
			chosenRepeats = order(as.numeric(cellDist),decreasing = F)[1:10]
			chosenRepeats = chosenRepeats[!is.na(chosenRepeats)]
			names(currType[chosenRepeats,])
		})

		# for (id in 1: nrow(deepcp)){
		# 	for (cell in 1:ncol(deepcp)){
		# 		temp = deepcp[id,chosenNeigList[[cell]]]
		# 		m =lm(temp~currSpace[names(temp),])
  # 			o = outlierTest(m)
		# 		h =  mean(temp)+2*sd(temp)
		# 		l = mean(temp)-2*sd(temp)
		# 		if(names(temp)[1] %in% names(o$rstudent)){
		# 			deepcp[id,cell] = mean(temp)
		# 		}
		# 	}
		# }
		# colnames(deepcp) = colnames(p) = colnames(currExpS)
		# rownames(deepcp) = rownames(p) = colnames(bulk)

		# deepcp = t(do.call(rbind,lapply(1:length(chosenNeigList),function(cell){
		#   rowMeans(deepcp[,chosenNeigList[[cell]]])
		#   })))

		if(smoothing){
		# deepcp = t(do.call(rbind,sapply(unique(mode),function(cell){
  # 		list(sapply(1:nrow(deepcp),function(r){
  #   	predict(loess(deepcp[r,which(mode==cell)]~currSpace[which(mode==cell),]))
  # 		}))
		# })))
		# colnames(deepcp)  = colnames(currExpS)
		# rownames(deepcp)  = colnames(bulk)
		deepcp = t(sapply(1:ncol(bulk),function(i){predict(loess(deepcp[i,]~currSpace,span = 0.75))}))

		}
		colnames(deepcp) = colnames(currExpS)
		rownames(deepcp) = colnames(bulk)




		#Scaling the result
		# deepcp = deepcp + abs(min(deepcp))
		# deepcp = sweep(deepcp,1,rowSums(deepcp),'/')

		result[[i]] = deepcp
		i = i+1
		}

	return(result)
}
