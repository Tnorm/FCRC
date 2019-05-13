library(minfi)
library(minfiData)
library(sva)
library(devtools)
library(MASS)
library(cluster)
library(mclust)
library(limma)
library(caret)
library(kernlab)
library(sparsereg)
library(randomForest)
library(rpart)
library(e1071)
library(grid)
library(gridExtra)
library(psych)
library(uniftest)
library(FisherEM)
library(vioplot)
library(entropy)
library(fastICA)
library(tclust)
library(FNN)
library(EntropyEstimation)
library(ClustOfVar)
library(corrplot)

FCRC <- function(
  data,
  annotation,
  phenoData,
  data_dir='',
  preprocess_method = 'raw',
  preprocess_funcs = list(),
  num_var_clusts = 50,
  num_clusts = 0,
  absorption_thresh = 1.0,
  IQR_filtered_thresh = 1.0,
  final_clust_method = 'GMM',
  dim_reduction_method = 'ICA',
  reduction_dim = 1,
  loop = 1,
  clus_drop_size = 5,
  approx_absorption_sample_size = 0,
  abs_pval_rejection_thresh = 1e-2,
)
{
  if(missing(data)){
  info = load_data(base_dir=data_dir)
  data = info$metadata
  annotation = info$annot
  phenoData = info$phenoData
  }
  
  GRset=  preprocess(data)
  mval = getM(GRset)
  
  train_inds = dim(mval)[1]
  
  mval_train = mval[,train_inds]
  mval_train = additional_preprocess(mval_train, annotation)
  
  BER = batch_adjust(mval_train, phenoData)
  phenoData = BER$phenoData
  mval_train_BER = BER$mval_train_BER
  
  DATA = t(mval_train_BER)
  
  DATAt = iqr_filter(DATA)
  DATAt_2 = scale_data(DATAt)
  
  for(reduce_count in 1:loop)
  {
    tree2 = kmeansvar(DATAt_2, init=num_var_clusts)
    processed_info = absorb_filter_reduce(DATAt_2, tree2, num_clusts=2, num_var_clusts=num_var_clustss, clus_drop_size=clus_drop_size, 
                                    approx_absorption_sample_size=approx_absorption_sample_size, 
                                    pval_rejection_thresh=abs_pval_rejection_thresh)
    X = processed_info$signals
    DATAt_2 = processed_info$clusts_filtered_daa
  }
  
  if(final_clust_method == 'GMM')
  {
    final_clusts = Mclust(X, num_clusts)$classification
  }
  else if(final_clust_method == 'kmeans')
  {
    final_clusts = kmeans(X, num_clusts)$cluster
  }
  else{ print('not a valid clustering method.')}
  return(final_clusts)
}
  
load_data <- function(base_dir) 
{
  targets = read.metharray.sheet(baseDir)
  RGSet <- read.metharray.exp(targets = targets)
  annotation <- getAnnotation(RGSet)
  phenoData <- pData(RGSet)
  return(list(metadata=RGSet, annot=annotation, phenoData=phenoData))
}

preprocess = function(data ,method='raw')
{
  switch(method,
         raw={return(preprocessRaw(data))}, {print("not valid preprocessing method")})
}


additional_preprocess = function(data, annotation, preprocess_funcs=list())
{
  if(length(preprocess_funcs) == 0){
    mt = mval_train[-which(rownames(data) %in%
                                     annotation$Name[which(annotation$chr == "chrY" | annotation$chr == "chrX")]),]
    RemInfInd = which(is.finite(rowSums(mt)))
    mt <- mt[RemInfInd,]
    return(mt)
  }
}


batch_adjust = function(data, phenoData)
{
  batch = phenoData$Batch
  
  if (length(which(is.na(batch)))>0){
    data <- data[,-which(is.na(batch))]
    phenoData = phenoData[-which(is.na(batch)),]
    batch = phenoData$Batch
  }
  
  modcombat = model.matrix(~ 1, data = phenoData)
  mval_train_BER = ComBat(
    dat = data,
    batch = batch,
    mod = modcombat,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  
  return(list(mval_train_BER=mval_train_BER, phenoData=phenoData))
}


iqr_filter = function(data, IQR_filtered_thresh)
{
  iqrvals = colIQRs(data)
  datat = data[,which(iqrvals > IQR_filtered_thresh)]
  return(datat)
}

scale_data = function(data, annotation=NA)
{
  #datat = data[,order(annotation$pos[which(annotation$Name %in% colnames(DATAt))])]
  datat = scale(data)
  return(datat)
}


absorb_filter_reduce = function(data, tree, num_clusts=2, num_var_clusts, 
                                clus_drop_size=5, approx_absorption_sample_size=0, 
                                pval_rejection_thresh=1e-2, absorption_thresh=1.0, dim_reduction_method = 'ICA',
                                reduction_dim=1)
{
  # num_clust > 2 is not supported.
  TOTAL_NUM = dim(data)[1]
  clus_clus = matrix(0, TOTAL_NUM, num_var_clusts)
  for(i in 1:num_var_clusts)
  {
    MDL = Mclust(data[, which(colnames(data) %in% names(tree$cluster[tree$cluster == i]))], num_clusts, verbose= FALSE)
    clus_clus[,i] = MDL$classification
  }
  
  if (approx_absorption_sample_size == 0)
  {
    approx_absorption_sample_size = dim(data)[2]
  }
  
  vals = apply(DATA[,sample(1:dim(DATA)[2], approx_absorption_sample_size)], 2, function(col){
    f = apply(clus_clus, 2, function(clcl_col){
      if(sum(clcl_col == 1) >= (TOTAL_NUM - clus_drop_size) | sum(clcl_col == 1) <= clus_drop_size)
      {
        10 # don't panic. This is just an invalid value
      }
      else{
        ks.test(col[which(clcl_col == 1)], col[which(clcl_col == 2)])$p.value
        # other tests may be used here. Based on our assumptions on DNAm data on the detection of rare diseases we used ks test.
      }
    })
    ret_val = c(min(f), which(f == min(f))[1])
  })
  
  goodval_inds = vals[2,which(vals[1,] < pval_rejection_thresh)]
  absorb = rep(0, num_var_clusts)
  for(i in 1:num_var_clusts)
  {
    absorb[i] = length(which(vals[2,goodval_inds] == i))
  }
  # absorb = as.numeric(table(vals[2,]))
  
  X = c()
  selected_clusts = c()
  for(i in 1:num_var_clust)
  {
    if(absorb[i] >= 0 & absorb[i] < absorption_thresh)
    {
      if(((length(which(clus_clus[,i] == 1)) <= clus_drop_size) | (length(which(clus_clus[,i] == 1)) >= TOTAL_NUM - clus_drop_size)))
      {
        next
      }
      if(dim_reduction_method == 'ICA'){
        X = cbind(X, fastICA(data[, which(colnames(data) %in% names(tree2$cluster[tree2$cluster == i]))], 2)$S[,1:reduction_dim])
      }
      else if (dim_reduction_method == 'PCA')
      {
        X = cbind(X, prcomp(data[, which(colnames(data) %in% names(tree2$cluster[tree2$cluster == i]))])$x[,1:reduction_dim]) 
      }
      else if (is.na(dim_reduction_method))
      {
        X = cbind(X, data[, which(colnames(data) %in% names(tree2$cluster[tree2$cluster == i]))])
      }
      selected_clusts = c(selected_clusts, i)
    }
  }
  clusts_filtered_data = data[, which(colnames(data) %in% names(tree2$cluster[tree2$cluster %in% selected_clusts]))]
  return(list(signals=X, clusts_filtered_data=clusts_filtered_data))
}






