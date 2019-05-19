setwd("/Users/David/Desktop/MoTrPAC/data/pass_1a/rnaseq/")
library(data.table)
library(DESeq2)
library(preprocessCore)
library(ggplot2)

# Data paths
site2fpkm_path = list(
  stanford = "./stanford/rsem_genes_fpkm_pass1a_batch1_Stanford.csv",
  sinai = "./sinai/rsem_genes_fpkm_pass1a_batch1_Sinai.csv"
)
site2genecount_path = list(
  stanford = "./stanford/rsem_genes_count_pass1a_batch1_Stanford.csv",
  sinai = "./sinai/rsem_genes_count_pass1a_batch1_Sinai.csv"
)

# load the metadata
# this is a data frame called rnaseq_meta that contains 
# the qc and sample metadata from both sites 
load("./rnaseq_meta.RData")
rnaseq_meta$Tissue = tolower(rnaseq_meta$Tissue)
rnaseq_meta$Tissue = gsub(" powder","",rnaseq_meta$Tissue)
rnaseq_meta[rnaseq_meta=="N/A"] = NA

# read the data in
site2fpkm = list()
site2counts = list()
for(site in names(site2fpkm_path)){
  currfpkm = fread(site2fpkm_path[[site]],header = T,
                   stringsAsFactors = F,data.table = F)
  rownames(currfpkm) = currfpkm[,1]
  currfpkm = currfpkm[,-1]
  site2fpkm[[site]] = currfpkm
  
  currcounts = fread(site2genecount_path[[site]],
                     header = T,stringsAsFactors = F,data.table = F)
  rownames(currcounts) = currcounts[,1]
  currcounts = currcounts[,-1]
  site2counts[[site]] = currcounts
}
###########################################################
###########################################################
###########################################################
# Pipeline 1: work with FPKMs
# Pipeline 1: work with FPKMs
#' Takes an FPKM matrix, removes lowly expressed genes and log transform
#' the remaining matrix
#' @return A matrix of log transformed FPKMs
process_fpkm1 <-function(fpkm_matrix){
  lowly_expressed_genes = rowSums(fpkm_matrix==0)/ncol(fpkm_matrix) > 0.2
  fpkm_matrix = fpkm_matrix[!lowly_expressed_genes,]
  fpkm_matrix = log(fpkm_matrix+1,base = 2)
  return(fpkm_matrix)
}
run_quantile_normalization<-function(x){
  x = as.matrix(x)
  mode(x) = "numeric"
  newx = preprocessCore::normalize.quantiles.robust(x)
  rownames(newx) = rownames(x)
  colnames(newx) = colnames(x)
  return(newx)
}

site_proc_fpkms = lapply(site2fpkm,process_fpkm1)
sapply(site_proc_fpkms,dim)
shared_genes = intersect(rownames(site_proc_fpkms[[1]]),rownames(site_proc_fpkms[[2]]))
proc_fpkms = cbind(site_proc_fpkms[[1]][shared_genes,],site_proc_fpkms[[2]][shared_genes,])
all(colnames(proc_fpkms) %in% rownames(rnaseq_meta))

# Do the PCA
fpkm_pca = prcomp(t(proc_fpkms),retx = T)
fpkm_pcax = fpkm_pca$x
# Without filtering the FPKMs
fpkm_all = cbind(site2fpkm[[1]],site2fpkm[[2]])
fpkm_all = log(fpkm_all+1,base=2)
fpkm_pca = prcomp(t(fpkm_all),retx = T)
fpkm_pcax = fpkm_pca$x

df = data.frame(fpkm_pcax[,1:10],
                tissue = rnaseq_meta[rownames(fpkm_pcax),"Tissue"],
                site = rnaseq_meta[rownames(fpkm_pcax),"site"])
ggplot(df,aes(x=PC1, y=PC2,shape=site, color=tissue)) + 
  geom_point(size=2) + ggtitle("PCA using FPKMs: PCs 1 and 2") + 
  theme(plot.title = element_text(hjust = 0.5))
ggplot(df,aes(x=PC4, y=PC3,shape=site, color=tissue)) + 
  geom_point(size=2) + ggtitle("PCA using FPKMs: PCs 3 and 4") + 
  theme(plot.title = element_text(hjust = 0.5))


# library(plotly)
# currcols = get_cols_vector_from_names(rnaseq_meta[rownames(fpkm_pcax),"site"],rainbow)
# curr_pchs = get_pch_vector_from_names(rnaseq_meta[rownames(fpkm_pcax),"site"])
# p <-plot_ly(x=fpkm_pcax[,1],y=fpkm_pcax[,2],color=names(currcols[[1]]),
#             ylab = "PC2",xlab="PC1",symbol = names(curr_pchs[[1]]),
#             cex=1.3,main="RNA-seq PCA by site",marker = list(size = 10))
# p
# 
# currcols = get_cols_vector_from_names(rnaseq_meta[rownames(fpkm_pcax),"Tissue"],rainbow)
# curr_pchs = get_pch_vector_from_names(rnaseq_meta[rownames(fpkm_pcax),"Tissue"])
# p <-plot_ly(x=fpkm_pcax[,1],y=fpkm_pcax[,2],color=names(currcols[[1]]),
#             ylab = "PC2",xlab="PC1",symbol = names(curr_pchs[[1]]),
#             cex=1.3,main="RNA-seq PCA by site",marker = list(size = 5))
# p

###########################################################
###########################################################
###########################################################
# Pipeline 2: work with count data
count_matrix = as.matrix(cbind(site2counts[[1]],site2counts[[2]]))
process_counts<-function(count_matrix){
  mode(count_matrix) = "integer"
  se <- SummarizedExperiment(count_matrix)
  dds <- DESeqDataSet(se, design = ~ 1 )
  #Estimate size factors
  dds <- estimateSizeFactors( dds )
  # Plot the size factors
  plot(sizeFactors(dds), colSums(counts(dds)))
  abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
  dds <- estimateDispersions(dds)
  return(dds)
}
dds = process_counts(count_matrix)
#The argument normalized equals true, divides each column by its size factor.
logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
pc <- prcomp( t( logcounts ) )
counts_pcax = pc$x
# VSD 
vsd <- varianceStabilizingTransformation(dds)
boxplot(assay(vsd)[,1:10])
pc2 <- prcomp( t( assay(vsd) ) )
counts_pcax = pc2$x

df = data.frame(counts_pcax[,1:10],
                tissue = rnaseq_meta[rownames(counts_pcax),"Tissue"],
                site = rnaseq_meta[rownames(counts_pcax),"site"])
ggplot(df,aes(x=PC1, y=PC2,shape=site, color=tissue)) + 
  geom_point(size=2) + ggtitle("PCA using normalized counts (DESeq2)") + 
  theme(plot.title = element_text(hjust = 0.5))
ggplot(df,aes(x=PC3, y=PC4,shape=site, color=tissue)) + 
  geom_point(size=2) + ggtitle("PCA using normalized counts (DESeq2)") + 
  theme(plot.title = element_text(hjust = 0.5))

###########################################################
###########################################################
###########################################################
# Site comparison: gastro tissue using FPKM data
tissue = rnaseq_meta$Tissue
names(tissue) = rownames(rnaseq_meta)
gastro_fpkm = lapply(site2fpkm,
    function(x,y)x[,grepl("gastro",y[colnames(x)],ignore.case = T)],y=tissue)
sapply(gastro_fpkm,dim)
corrs = cor(gastro_fpkm[[1]],gastro_fpkm[[2]],method="spearman")
colMeans(corrs)
rowMeans(corrs)
gastro_fpkm_processed = lapply(gastro_fpkm,process_fpkm1)
sapply(gastro_fpkm_processed,dim)
shared_genes = intersect(rownames(gastro_fpkm_processed[[1]]),
                         rownames(gastro_fpkm_processed[[2]]))
gastro_fpkm_mat = cbind(gastro_fpkm_processed[[1]][shared_genes,],
                        gastro_fpkm_processed[[2]][shared_genes,])

# Analysis of the metadata
gastro_metadata = rnaseq_meta[colnames(gastro_fpkm_mat),]
sample_id = paste(gastro_metadata$BID,gastro_metadata$PID,sep=";")
names(sample_id) = rownames(gastro_metadata)
table(sample_id)
gastro_metadata = gastro_metadata[order(gastro_metadata$site,sample_id),]
sample_id = sample_id[rownames(gastro_metadata)]

metadata2site_pval = c()
site_ind = gastro_metadata$site==gastro_metadata$site[1]
for(col in names(gastro_metadata)){
  x = gastro_metadata[[col]]
  if(! mode(x)=="numeric"){next}
  x1 = x[site_ind];x2 = x[!site_ind]
  try({metadata2site_pval[col] = wilcox.test(x1,x2,paired=T)$p.value})
}
selected_qc_comparisons = sort(metadata2site_pval)[1:30]
selected_qc_comparisons[grepl("pct_",names(selected_qc_comparisons))]

# Comparison 1: qc scores
par(mfrow=c(2,2))
boxplot(pct_umi_dup~site,data=gastro_metadata,main="% UMI dup")
boxplot(r_260_230~site,data=gastro_metadata,main = "260/230");abline(h = 1.8,lty=2,col="red")
boxplot(pct_rRNA~site,data=gastro_metadata,main="%rRNA")
boxplot(pct_globin~site,data=gastro_metadata,main="%Globin")

# Comparison 2: boxplots
currcols = get_cols_vector_from_names(rnaseq_meta[colnames(gastro_fpkm_mat),"site"],rainbow)
curr_pchs = get_pch_vector_from_names(rnaseq_meta[colnames(gastro_fpkm_mat),"site"])
inds_for_boxplot = c(1:10,81:90)
boxplot(gastro_fpkm_mat[,inds_for_boxplot],col=currcols[[1]][inds_for_boxplot],
        ylim = c(0,20))
legend(x="top",names(currcols[[2]]),fill=currcols[[2]])

# Comparison 3: PCA
gastro_fpkm_pca = prcomp(t(gastro_fpkm_mat))
gastro_fpkm_pcax = gastro_fpkm_pca$x
plot(x=gastro_fpkm_pcax[,1],y=gastro_fpkm_pcax[,2],col=currcols[[1]],
     ylab = "PC2",xlab="PC1",pch = curr_pchs[[1]],
     main="PCA, normalized counts",cex=0.7,lwd=2)
legend(x="bottomright",names(currcols[[2]]),fill=currcols[[2]])
# Retry - quantile normalization before PCA
gastro_fpkm_mat_q = run_quantile_normalization(gastro_fpkm_mat)
gastro_fpkm_pca = prcomp(t(gastro_fpkm_mat_q))
gastro_fpkm_pcax = gastro_fpkm_pca$x
df = data.frame(gastro_fpkm_pcax[,1:10],
                shape=gastro_metadata$site, color=gastro_metadata$site,
                sample = sample_id[rownames(gastro_fpkm_pcax)])
df = df[order(df$sample),]
ggplot(df,aes(x=PC1, y=PC2, shape=shape, color=color,group=sample)) +
  geom_point(size=2) + geom_path(size=0.1,color="blue")

# Comparison 4: pairwise correlation 
# Raw FPKMs
x1 = site2fpkm[[1]][,colnames(gastro_fpkm_processed[[1]])]
x2 = site2fpkm[[2]][,colnames(gastro_fpkm_processed[[2]])]
x1 = x1[,order(sample_id[colnames(x1)])]
x2 = x2[,order(sample_id[colnames(x2)])]
all(sample_id[colnames(x1)] == sample_id[colnames(x2)])
corrs = cor(x1[shared_genes,],x2[shared_genes,],method="spearman")
# corrs = cor(x1,x2,method="spearman")
l = list(
  within = diag(corrs),
  between = corrs[lower.tri(corrs,diag = F)]
)
boxplot(l)

# Processed FPKMs
x1 = gastro_fpkm_processed[[1]]
x2 = gastro_fpkm_processed[[2]]
x1 = x1[,order(sample_id[colnames(x1)])]
x2 = x2[,order(sample_id[colnames(x2)])]
all(sample_id[colnames(x1)] == sample_id[colnames(x2)])
corrs = cor(x1[shared_genes,],x2[shared_genes,],method="spearman")
# corrs = cor(x1,x2,method="spearman")
l = list(
  within = diag(corrs),
  between = corrs[lower.tri(corrs,diag = F)]
)
boxplot(l)


#' #'''
#' #' @assumption: first column in each results matrix is the row name
#' get_matrix_from_results_list<-function(l,name){
#'   m = as.matrix(l[[1]][,name],ncol=1)
#'   rownames(m) = l[[1]][,1]
#'   colnames(m)[1] = names(l)[1]
#'   if(length(l)<2){return(m)}
#'   for(j in 2:length(l)){
#'     if(! all(l[[j]][,1]==rownames(m))){
#'       stop(paste("results row names do not match at item number",j))
#'     }
#'     m = cbind(m,l[[j]][,name])
#'     colnames(m)[ncol(m)] = names(l)[[j]]
#'   }
#'   return(m)
#' }
#' # Read RSEM results, get FPKMs and analyze the matrix
#' rsem_results = list()
#' for(ff in rsem_files){
#'   curr_res = fread(paste(rsem_path,ff,sep=""),data.table = F,stringsAsFactors = F)
#'   curr_name = strsplit(ff,split="\\.")[[1]][1]
#'   rsem_results[[curr_name]] = curr_res
#' }
#' sapply(rsem_results,dim)
# fpkm_matrix = get_matrix_from_results_list(rsem_results,"FPKM")


