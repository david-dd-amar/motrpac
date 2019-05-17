setwd("/Users/David/Desktop/MoTrPAC/data/pass_1a/rnaseq/")
library(data.table)
library(DESeq2)
library(preprocessCore)

rsem_path = "./stanford/gene_rsem/"
rsem_files = list.files(rsem_path)
load("./stanford_qc_and_meta.RData")

#'''
#' @assumption: first column in each results matrix is the row name
get_matrix_from_results_list<-function(l,name){
  m = as.matrix(l[[1]][,name],ncol=1)
  rownames(m) = l[[1]][,1]
  colnames(m)[1] = names(l)[1]
  if(length(l)<2){return(m)}
  for(j in 2:length(l)){
    if(! all(l[[j]][,1]==rownames(m))){
      stop(paste("results row names do not match at item number",j))
    }
    m = cbind(m,l[[j]][,name])
    colnames(m)[ncol(m)] = names(l)[[j]]
  }
  return(m)
}

# Read RSEM results, get FPKMs and analyze the matrix
rsem_results = list()
for(ff in rsem_files){
  curr_res = fread(paste(rsem_path,ff,sep=""),data.table = F,stringsAsFactors = F)
  curr_name = strsplit(ff,split="\\.")[[1]][1]
  rsem_results[[curr_name]] = curr_res
}
sapply(rsem_results,dim)

fpkm_matrix = get_matrix_from_results_list(rsem_results,"FPKM")
lowly_expressed_genes = rowSums(fpkm_matrix==0)/ncol(fpkm_matrix) > 0.2
fpkm_matrix = fpkm_matrix[!lowly_expressed_genes,]
fpkm_matrix = log(fpkm_matrix+1,base = 2)

stanford_qc_and_meta = stanford_qc_and_meta[colnames(fpkm_matrix),]
length(intersect(colnames(fpkm_matrix),rownames(stanford_qc_and_meta))) == nrow(stanford_qc_and_meta)
tissue = as.character(stanford_qc_and_meta$Tissue)
tissue_col = rainbow(length(unique(tissue)))
names(tissue_col) = unique(tissue)
is_flagged_pch = 20+as.numeric(stanford_qc_and_meta$IsFlagged)*5

boxplot(fpkm_matrix[,1:50],col=tissue_col[tissue[1:50]],las=2,cex.axis=0.5,
        ylab="log2FPKM")
legend(x="top",names(tissue_col),fill=tissue_col)
boxplot(fpkm_matrix[,tissue==tissue[10]],las=2,cex.axis=0.5,ylab="log2FPKM")

fpkm_pca = prcomp(t(fpkm_matrix),retx = T)
fpkm_pcax = fpkm_pca$x
plot(x=fpkm_pcax[,1],y=fpkm_pcax[,2],col=tissue_col[tissue],
     ylab = "PC2",xlab="PC1",pch=is_flagged_pch,
     cex=1.3,main="RNA-seq Stanford",
     lwd=1+2*as.numeric(stanford_qc_and_meta$IsFlagged))
legend(x="topleft",names(tissue_col),fill=tissue_col)

library(plotly)

