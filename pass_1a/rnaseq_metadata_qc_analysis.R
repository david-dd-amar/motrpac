setwd("/Users/David/Desktop/MoTrPAC/data/pass_1a/rnaseq/")
library(corrplot)

metadata_files = list(
  stanford = "./stanford/metadata_PASS1A_RNA_Batch1.csv",
  sinai = "./sinai/batch1_meta_info.txt"
)
qc_files = list(
  stanford = "./stanford/Stanford_PASS1A_rnaseq_qc_info_bic.csv",
  sinai = "./sinai/Sinai_batch1_qc_info_bic.csv"
)

cols_for_qc_analysis = c("RIN","reads","pct_GC","pct_rRNA","pct_globin",
                      "pct_umi_dup","median_5_3_bias","pct_trimmed_bases",
                      "pct_multimapped")

metadata = list()
for(site in names(metadata_files)){
  currfile = metadata_files[[site]]
  currmeta = NULL
  try({currmeta = read.csv(currfile,stringsAsFactors = F,header=T)})
  if(is.null(currmeta)){
    try({currmeta = read.delim(currfile,stringsAsFactors = F,header=T)})
  }
  if(is.null(currmeta)){
    try({currmeta = read.table(currfile,stringsAsFactors = F,header=T,sep=" ")})
  }
  if(is.null(currmeta)){next}
  rownames(currmeta) = currmeta[,1]
  metadata[[site]] = currmeta
}
sapply(metadata,dim)

qc_info = list()
for(site in names(metadata_files)){
  currfile = qc_files[[site]]
  currqc = NULL
  try({currqc = read.csv(currfile,stringsAsFactors = F,header=T)})
  if(is.null(currmeta)){next}
  rownames(currqc) = currqc[,1]
  qc_info[[site]] = currqc
}
sapply(qc_info,dim)

# QC - should not print anyting
for(site in names(qc_info)){
  currmeta = metadata[[site]]
  currqc = qc_info[[site]]
  if(length(intersect(rownames(currqc),rownames(currmeta))) != nrow(currqc)){
    print(paste(site," qc and meta tables do not have the same row names"))
  }
}
meta_colnames = colnames(metadata[[1]])
qc_colnames = colnames(qc_info[[1]])
for(j in 2:length(qc_info)){
  if(!all(meta_colnames==colnames(metadata[[j]]))){
    print("Site metadata do not have the same columns")
  }
  if(!all(qc_colnames==colnames(qc_info[[j]]))){
    print("Site qc tables do not have the same columns")
  }
}

# Assuming QC above did not fail - put all metadata in a single file
rnaseq_meta = c()
site_metadata_with_qc = list()
for(site in names(qc_info)){
  currmeta = metadata[[site]]
  currmeta = cbind(currmeta,rep(site,nrow(currmeta)))
  colnames(currmeta)[ncol(currmeta)] = "site"
  currqc = qc_info[[site]]
  currqc = currqc[rownames(currmeta),]
  site_metadata_with_qc[[site]] = cbind(currmeta,currqc)
  rnaseq_meta = rbind(rnaseq_meta,site_metadata_with_qc[[site]])
}
dim(rnaseq_meta)

# Look at the correlation between some selected columns
x = site_metadata_with_qc$stanford[,cols_for_qc_analysis]
corrplot(cor(x,method="spearman"),tl.cex = 0.8,type = "upper")

# Implementation of the MOP's section 9 flagging of problematic samples
get_problematic_samples<-function(qc_and_meta_data){
  pct_cols = colnames(qc_and_meta_data)[grepl("pct_",colnames(qc_and_meta_data))]
  pct_mapped_cols = pct_cols[grepl("mapped",pct_cols)]
  pct_unmapped_cols = pct_cols[grepl("unmapped",pct_cols)]
  total_unmapped = rowSums(qc_and_meta_data[,pct_unmapped_cols])
  total_exonic = qc_and_meta_data$pct_coding + qc_and_meta_data$pct_utr
  # read_cols = colnames(qc_and_meta_data)[grepl("read",colnames(qc_and_meta_data))]
  is_problematic = qc_and_meta_data$RIN < 6 |
    qc_and_meta_data$reads < 20000000 |
    qc_and_meta_data$pct_GC > 80 |
    qc_and_meta_data$pct_GC < 20 |
    qc_and_meta_data$pct_rRNA > 20 |
    qc_and_meta_data$pct_uniquely_mapped < 60 | 
    total_exonic < 50
  return(is_problematic)
}
rnaseq_meta = cbind(rnaseq_meta,get_problematic_samples(rnaseq_meta))
colnames(rnaseq_meta)[ncol(rnaseq_meta)] = "IsFlagged"
table(rnaseq_meta$IsFlagged)
save(rnaseq_meta,file="rnaseq_meta.RData")

