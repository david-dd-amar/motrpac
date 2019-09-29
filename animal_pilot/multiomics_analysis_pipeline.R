setwd("/Users/David/Desktop/MoTrPAC/march_2019/pilot_data/")

# Read in the RNAseq data
rnaseq_files = list.files("rnaseq/")
rnaseq_files = rnaseq_files[grepl("genes",rnaseq_files)]
rnaseq_files = rnaseq_files[grepl("M_",rnaseq_files) | grepl("A_",rnaseq_files)]
rnaseq_data = c()
rnaseq_data_counts = c(); rnaseq_data_gene_length=c()
for(ff in rnaseq_files){
  currd = read.table(paste("rnaseq/",ff,sep=""),header=T,row.names = 1)
  if(length(rnaseq_data)==0){
    rnaseq_data = matrix(currd$FPKM,ncol=1)
    rownames(rnaseq_data) = rownames(currd)
    rnaseq_data_counts = matrix(currd$expected_count,ncol=1)
    rownames(rnaseq_data_counts) = rownames(rnaseq_data)
    rnaseq_data_gene_length = matrix(currd$effective_length,ncol=1)
    rownames(rnaseq_data_gene_length) = rownames(rnaseq_data)
  }
  else{
    if(!all(rownames(currd)==rownames(rnaseq_data))){
      print("Error, row names do not match")
      print(ff)
    }
    rnaseq_data = cbind(rnaseq_data,currd$FPKM)
    rnaseq_data_counts = cbind(rnaseq_data_counts,currd$expected_count)
    rnaseq_data_gene_length = cbind(rnaseq_data_gene_length,currd$effective_length)
  }
  colnames(rnaseq_data)[ncol(rnaseq_data)] = ff
  colnames(rnaseq_data_counts)[ncol(rnaseq_data_counts)] = ff
  colnames(rnaseq_data_gene_length)[ncol(rnaseq_data_gene_length)] = ff
}
dim(rnaseq_data)
colnames(rnaseq_data) = gsub("genes.results","",colnames(rnaseq_data))
to_rem = apply(rnaseq_data<=1,1,all)
table(to_rem)
rnaseq_data = rnaseq_data[!to_rem,]
rnaseq_data = log(1+rnaseq_data,base=2)
rnaseq_data_counts = rnaseq_data_counts[rownames(rnaseq_data),]
rnaseq_data_gene_length = rnaseq_data_gene_length[rownames(rnaseq_data),]
par(mar=c(8,4,2,2))
boxplot(rnaseq_data,las=2)

# helper functions for data formatting
get_scores_vec<-function(sample,dframe,col_sample,col_molecule_names,col_scores){
  m = dframe[as.character(dframe[,col_sample])==sample,c(col_molecule_names,col_scores)]
  v = m[,2]
  names(v) = as.character(m[,1])
  return(v)
}
scores_list_to_matrix<-function(l){
  names_intersect = names(l[[1]])
  for(nn in names(l)){
    print(length(names_intersect))
    names_intersect = intersect(names_intersect,names(l[[nn]]))
  }
  m = matrix(NA,nrow=length(names_intersect),ncol=length(l),
             dimnames = list(names_intersect,names(l)))
  for(samp in colnames(m)){
    m[,samp] = l[[samp]][names_intersect]
  }
  return(m)
}

#'
#' Create a list of dataset by site.
#' For each site the data are:
#' 1. A matrix of the protein levels (intensity or ph) for each sample (named)
#' 2. sample2plex (batch)
#' 3. sample2condition
format_prot_data<-function(d,col_sample = "BioReplicate",col_batch="Run",
                           col_cond = "Condition",col_site="Site",
                           col_molecule_names="Protein",col_scores="Abundance"){
  l = list()
  sites = unique(as.character(d[,col_site]))
  for(si in sites){
    currd = d[d[,col_site]==si,]
    tmp = unique(currd[,c(col_batch,col_sample)])
    sample2batch = as.character(tmp[,1])
    names(sample2batch) = as.character(tmp[,2])
    tmp = unique(currd[,c(col_cond,col_sample)])
    sample2condition = as.character(tmp[,1])
    names(sample2condition) = as.character(tmp[,2])
    curr_l = sapply(names(sample2condition),get_scores_vec,dframe=currd,
                    col_sample=col_sample,col_molecule_names=col_molecule_names,
                    col_scores=col_scores)
    if(class(curr_l)!="matrix"){
      curr_m = scores_list_to_matrix(curr_l)
    }
    else{
      curr_m = curr_l
    }
    l[[si]] = list(sample2batch=sample2batch,
                   sample2condition=sample2condition,m=curr_m)
  }
  return(l)
}

# Proteomics - intensity
load("proteomics/pilot.cas.prot.global_oneid_intensity_norm.Rdata")
prot_global_intensity = format_prot_data(pilot.cas.prot.global_oneid_intensity_norm)
boxplot(prot_global_intensity[[1]]$m)
# boxplot(run_quantile_norm(prot_global_intensity[[1]]$m))

# Proteomics - phospho
load("proteomics/pilot.cas.prot.phglobal_oneid_intensity_norm.Rdata")
prot_global_ph = format_prot_data(pilot.cas.prot.phglobal_oneid_intensity_norm)
boxplot(prot_global_ph[[1]]$m)

# Metabolomics
load("metabolomics/pilot.duke.metab.targeted.quantifications.Rdata")
met_duke = format_prot_data(pilot.duke.metab.targeted.quantifications.long,
                            col_sample = "BioReplicate",col_batch="Animal",
                            col_cond = "Condition",col_site="Tissue",
                            col_molecule_names="Metabolite",
                            col_scores="Calculated.Concentration")
boxplot(log(1+met_duke[[1]]$m,base=2))
load("metabolomics/pilot.metab.emory.target.oxylipids.long.Rdata")
met_emory = format_prot_data(pilot.metab.emory.target.oxylipids.long,
                            col_sample = "BioReplicate",col_batch="Sample",
                            col_cond = "Condition",col_site="Tissue",
                            col_molecule_names="Component.Name",
                            col_scores="Calculated.Concentration")
boxplot(log(1+met_emory[[1]]$m,base=2))
dim(met_emory$Adipose$m)
load("metabolomics/pilot.metab.umichigan.untarget.long.Rdata")
x1 = pilot.metab.umichigan.untarget.muscle.positive.quant.long
x2 = pilot.metab.umichigan.untarget.plasma.positive.quant.long
x1$Tissue = rep("Muscle",nrow(x1))
x2$Tissue = rep("Plasma",nrow(x2))
pilot.metab.umichigan.untarget = rbind(x1,x2)
met_umich = format_prot_data(pilot.metab.umichigan.untarget,
                               col_sample = "BioReplicate",col_batch="Animal",
                               col_cond = "Condition",col_site="Tissue",
                               col_molecule_names="HMDB_name",
                               col_scores="Intensity")
boxplot(log(1+met_umich[[1]]$m,base=2))
boxplot(log(1+met_umich[[2]]$m,base=2))

# Prepare datasets for integrative analysis

# Helper functions and useful packages
source("https://bioconductor.org/biocLite.R")
biocLite("org.Rn.eg.db")
library("preprocessCore");library("org.Rn.eg.db")
library(corrplot)
# Get some rat data
# Map to gene names
refseq2gene = as.list(org.Rn.egREFSEQ2EG)
hist(sapply(refseq2gene,length))
ensembl2gene = as.list(org.Rn.egENSEMBL2EG)
# just for this analysis, we take the first mapping
# when there is more than one
ensembl2gene = sapply(ensembl2gene,function(x)x[[1]])

# 0. Make sure all sample names are comparable
# rnaseq
rnaseq_sample_names = colnames(rnaseq_data)
rnaseq_sample_names = unname(sapply(rnaseq_sample_names,function(x)
  paste(strsplit(x,split="_")[[1]][1:2],collapse="_")))
colnames(rnaseq_data) = rnaseq_sample_names
colnames(rnaseq_data_counts) = rnaseq_sample_names
colnames(rnaseq_data_gene_length) = rnaseq_sample_names
# prot1
parse_prot_colnames<-function(x,tissue_prefix="M"){
  x = gsub("RUN.0","R",x)
  x = gsub("SED.0","S",x)
  x = gsub("RUN.1","R1",x)
  x = gsub("SED.1","S1",x)
  x = paste(tissue_prefix,x,sep="_")
  return(x)
}
for(i in 1:length(prot_global_intensity)){
  colnames(prot_global_intensity[[i]]$m) = parse_prot_colnames(
    colnames(prot_global_intensity[[i]]$m))
  print(colnames(prot_global_intensity[[i]]$m))
}
# prot2
for(i in 1:length(prot_global_intensity)){
  colnames(prot_global_ph[[i]]$m) = parse_prot_colnames(
    colnames(prot_global_ph[[i]]$m))
  print(colnames(prot_global_ph[[i]]$m))
}
# metabolomics
colnames(met_duke[[1]]$m) = parse_prot_colnames(colnames(met_duke[[1]]$m))
for(i in 1:length(met_emory)){
  pref = substr(names(met_emory)[i],start = 1,stop = 1)
  if(names(met_emory)[i]=="Standard"){perf="Standard"}
  colnames(met_emory[[i]]$m) = parse_prot_colnames(
    colnames(met_emory[[i]]$m),pref)
  print(colnames(met_emory[[i]]$m))
}
for(i in 1:length(met_umich)){
  pref = substr(names(met_umich)[i],start = 1,stop = 1)
  if(names(met_umich)[i]=="Standard"){perf="Standard"}
  colnames(met_umich[[i]]$m) = parse_prot_colnames(
    colnames(met_umich[[i]]$m),pref)
  print(colnames(met_umich[[i]]$m))
}

# Put all omics in one list. Entry name is omics;tissue;site
pilot_datasets = list()
pilot_datasets[["rnaseq;muscle;get"]] = rnaseq_data[,grepl("M_",colnames(rnaseq_data))]
pilot_datasets[["rnaseq;adipose;get"]] = rnaseq_data[,grepl("A_",colnames(rnaseq_data))]
pilot_datasets[["rnaseq_counts;muscle;get"]] = rnaseq_data_counts[,grepl("M_",colnames(rnaseq_data))]
pilot_datasets[["rnaseq_counts;adipose;get"]] = rnaseq_data_counts[,grepl("A_",colnames(rnaseq_data))]
pilot_datasets[["proteomics_abundance;muscle;pnnl"]] = prot_global_intensity$PNNL$m
pilot_datasets[["proteomics_abundance;muscle;broad"]] = prot_global_intensity$BROAD$m
pilot_datasets[["proteomics_ph;muscle;pnnl"]] = prot_global_ph$PNNL$m
pilot_datasets[["proteomics_ph;muscle;broad"]] = prot_global_ph$BROAD$m
pilot_datasets[["tar_met;muscle;duke"]] = met_duke$Muscle$m
pilot_datasets[["tar_met;muscle;emory"]] = met_emory$Muscle$m
pilot_datasets[["tar_met;adipose;emory"]] = met_emory$Adipose$m
pilot_datasets[["tar_met;plasma;emory"]] = met_emory$Plasma$m
pilot_datasets[["untar_met;muscle;umich"]] = met_umich$Muscle$m
pilot_datasets[["untar_met;plasma;umich"]] = met_umich$Plasma$m
metadata = list()
metadata[["proteomics_abundance;muscle;pnnl"]] = 
  list(batch=prot_global_intensity$PNNL$sample2batch)
metadata[["proteomics_abundance;muscle;broad"]] = 
  list(batch=prot_global_intensity$BROAD$sample2batch)
metadata[["proteomics_ph;muscle;pnnl"]] = 
  list(batch=prot_global_ph$PNNL$sample2batch)
metadata[["proteomics_ph;muscle;broad"]] = 
  list(batch=prot_global_ph$BROAD$sample2batch)
save(pilot_datasets,metadata,file="pilot_data_matrices.RData")

run_quantile_norm<-function(xx){
  qx = normalize.quantiles(xx)
  dimnames(qx) = dimnames(xx)
  return(qx)
}

# Some plots
xx = pilot_datasets$`rnaseq;muscle;get`
xx = pilot_datasets$`proteomics_ph;muscle;pnnl`
xx = xx[!apply(is.na(xx),1,any),]
qx = normalize.quantiles(xx)
dimnames(qx) = dimnames(xx)
qx = qx[apply(qx,1,sd)>0,]
conds = grepl("_R",colnames(xx))
boxplot(qx,las=2)
pcax = prcomp(t(qx),retx = T,center = T,scale. = T)
plot(pcax$x[,1],pcax$x[,2],col=as.numeric(conds)+10)
text(pcax$x[,1],pcax$x[,2],labels=rownames(pcax$x))

# 1. Differential abundance analysis

library(lme4);library(lmerTest)
run_stat_test<-function(x,y,ctrl=NULL,func=t.test,...){
  if(is.null(ctrl)){ctrl=y[1]}
  x1 = x[y==ctrl]
  x2 = x[y!=ctrl]
  tt = func(x1,x2,...)
  d = mean(x2)-mean(x1)
  return(c(statistic=tt$statistic,d=d,p=tt$p.value))
}
run_lmm_test<-function(x,y,batch,func=t.test,...){
  model = lmer(x~y + (1|batch))
  model_test = summary(model)
  return(model_test$coefficients[2,])
}

# Use a simple t-test by default
diff_results = list()
for(nn in names(pilot_datasets)){
  if(is.element(nn,set=names(diff_results))){next}
  x = pilot_datasets[[nn]]
  x = run_quantile_norm(x)
  x = x[apply(x,1,sd)>0,]
  x = x[!apply(is.na(x),1,any),]
  inds = grepl("_R",colnames(x))|grepl("_S",colnames(x))
  x = x[,inds]
  curr_y = rep("R",ncol(x))
  curr_y[grepl("_S",colnames(x))] = "S"
  diff_results[[nn]] = t(apply(x,1,run_stat_test,y=curr_y,ctrl="S"))
}

# When we have a batch, try LMMs
for(nn in names(metadata)){
  x = pilot_datasets[[nn]]
  x = run_quantile_norm(x)
  x = x[apply(x,1,sd)>0,]
  x = x[!apply(is.na(x),1,any),]
  inds = grepl("_R",colnames(x))|grepl("_S",colnames(x))
  x = x[,inds]
  y = rep("R",ncol(x))
  y[grepl("_S",colnames(x))] = "S"
  batch = metadata[[nn]]$batch
  names(batch) = colnames(pilot_datasets[[nn]])
  batch = batch[colnames(x)]
  diff_results[[nn]] = t(apply(x,1,run_lmm_test,y=y,batch=batch))
}

# For RNAseq with counts use DESeq2
library(DESeq2)
for(nn in names(pilot_datasets)){
  if(!grepl("rnaseq_count",nn)){next}
  x = pilot_datasets[[nn]]
  x = floor(x)
  y = rep("R",ncol(x))
  y[grepl("_S",colnames(x))] = "S"
  y = factor(y)
  dds <- DESeqDataSetFromMatrix(x, DataFrame(y), ~ y)
  # standard analysis
  dds <- DESeq(dds)
  res <- results(dds)
  res = res[,c(1:4,6:5)]
  diff_results[[nn]] = res
}

length(diff_results)
par(mfrow=c(4,3))
for(nn in names(diff_results)[!grepl("rnaseq;",names(diff_results))]){
  res = diff_results[[nn]]
  hist(res[,ncol(res)],main=nn,xlab="p-value")
  ps = res[,ncol(res)]
  print(paste(nn,sum(p.adjust(ps,method="fdr")<0.1,na.rm=T)))
}
dev.off()

# 2. Multiomics data

# 3. Site differences - proteomics

# Abundance data sample correlations
x1 = pilot_datasets$`proteomics_abundance;muscle;pnnl`
x1 = x1[!apply(is.na(x1),1,any),]
x2 = pilot_datasets$`proteomics_abundance;muscle;broad`

# Phospho data sample correlations
x1 = pilot_datasets$`proteomics_ph;muscle;pnnl`
x1 = x1[!apply(is.na(x1),1,any),]
x2 = pilot_datasets$`proteomics_ph;muscle;broad`

x1 = run_quantile_norm(x1)
x2 = run_quantile_norm(x2)

inds=intersect(rownames(x1),rownames(x2))
s_inds = intersect(colnames(x1),colnames(x2))
corrs = cor(x1[inds,s_inds],x2[inds,s_inds])
corrs = cor(x1[inds,s_inds])
corrplot(corrs)
diag(corrs)

# Compare proteomics datasets by their differential abundance
x1 = diff_results$`proteomics_abundance;muscle;broad`
x2 = diff_results$`proteomics_abundance;muscle;pnnl`
inds=intersect(rownames(x1),rownames(x2))
s_inds = intersect(colnames(x1),colnames(x2))
corrs = cor(x1[inds,s_inds],x2[inds,s_inds])
plot(y=x1[inds,1],x=x2[inds,1],ylab="LMM estimate, Broad",
     xlab="LMM estimate, PNNL",pch=20,cex=0.7,
     main=paste("Proteomics differential abundance, rho="
                ,format(corrs[1,1],digits=2),sep=""))
x1_bin = abs(x1[inds,1])>=sort(abs(x1[inds,1]),decreasing=T)[100]
x2_bin = abs(x2[inds,1])>=sort(abs(x2[inds,1]),decreasing=T)[100]
plot(table(x1_bin,x2_bin),main="Is in top 100?")
fisher.test(table(x1_bin,x2_bin))$p.value
cor.test(x1[inds,1],x2[inds,1],method = "spearman")

# Compare proteomics ph datasets by their differential abundance
x1 = diff_results$`proteomics_ph;muscle;pnnl`
x2 = diff_results$`proteomics_ph;muscle;broad`
inds=intersect(rownames(x1),rownames(x2))
s_inds = intersect(colnames(x1),colnames(x2))
corrs = cor(x1[inds,s_inds],x2[inds,s_inds])
plot(y=x1[inds,1],x=x2[inds,1],ylab="LMM estimate, Broad",
     xlab="LMM estimate, PNNL",pch=20,cex=0.7,
     main=paste("Proteomics ph differential abundance, rho="
                ,format(corrs[1,1],digits=2),sep=""))
x1_bin = abs(x1[inds,1])>sort(abs(x1[inds,1]),decreasing=T)[100]
x2_bin = abs(x2[inds,1])>sort(abs(x2[inds,1]),decreasing=T)[100]
plot(table(x1_bin,x2_bin),main="Is in top 100?")
fisher.test(table(x1_bin,x2_bin))
cor.test(x1[inds,1],x2[inds,1],method = "spearman")

# Compare with RNAseq
x1 = diff_results$`rnaseq_counts;muscle;get`
x1 = x1[!is.na(ensembl2gene[rownames(x1)]),]
rownames(x1) = ensembl2gene[rownames(x1)]
x2 = diff_results$`proteomics_abundance;muscle;pnnl`
rownames(x2) = sapply(rownames(x2),function(x)strsplit(x,split="\\.")[[1]][1])
rownames(x2) = refseq2gene[rownames(x2)]
inds=intersect(rownames(x1),rownames(x2))
corr = cor(x1[inds,ncol(x1)],x2[inds,ncol(x2)])
plot(y=x1[inds,2],x=x2[inds,1],ylab="RNAseq fchange",
     xlab="LMM estimate, PNNL",pch=20,cex=0.7,
     main=paste("Proteomics vs RNAseq, rho="
                ,format(corr,digits=2),sep=""))

x1 = diff_results$`rnaseq;muscle;get`
rownames(x1) = ensembl2gene[rownames(x1)]
x2 = diff_results$`proteomics_abundance;muscle;broad`
rownames(x2) = sapply(rownames(x2),function(x)strsplit(x,split="\\.")[[1]][1])
rownames(x2) = refseq2gene[rownames(x2)]
inds=intersect(rownames(x1),rownames(x2))
corr = cor(x1[inds,ncol(x1)],x2[inds,ncol(x2)])
plot(y=x1[inds,2],x=x2[inds,1],ylab="RNAseq fchange",
     xlab="LMM estimate, PNNL",pch=20,cex=0.7,
     main=paste("Proteomics vs RNAseq, rho="
                ,format(corr,digits=2),sep=""))

# 4. Site differences - metabolomics correlation analysis
cor(met_duke$Muscle$m,met_emory$Muscle$m)

# 5. Multiomics factorization
library(devtools)
source("https://bioconductor.org/biocLite.R")
biocLite("GenomeInfoDbData",dependencies=T)
devtools::install_github("bioFAM/MOFA", build_opts = c("--no-resave-data"))
source("https://bioconductor.org/biocLite.R")
biocLite("MultiAssayExperiment")

motrpac_multiomics_data1 = c()
# prot 1
x1 = pilot_datasets$`proteomics_abundance;muscle;broad`
x1 = x1[,grepl("_S",colnames(x1)) | grepl("_R",colnames(x1))]
rownames(x1) = sapply(rownames(x1),function(x)strsplit(x,split="\\.")[[1]][1])
rownames(x1) = refseq2gene[rownames(x1)]
rownames(x1) = paste("proteomics.",rownames(x1),sep="")
motrpac_multiomics_data1 = cbind(motrpac_multiomics_data1,t(x1))
# prot 2
x1 = pilot_datasets$`proteomics_ph;muscle;broad`
x1 = x1[,grepl("_S",colnames(x1)) | grepl("_R",colnames(x1))]
rownames(x1) = sapply(rownames(x1),function(x)strsplit(x,split="\\.")[[1]][1])
rownames(x1) = refseq2gene[rownames(x1)]
rownames(x1) = paste("proteomics_ph.",rownames(x1),sep="")
x1 = x1[,rownames(motrpac_multiomics_data1)]
motrpac_multiomics_data1 = cbind(motrpac_multiomics_data1,t(x1))
# untar met
x1 = pilot_datasets$`untar_met;muscle;umich`
x1 = x1[,grepl("_S",colnames(x1)) | grepl("_R",colnames(x1))]
rownames(x1)[is.na(rownames(x1))] = paste("met",1:sum(is.na(rownames(x1))),sep="")
rownames(x1) = paste("untar_met.",rownames(x1),sep="")
x1 = x1[,rownames(motrpac_multiomics_data1)]
motrpac_multiomics_data1 = cbind(motrpac_multiomics_data1,t(x1))
dim(motrpac_multiomics_data1)

# with rnaseq
motrpac_multiomics_data2 = c()
# rnaseq
x1 = pilot_datasets$`rnaseq;muscle;get`
x1 = x1[!is.na(ensembl2gene[rownames(x1)]),]
rownames(x1) = ensembl2gene[rownames(x1)]
rownames(x1) = paste("rnaseq.",rownames(x1),sep="")
motrpac_multiomics_data2 = cbind(motrpac_multiomics_data2,t(x1))
# prot 1
x1 = pilot_datasets$`proteomics_abundance;muscle;broad`
x1 = x1[,grepl("_S",colnames(x1)) | grepl("_R",colnames(x1))]
rownames(x1) = sapply(rownames(x1),function(x)strsplit(x,split="\\.")[[1]][1])
rownames(x1) = refseq2gene[rownames(x1)]
rownames(x1) = paste("proteomics.",rownames(x1),sep="")
x1 = x1[,rownames(motrpac_multiomics_data2)]
motrpac_multiomics_data2 = cbind(motrpac_multiomics_data2,t(x1))
# prot 2
x1 = pilot_datasets$`proteomics_ph;muscle;broad`
x1 = x1[,grepl("_S",colnames(x1)) | grepl("_R",colnames(x1))]
rownames(x1) = sapply(rownames(x1),function(x)strsplit(x,split="\\.")[[1]][1])
rownames(x1) = refseq2gene[rownames(x1)]
rownames(x1) = paste("proteomics_ph.",rownames(x1),sep="")
x1 = x1[,rownames(motrpac_multiomics_data2)]
motrpac_multiomics_data2 = cbind(motrpac_multiomics_data2,t(x1))
# untar met
x1 = pilot_datasets$`untar_met;muscle;umich`
x1 = x1[,grepl("_S",colnames(x1)) | grepl("_R",colnames(x1))]
rownames(x1)[is.na(rownames(x1))] = paste("met",1:sum(is.na(rownames(x1))),sep="")
rownames(x1) = paste("untar_met.",rownames(x1),sep="")
x1 = x1[,rownames(motrpac_multiomics_data2)]
motrpac_multiomics_data2 = cbind(motrpac_multiomics_data2,t(x1))
dim(motrpac_multiomics_data2)

save(motrpac_multiomics_data1,motrpac_multiomics_data2,file="motrpac_pilot_multiomics.RData")
getwd()

# Simple PCAs
run_pca_and_plot<-function(x,transpose=F,num_pcs=3,plot=T,...){
  if(transpose){x = t(x)}
  x = x[,!apply(is.na(x),2,any)]
  x = x[,!apply(x,2,sd)==0]
  x = scale(x)
  pcx = prcomp(x,retx=T)
  pcs = pcx$rotation[,1:num_pcs]
  eigvals = pcx$sdev
  pcs2_var = sum(eigvals[1:2])/sum(eigvals)
  pcx = pcx$x[,1:num_pcs]
  labs = rep("R",nrow(pcx))
  labs[grepl("_S",rownames(pcx))]="S"
  labs = as.factor(labs)
  cols = rep("blue",length(labs))
  cols[labs=="S"] = "red"
  if(plot){
    plot(x=pcx[,1],y=pcx[,2],xlab="PC1",ylab="PC2",pch=20,col=cols,lwd=2,...)
  }
  return(list(pcx=pcx,pcs=pcs,eigvals=eigvals,pcs2_var=pcs2_var))
}

par(mfrow=c(2,2))
run_pca_and_plot(pilot_datasets$`rnaseq;muscle;get`,transpose = T,
                 main="PCA, RNAseq, GET")
legend("bottomright",legend = c("R","S"),fill=c("red","blue"),cex=1.2)
run_pca_and_plot(pilot_datasets$`untar_met;muscle;umich`,transpose = T,
                 main="PCA, met untar, Umich")
run_pca_and_plot(pilot_datasets$`proteomics_ph;muscle;broad`,transpose = T,
                 main="PCA, prot ph, Broad")
run_pca_and_plot(pilot_datasets$`tar_met;muscle;duke`,transpose = T,
                 main="PCA, targeted, Duke")
dev.off()

pca2 = run_pca_and_plot(motrpac_multiomics_data2)
pca1 = run_pca_and_plot(motrpac_multiomics_data1)

omic_types1 = get_omics_vec(motrpac_multiomics_data1)
names(omic_types1) = colnames(motrpac_multiomics_data1)

get_loadings_dist_across_omics<-function(pcs,range=1:3,num=100,omic_types){
  pcts = c()
  omic_types = as.factor(omic_types)
  for(j in 1:ncol(pcs)){
    thr = sort(abs(pcs[,j]),decreasing=T)[num]
    selected = which(abs(pcs[,j])>thr)
    curr_table = table(omic_types[selected])/length(selected)
    pcts = rbind(pcts,curr_table)
  }
  return(pcts)
}

names(pca1)

# Get the different loadings and compare the methods
setwd("multi_loadings/")
pca_l = read.table("./aer_loadings.txt",sep=",")

# between omics correlation
# a series of functions for analyzing a multiomics matrix
get_omics_vec<-function(x){
  omic_types = unname(sapply(colnames(x),function(x)strsplit(x,split="\\.")[[1]][1]))
  return(omic_types)
}
get_between_omics_corr_sample<-function(x,n=100,...){
  o_types = get_omics_vec(x)
  omic2selected_ind = list()
  for(om in unique(o_types)){
    curr_log_inds = which(o_types==om)
    omic2selected_ind[[om]] = sample(curr_log_inds)[1:min(n,length(curr_log_inds))]
  }
  corrs = c()
  for(i in 2:length(omic2selected_ind)){
    om1 = names(omic2selected_ind)[i]
    for(j in 1:(i-1)){
      om2 = names(omic2selected_ind)[j]
      curr_corrs = cor(x[,omic2selected_ind[[i]]],x[,omic2selected_ind[[j]]],...)
      corrs = rbind(corrs,c(curr_corrs))
      rownames(corrs)[nrow(corrs)] = paste(om1,om2,sep=";")
    }
  }
  return(corrs)
}

get_between_omics_corr<-function(x,...){
  o_types = get_omics_vec(x)
  omic2selected_ind = list()
  for(om in unique(o_types)){
    curr_log_inds = which(o_types==om)
    omic2selected_ind[[om]] = curr_log_inds
  }
  corrs = list()
  for(i in 2:length(omic2selected_ind)){
    om1 = names(omic2selected_ind)[i]
    for(j in 1:(i-1)){
      om2 = names(omic2selected_ind)[j]
      curr_corrs = cor(x[,omic2selected_ind[[i]]],x[,omic2selected_ind[[j]]],...)
      corrs[[paste(om1,om2,sep=";")]] = curr_corrs
    }
  }
  return(corrs)
}
get_significant_cor_pairs<-function(corrs,thr=0.95){
  m = c()
  for(n1 in rownames(corrs)){
    for(n2 in colnames(corrs)){
      if(!is.na(corrs[n1,n2]) && corrs[n1,n2]>thr){
        m = rbind(m,c(n1,n2,corrs[n1,n2]))
      }
    }
  }
  return(m)
}


corrs1 = get_between_omics_corr_sample(motrpac_multiomics_data1,200,method="spearman")
par(mfrow=c(3,1))
for(nn in rownames(corrs1)){
  hist(corrs1[nn,],main=nn,xlim=c(-0.8,1))
}
corrs2 = get_between_omics_corr_sample(motrpac_multiomics_data2,200,method="spearman")
dev.off()
par(mfrow=c(3,2))
for(nn in rownames(corrs2)){
  hist(corrs2[nn,],main=nn,xlim=c(-0.8,1))
}

all_cors = get_between_omics_corr(motrpac_multiomics_data1,method="spearman")
sig_corrs = lapply(all_cors,get_significant_cor_pairs)
sapply(all_cors,dim)
sapply(sig_corrs,nrow)
sapply(all_cors,function(x)sum(x>0.95,na.rm=T))


