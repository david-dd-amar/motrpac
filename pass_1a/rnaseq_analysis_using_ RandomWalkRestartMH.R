if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RandomWalkRestartMH")
library(igraph)
library("RandomWalkRestartMH")

source("~/Desktop/repos/motrpac/tools/gcp_functions.R")
source("~/Desktop/repos/motrpac/tools/rat_data.R")

pairwise_tests_results = load_from_bucket(
  "tissue2diff_ttest_and_mediation.RData",
  "gs://bic_data_analysis/pass1a/rnaseq/"
)
t_tests_results = pairwise_tests_results[["t_tests_results"]]

all_ttest_pvals = lapply(t_tests_results,function(x)x[,3])
all_pvals = unname(unlist(all_ttest_pvals))
all_qvals = p.adjust(all_pvals,method="BY")
pval_thr = max(all_pvals[all_qvals<0.1])
all_ttest_fcs = lapply(t_tests_results,function(x)x[,4])

merge_matrix<-function(l,missing_val = 1){
  objs = unique(unname(unlist(sapply(l,names))))
  m = matrix(missing_val,nrow=length(objs),ncol=length(l),
             dimnames = list(objs,names(l)))
  for(j in 1:length(l)){
    v = l[[j]]
    m[names(v),j] = v
  }
  return(m)
}

all_ttest_pvals = merge_matrix(all_ttest_pvals)
all_ttest_fcs = merge_matrix(all_ttest_fcs,0)

g = igraph::graph_from_edgelist(rat_ppi_entrez)
nodes = names(V(g))
res_mat = all_ttest_pvals < pval_thr
rownames(res_mat) = rat_ensembl2entrez_simple[rownames(res_mat)]
to_rem = is.na(rownames(res_mat)) | !is.element(rownames(res_mat),set=nodes)
res_mat = res_mat[!to_rem,]
currgenes = rownames(res_mat)
res_mat = apply(res_mat,2,function(x,y)tapply(x,INDEX = y,FUN = sum),y=currgenes)
mode(res_mat) = "numeric"
res_mat = as.data.frame(res_mat,check.names = F)
g = igraph::subgraph(g,rownames(res_mat))

# combine the data into one network
merged_net = igraph::as_edgelist(g)
gene_interactions = c()
for(nn in colnames(res_mat)){
  curr_genes = rownames(res_mat)[res_mat[,nn]>0]
  curr_edges = cbind(curr_genes ,rep(nn,length(curr_genes)))
  print(dim(curr_edges))
  merged_net = rbind(merged_net,curr_edges)
  gene_interactions = rbind(gene_interactions, curr_edges)
}
gene_interactions = data.frame(gene_interactions)
merged_g = igraph::graph_from_edgelist(merged_net,directed = F)
nodes = names(V(merged_g))
seeds = matrix(0,nrow=length(nodes),ncol=ncol(res_mat))
rownames(seeds) = nodes
colnames(seeds) = colnames(res_mat)
for(nn in names(res_mat)){seeds[nn,nn]=1}
seeds[names(res_mat)] = 1
# set attribute for node type
nodetype = rep("gene",length(seeds))
names(nodetype) = names(seeds)
nodetype[names(res_mat)] = "experiment"
nodetype = nodetype[names(V(merged_g))]
merged_g = set_vertex_attr(merged_g,"nodetype", index = V(merged_g), nodetype)

# Run the method
PPI_MultiplexObject <- create.multiplex(merged_g,Layers_Name=c("PPI"))
AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)
prop_results = c()
for(nn in colnames(res_mat)){
	seed_res = Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI, PPI_MultiplexObject,nn)
	seed_res = seed_res[[1]]
	v = seed_res[,2]
	names(v) = seed_res[,1]
	v = v[colnames(res_mat)]
	prop_results = rbind(prop_results,v)
	rownames(prop_results)[nrow(prop_results)] = nn
}
prop_results[is.na(prop_results)]=0
for(i in 1:nrow(prop_results)){
	prop_results[i,] = prop_results[i,]/(sum(prop_results[i,]))
}

library(corrplot)
corrplot(prop_results,is.corr=F,tl.cex=0.5)

#################################
# Run the multiplex approach
#################################
PPI_MultiplexObject <- create.multiplex(g,Layers_Name=c("PPI"))
AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)
disease_g = igraph::graph_from_edgelist(cbind(colnames(res_mat),colnames(res_mat)))

PPI_Disease_Net <- create.multiplexHet(PPI_MultiplexObject, disease_g, gene_interactions, c("Disease"))
PPIHetTranMatrix <- compute.transition.matrix(PPI_Disease_Net)

prop_results = c()
for(nn in colnames(res_mat)){
	currgenes = rownames(res_mat)[res_mat[,nn]>0]
	seed_res = Random.Walk.Restart.MultiplexHet(PPIHetTranMatrix,PPI_Disease_Net, currgenes,nn)
	seed_res = seed_res[[1]]
	v = seed_res[,2]
	names(v) = seed_res[,1]
	v = v[colnames(res_mat)]
	prop_results = rbind(prop_results,v)
	rownames(prop_results)[nrow(prop_results)] = nn
}
prop_results[is.na(prop_results)]=0
for(i in 1:nrow(prop_results)){
	prop_results[i,] = prop_results[i,]/(sum(prop_results[i,]))
}


# From the package's tutorial
# data(PPI_Network)
# class(PPI_Network)
# PPI_MultiplexObject <- create.multiplex(PPI_Network,Layers_Name=c("PPI"))
# AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
# AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)
# Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI, + PPI_MultiplexObject,SeedGene)


