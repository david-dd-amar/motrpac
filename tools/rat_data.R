library("STRINGdb")
library(org.Rn.eg.db)
rat_entrez2ensemble_prot = as.list(org.Rn.egENSEMBLPROT)
rat_ensemble_prot2entrez = as.list(org.Rn.egENSEMBLPROT2EG)
rat_ensembl2entrez = as.list(org.Rn.egENSEMBL2EG)
rat_entrez2ensembl = as.list(org.Rn.egENSEMBL)
rat_proteins = unique(unname(unlist(rat_entrez2ensemble_prot)))
rat_code = 10116
string_db <- STRINGdb$new( version="10", species=rat_code,
                           score_threshold=0, input_directory="" )
rat_proteins_string_ids =  string_db$mp(rat_proteins)
rat_string_net = string_db$get_interactions(rat_proteins_string_ids)
rat_ppi = rat_string_net[rat_string_net$combined_score>400,1:2]
rm(rat_string_net)
rat_ppi[,1] = gsub("10116.","",rat_ppi[,1])
rat_ppi[,2] = gsub("10116.","",rat_ppi[,2])

get_first<-function(x){return(sort(x,decreasing = F)[1])}
rat_ensemble_prot2entrez_simple = sapply(rat_ensemble_prot2entrez,get_first)

rat_ppi_entrez = cbind(
  rat_ensemble_prot2entrez_simple[rat_ppi[,1]],
  rat_ensemble_prot2entrez_simple[rat_ppi[,2]]
)
rat_ppi_entrez = rat_ppi_entrez[!apply(is.na(rat_ppi_entrez),1,any),]

# Test propagation/diffusion
source("http://bioconductor.org/biocLite.R")
biocLite(c("foreach","doParallel"))
biocLite("supraHex")
install.packages("dnet")
library(dnet);library(foreach);library(doParallel);library(parallel)
g = igraph::graph_from_edgelist(rat_ppi_entrez)
nodes = names(V(g))
w = rep(0,length(nodes))
w[1:100]=1
w = data.frame(w)
rownames(w) = nodes
test_rw = dRWR(g,normalise="none", setSeeds=w,
               restart=0.75, parallel=FALSE)
test_rw = as.matrix(test_rw)
plot(test_rw[,1],w[,1])
visHeatmapAdv(test_rw)
# an alternative: RandomWalkRestartMH

