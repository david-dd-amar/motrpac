################ GO enrichment analysis ################

# GO enrichment
library(topGO)
run_topgo_enrichment_fisher<-function(genesOfInterest,geneUniverse,
         gene2name = NULL,
         go_term_size=10,go_dags=c("BP","MF"),pthr=0.05,
         mapping="org.Hs.eg.db",ID = "entrez", ...){
  l = list()
  if(class(genesOfInterest)=="character"){
    l[["set1"]] = genesOfInterest
  }
  else{
    l = genesOfInterest
  }
  res = c()
  for(nn in names(l)){
    genesOfInterest = l[[nn]]
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    for(type in go_dags){
      myGOdata <- new("topGOdata", description="My project", ontology=type, 
                      nodeSize=go_term_size,
                      allGenes=geneList, annot = annFUN.org, mapping=mapping, ID = ID)
      allGOs = usedGO(myGOdata)
      resultFisher <- runTest(myGOdata, algorithm="classic", statistic="fisher")
      allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = 
                           "resultFisher", ranksOf = "classicFisher",
                         topNodes = length(score(resultFisher)))
      
      # add selected genes, only for enrichments with p < 0.05
      gene_rep = c()
      for(i in 1:nrow(allRes)){
        currp = as.numeric(allRes[i,"classicFisher"])
        if (is.na(currp)){
          allRes[i,"classicFisher"] = 1e-30
          currp = allRes[i,"classicFisher"] 
        }
        if ( is.na(currp) || currp > pthr){
          gene_rep[i]=""
        }
        else{
          curr_go = allRes[i,"GO.ID"]
          curr_genes = genesInTerm(myGOdata,curr_go)[[1]]
          curr_genes = intersect(curr_genes,l[[nn]])
          if(length(curr_genes)!=as.numeric(allRes[i,"Significant"])){
            print("Error in intersect, debug the topGO wrapper!")
            return(NULL)
          }
          if(!is.null(gene2name)){
            curr_genes = unlist(gene2name[curr_genes])
          }
          gene_rep[i] = paste(curr_genes,collapse=",")
        }
      }
      
      setname = rep(nn,nrow(allRes))
      typev = rep(type,nrow(allRes))
      allRes = cbind(setname,typev,allRes,gene_rep)
      res = rbind(res,allRes)
    }
  }
  go_pvals = as.numeric(res[,ncol(allRes)])
  go_pvals[is.na(go_pvals)] = 1
  go_qvals = p.adjust(go_pvals,method='fdr')
  res = cbind(res,go_qvals)
  return(res)
}
extract_top_go_results<-function(res,qval=0.1,maxsize=2000){
  res = res[,-which(colnames(res)=="typev")]
  res = res[res$Annotated<maxsize,]
  res = unique(res)
  res[res$GO.ID=="GO:0002532",]
  new_qvals = p.adjust(res$classicFisher,method='fdr')
  #plot(new_qvals,res$go_qvals);abline(0,1)
  res$go_qvals = new_qvals
  return(res[res$go_qvals<=qval,])
}
get_most_sig_enrichments_by_groups <- function(res,num=1,gcol=1,pcol="classicFisher"){
  gs = unique(as.character(res[,1]))
  m = c()
  for(g in gs){
    res0 = res[res[,gcol]==g,]
    ps = as.numeric(res0[,pcol])
    thr = sort(ps,decreasing = F)[min(num,length(ps))]
    m = rbind(m,res0[ps<=thr,])
  }
  return(m)
}

get_most_sig_enrichments_for_cluster <- function(res,cluster,...){
  res = res[res[,1]==cluster,]
  return(get_most_sig_enrichments_by_groups(res,...))
}

# GSEA
library(fgsea)
fgsea_wrapper <- function(pathways,scores,nperm=2000,run_nperm=1000,...){
  num_runs = nperm/run_nperm
  l = list()
  for(j in 1:num_runs){
    l[[j]] = fgsea(pathways,scores,nperm = run_nperm,...)
  }
  emp_pvals = sapply(l,function(x)x$pval)
  emp_pvals = emp_pvals*run_nperm
  min_to_add = min(emp_pvals)
  emp_pvals = emp_pvals-min_to_add
  new_pvals = rowSums(emp_pvals)+min_to_add
  new_pvals = new_pvals/nperm
  new_qvals = p.adjust(new_pvals,method='fdr')
  res = l[[1]]
  res[,"pval"] = new_pvals
  res[,"padj"] = new_qvals
  return(res)
}

get_matrix_p_adjust<-function(x,q=0.1,...){
  v = c(x)
  v = v[!is.na(x)]
  vq = p.adjust(v,...)
  thr = max(v[vq<=q])
  return(thr)
}

