path1 = "/Users/David/Desktop/MoTrPAC/march_2019/pheno_data/PASS1A.6M-RLS 0,01/3-Data Sets/"
path2 = "/Users/David/Desktop/MoTrPAC/april_2019/DMAQC_Transfer_Pass_1A.6M_1/3-Data_Sets/"
files1 = list.files(path1,full.names = T)
files2 = list.files(path2,full.names = T)
for(i in 1:length(files1)){
  d1 = read.csv(files1[i],stringsAsFactors = F)
  d2 = read.csv(files2[i],stringsAsFactors = F)
  if(sum(d1!=d2,na.rm = T)>0){
    rownames(d1) = d1$labelid
    rownames(d2) = d2$labelid
    d1 = d1[rownames(d2),]
    print(sort(colSums(d1!=d2)))
  }
  
}

