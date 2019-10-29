
tryCatch({library(cloudml)}, error = function(e) {
  print("Cannot load cloudml, please install")
})
tryCatch({library(data.table)}, error = function(e) {
  print("Cannot load data.table, please install")
})

tissue_code2name = c(
  "T30"="PaxGene_RNA",
  "T31"="Plasma",
  "T32"="Packed_Cells",
  "T33"="Hippocampus",
  "T34"="Cortex",
  "T35"="Hypothalamus",
  "T36"="Gastrocnemius",
  "T37"="Vastus_Lateralis",
  "T38"="Tibia",
  "T39"="Heart",
  "T40"="Kidney",
  "T41"="Adrenals",
  "T42"="Colon",
  "T43"="Spleen",
  "T44"="Testes",
  "T45"="Ovaries",
  "T46"="Brown_Adipose",
  "T47"="White_Adipose",
  "T48"="Aorta",
  "T49"="Lung",
  "T50"="Small_Intestine",
  "T51"="Liver",
  "T52"="Hippocampus_Powder",
  "T53"="Cortex_Powder",
  "T54"="Hypothalamus_Powder",
  "T55"="Gastrocnemius_Powder",
  "T56"="Vastus_Lateralis_Powder",
  "T57"="Tibia_Powder",
  "T58"="Heart_Powder",
  "T59"="Kidney_Powder",
  "T60"="Adrenal_Powder",
  "T61"="Colon_Powder",
  "T62"="Spleen_Powder",
  "T63"="Testes_Powder",
  "T64"="Ovaries_Powder",
  "T65"="Aorta_Powder",
  "T66"="Lung_Powder",
  "T67"="Small_Intestine_Powder",
  "T68"="Liver_Powder",
  "T69"="Brown_Adipose_Powder",
  "T70"="White_Adipose_Powder"
)

load_from_bucket<-function(file,bucket,delete=T){
  system(paste("~/google-cloud-sdk/bin/gsutil cp",
               paste(bucket,file,sep=""),
               getwd()))
  objects_before = ls()
  load(file)
  objects_after = ls()
  if(delete){
    system(paste("rm",file))
  }
  added_objects = setdiff(objects_after,objects_before)
  added_objects = setdiff(added_objects,"objects_before")
  res = list()
  for(obj_name in added_objects){
    res[[obj_name]] = get(obj_name)
  }
  return(res)
}

download_bucket_files_to_local_dir<-function(bucket,local_path=NULL){
  if(is.null(local_path)){local_path=paste(getwd(),"gs_download",sep="/")}
  # download files
  # remove old data if exists
  if(file.exists(local_path)){system(paste("rm -r",local_path))}
  # download
  system(paste("mkdir",local_path))
  cmd = paste("~/google-cloud-sdk/bin/gsutil -m cp -r",paste(bucket,"*",sep=""),local_path)
  system(cmd)
  # get the data
  downloaded_files = list.files(local_path,full.names = T)
  return(list(downloaded_files=downloaded_files,local_path=local_path))
}

save_to_bucket<-function(...,file,bucket){
  save(...,file = file)
  system(paste("~/google-cloud-sdk/bin/gsutil cp",
               file, bucket))
  system(paste("rm",file))
}

get_files_in_bucket<-function(bucket){
  tmp_name = paste("tmp",runif(1,1,1000),".txt",sep="")
  system(paste("~/google-cloud-sdk/bin/gsutil ls",
               bucket,">",tmp_name))
  files = readLines(tmp_name)
  system(paste("rm",tmp_name))
  return(files)
}

read_metabolomics_datasets_from_release_bucket<-function(bucket,local_path=NULL,delete=T){
  metabolomics_parsed_datasets = list()
  site_buckets = get_files_in_bucket(bucket)
  names(site_buckets) = gsub(bucket,replacement = "",site_buckets)
  names(site_buckets) = gsub("/",replacement = "",names(site_buckets))
  for(site in names(site_buckets)){
    site_bucket = site_buckets[site]
    tissues_in_bucket = get_files_in_bucket(site_bucket)
    tissues_in_bucket = tissues_in_bucket[!grepl("qc_report",tissues_in_bucket)]
    names(tissues_in_bucket) = gsub(site_bucket,replacement = "",tissues_in_bucket)
    names(tissues_in_bucket) = gsub("/",replacement = "",names(tissues_in_bucket))
    for(tissue in names(tissues_in_bucket)){
      tissue_bucket = tissues_in_bucket[tissue]
      platforms_in_bucket = get_files_in_bucket(tissue_bucket)
      names(platforms_in_bucket) = gsub(tissue_bucket,replacement = "",platforms_in_bucket)
      names(platforms_in_bucket) = gsub("/",replacement = "",names(platforms_in_bucket))
      for(platform in names(platforms_in_bucket)){
        platform_bucket = platforms_in_bucket[platform]
        platform_bucket = get_files_in_bucket(platform_bucket)[1]
        platform_datasets = get_files_in_bucket(platform_bucket)
        named_bucket = paste(platform_bucket,"NAMED/",sep="")
        named_data = read_single_metabolomics_dataset(named_bucket)
        unnamed_bucket = paste(platform_bucket,"UNNAMED/",sep="")
        dataset_name = paste(c(site,tissue_code2name[tissue],platform),collapse=",")
        if(is.element(unnamed_bucket,set=platform_datasets)){
          unnamed_data = read_single_metabolomics_dataset(unnamed_bucket)
          merged_data = merge_metabolomics_datasets(named_data,unnamed_data)
          metabolomics_parsed_datasets[[dataset_name]] = merged_data
        }
        else{
          metabolomics_parsed_datasets[[dataset_name]] = named_data
        }
        print(dataset_name)
        print(dim(metabolomics_parsed_datasets[[dataset_name]]$sample_data))
      }
    }
  }
  return(metabolomics_parsed_datasets)
}
# # test
# bucket = "gs://motrpac-internal-release1-1-results/metabolomics/"

#' Read the data from the release
#' @description read the data of a single folder in the data release; remove data rows that are all NAs or missing
#' @description assumes that bucket has exactly four files: experimental details, sample metadata, metabolite metadata, and the dataset called the results file
read_single_metabolomics_dataset<-function(bucket,local_path=NULL,delete=T){
  download_bucket = download_bucket_files_to_local_dir(bucket,local_path)
  local_path = download_bucket[[2]]
  downloaded_files = download_bucket[[1]]
  downloaded_files = downloaded_files[!(
    grepl("experimental",downloaded_files,ignore.case = T) &
      grepl("details",downloaded_files,ignore.case = T))]
  sample_meta_file = downloaded_files[
    grepl("metadata",downloaded_files,ignore.case = T) & 
      grepl("sample",downloaded_files)]
  downloaded_files = setdiff(downloaded_files,sample_meta_file)
  results_file = downloaded_files[grepl("results",downloaded_files,ignore.case = T)]
  downloaded_files = setdiff(downloaded_files,results_file)
  metabolite_meta_file = downloaded_files[1]
  
  sample_meta = fread(sample_meta_file,data.table = F,stringsAsFactors = F)
  row_annot = fread(metabolite_meta_file,data.table = F,stringsAsFactors = F)
  to_rem = apply(row_annot=="" | is.na(row_annot),1,all)
  row_annot = row_annot[!to_rem,]
  to_rem = apply(row_annot=="" | is.na(row_annot),2,all)
  row_annot = row_annot[,!to_rem]
  possible_row_names = paste(row_annot[,"mz"],row_annot[,"rt"],sep="_")
  if(all(table(possible_row_names)==1)){
    rownames(row_annot) = possible_row_names
  }
  sample_data = fread(results_file,data.table = F,stringsAsFactors = F)
  to_rem = apply(sample_data=="" | is.na(sample_data),1,all)
  sample_data = sample_data[!to_rem,]
  to_rem = apply(sample_data=="" | is.na(sample_data),2,all)
  sample_data = sample_data[,!to_rem]
  
  if(nrow(sample_data)!=nrow(row_annot)){
    print("Error: row annotation does not have the same number of rows as the data table!")
    return(NULL)
  }
  if(!all(sample_data[,1]==row_annot[,1])){
    print("Error: row annotation data order is not the same as in the dataset")
    return(NULL)
  }
  
  rownames(sample_data) = rownames(row_annot)
  sample_data = sample_data[,-1]
  
  # reorganaize the data
  sample_ids = sample_meta[sample_meta[,2]=="Sample",1]
  # QA: all samples in meta should be in the data
  if(!all(is.element(sample_ids,set=colnames(sample_data)))){
    print("Error: not all samples in metadata appear in the dataset")
    return(NULL)
  }
  control_data = sample_data[,!is.element(colnames(sample_data),set=sample_ids)]
  sample_data = sample_data[,sample_ids]
  is_untargeted = grepl("unnamed",bucket,ignore.case = T)
  # delete
  if(delete){system(paste("rm -r",local_path))}
  
  res = list(sample_data=sample_data,control_data=control_data,row_annot=row_annot,
             sample_meta = sample_meta,is_untargeted = is_untargeted)
  return(res)
}
# # test
# bucket = "gs://motrpac-internal-release1-1-results/metabolomics/gtech/T68/RPNEG/PROCESSED/UNNAMED/"
# curr_dataset = read_single_metabolomics_dataset(bucket)
# boxplot(log(1+curr_dataset$control_data))

#' Merge two datasets over the same sample set
merge_metabolomics_datasets<-function(named_data,unnamed_data){
  if(!all(is.element(colnames(named_data$sample_data),
                     set = colnames(unnamed_data$sample_data)))){
    print("Error: not all samples in named data appear in unnamed data")
    return(NULL)
    
  }
  if(!all(is.element(colnames(unnamed_data$sample_data),
                     set = colnames(named_data$sample_data)))){
    print("Error: not all samples in unnamed data appear in named data")
    return(NULL)
  }
  if(!all(is.element(colnames(named_data$control_data),
                     set = colnames(unnamed_data$control_data)))){
    print("Error: not all control samples in named data appear in unnamed data")
    return(NULL)
  }
  if(!all(is.element(colnames(unnamed_data$control_data),
                     set = colnames(named_data$control_data)))){
    print("Error: not all control samples in unnamed data appear in named data")
    return(NULL)
  }
  
  s1 = named_data$sample_data
  s2 = unnamed_data$sample_data
  sample_data = rbind(s1[,colnames(s2)],s2)
  
  r1 = named_data$row_annot
  r2 = unnamed_data$row_annot
  for(colname in setdiff(colnames(r1),colnames(r2))){
    r2[colname] = NA
  }
  row_annot = rbind(r1,r2[,colnames(r1)])
  
  c1 = named_data$control_data
  c2 = unnamed_data$control_data
  control_data = rbind(c1[,colnames(c2)],c2)
  
  sample_meta1 = named_data$sample_meta
  sample_meta2 = unnamed_data$sample_meta
  if(!all(sample_meta1==sample_meta2)){
    print("Error: sample metadata is not the same in both datasets, cannot merge")
    return(NULL)
  }
  
  res = list(sample_data=sample_data,control_data=control_data,row_annot=row_annot,
             sample_meta = sample_meta1,is_untargeted = T)
  return(res)
}

# # test
# bucket1 = "gs://motrpac-internal-release1-1-results/metabolomics/gtech/T68/RPNEG/PROCESSED/UNNAMED/"
# unnamed_data = read_single_metabolomics_dataset(bucket1)
# bucket2 = "gs://motrpac-internal-release1-1-results/metabolomics/gtech/T68/RPNEG/PROCESSED/NAMED/"
# named_data = read_single_metabolomics_dataset(bucket2)
# merged_data = merge_metabolomics_datasets(named_data,unnamed_data)
# boxplot(log(1+merged_data$control_data))

read_proteomics_datasets_from_release_bucket<-function(bucket,local_path=NULL,delete=T){
  proteomics_parsed_datasets = list()
  site_buckets = get_files_in_bucket(bucket)
  names(site_buckets) = gsub(bucket,replacement = "",site_buckets)
  names(site_buckets) = gsub("/",replacement = "",names(site_buckets))
  for(site in names(site_buckets)){
    site_bucket = site_buckets[site]
    tissues_in_bucket = get_files_in_bucket(site_bucket)
    tissues_in_bucket = tissues_in_bucket[!grepl("qc_report",tissues_in_bucket)]
    names(tissues_in_bucket) = gsub(site_bucket,replacement = "",tissues_in_bucket)
    names(tissues_in_bucket) = gsub("/",replacement = "",names(tissues_in_bucket))
    for(tissue in names(tissues_in_bucket)){
      tissue_bucket = tissues_in_bucket[tissue]
      platforms_in_bucket = get_files_in_bucket(tissue_bucket)
      names(platforms_in_bucket) = gsub(tissue_bucket,replacement = "",platforms_in_bucket)
      names(platforms_in_bucket) = gsub("/",replacement = "",names(platforms_in_bucket))
      for(platform in names(platforms_in_bucket)){
        platform_bucket = platforms_in_bucket[platform]
        platform_bucket = get_files_in_bucket(platform_bucket)[1]
        print(length(get_files_in_bucket(platform_bucket)))
        dataset_name = paste(c(site,tissue_code2name[tissue],platform),collapse=",")
        dataset = read_single_proteomics_dataset(platform_bucket)
        proteomics_parsed_datasets[[dataset_name]] = dataset
        print(dataset_name)
        print(dim(dataset$ratio_data))
      }
    }
  }
  return(proteomics_parsed_datasets)
}
# # test
# bucket = "gs://motrpac-internal-release1-1-results/proteomics"
# save_to_bucket(proteomics_parsed_datasets,file="proteomics_parsed_datasets.RData",
#                 bucket = "gs://bic_data_analysis/pass1a/")

read_single_proteomics_dataset<-function(bucket,local_path=NULL,delete=T){
  download_bucket = download_bucket_files_to_local_dir(bucket,local_path)
  local_path = download_bucket[[2]]
  downloaded_files = download_bucket[[1]]
  
  # ratio data
  l1 = read_proteomics_subdataset_by_regex(downloaded_files,"ratio")
  if(!is.null(l1)){names(l1) = paste("ratio",names(l1),sep="_")}
  # intensity data
  l2 = read_proteomics_subdataset_by_regex(downloaded_files,"rii")
  if(!is.null(l2)){names(l2) = paste("rii",names(l2),sep="_")}
  
  if(delete){system(paste("rm -r",local_path))}
  
  merged = c(l1,l2)
  return(merged)
}

read_proteomics_subdataset_by_regex<-function(downloaded_files,regex="ratio"){
  data_file = downloaded_files[grepl(regex,downloaded_files) &
                                 grepl("results",downloaded_files)]
  metadata_file = downloaded_files[grepl("vial",downloaded_files) &
                                     grepl("metadata",downloaded_files)]
  if(length(metadata_file)==0){
    metadata_file = downloaded_files[grepl("metadata",downloaded_files)]
  }
  
  data = fread(data_file,data.table = F,stringsAsFactors = F,header = T)
  row_annot_cols = which(!sapply(data,is.numeric))
  row_annot = as.matrix(data[,row_annot_cols])
  data = data[,-row_annot_cols]
  
  sample_meta = fread(metadata_file,data.table = F,stringsAsFactors = F,header = T)
  rownames(sample_meta) = sample_meta[,1]
  if(!all(is.element(colnames(data)[-c(1:2)],set=rownames(sample_meta)))){
    print("Error in ratio data, not all samples are in the metadata file")
    return(NULL)
  }
  return(list(data=data,sample_meta=sample_meta,row_annot=row_annot))
}
# #test
# bucket =  "gs://motrpac-internal-release1-1-results/proteomics/broad/T58/pr/results/"
# res = read_single_proteomics_dataset(bucket)




