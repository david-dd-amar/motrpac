#' This pipeline goes through a set of phenotypic data buckets and creates
#' a merged annotation across the buckets together with a merged dataset in
#' each bucket. 
#' 
#' NOTE: We currently ignore the site-level tables. 
#' These include additional information such as temperature and humidity.
#' These tables will be provided as additional supplementary tables.
#' Interested users can use them and match their data to animal ids
#' by looking at the date of the visit (provided in the raw data).
#'
#' The pipeline was designed with the following assumptions:
#' 1. each bucket has two subdirectories: "3-Data_Sets", and "1-Data_Dictionary"
#' 2. files in each subdirectory are csvs with the same names
#' 3. field guids are unique cross-form/file/csc static ids
#' 4. Each dictionary satisfies the following:
#'    5.1 the first column in each dictionary specifies the security level
#'    of each field with 1 being safe and 4 unsafe for external release
#'    5.2 a Field.Name column exists
#' 5. The provided bucket list is ordered such that if bucket j (j = 2,3,...)
#'    has a field x then it must appear in the same csv form name as bucket j-1.
#'    In other words, a new bucket makes sure that previous field name appear
#'    in the same dictionary as in previous buckets.
#' 6. Data tables that specify longitudinal data have a field called days_visit.
#' 
#' Additional details
#' 1. Dictionary data are taken from the last form version number by default.
#' 2. BIC creates a set of field IDs, where the ID columns (e.g., viallabel, pid,
#'    bid, labelid) are enumerated first.
#' 3. Longitudinal data typically provided to BIC in a "long" format - tables in
#'    which each row corresponds to a visit of a specific animal. When creating
#'    a merged table, these data are converted to a wide format and field ids are
#'    extended to have the visit number of the animal as a suffix. For example, the field
#'    for vo2max score can have suffix "_k" where k specifies the visit number. 
#'    By taking all vo2max scores of an animal we can draw the time trajectory.
#'    This is similar to the way the UK-Biobank stores longitudinal data.

#' QA/QC tests made:
#' 1. Make sure that the number of fields in a dictionary fit the data table
#' 2. Make sure that each data csv has a dictionary

#' Specify meta-data and constants
#' -------------------------------
#' Order is important here, it should go from animal/subject id
#' into lower level ids. In addition this order specifies the
#' order of field enumeration. Not all form/file names below must
#' be available in the data transfer bucket.
#' 
#' This variable sets the id name in each csv that will be used
#' as the identifier field when the table is merged with other 
#' tables. In addition, this will also be used to check whether the
#' rows in each table are unique. As explained above, if the rows are 
#' not unique then this will trigger a procedure to transform the csv to a "wider"
#' version by looking at the days_visit. 
csv_name_main_field = c(
  "Key" = "pid",
  "Registration" = "pid",
  "Familiarization" = "pid",
  "Training" = "pid",
  "VO2.Max.Test" = "pid",
  "NMR.Testing" = "pid",
  "Terminal.Weight" = "pid",
  "Specimen.Collection" = "pid",
  "Specimen.Processing" = "labelid",
  "Calculated.Variables" = "labelid"
)


# Specify input buckets and field input
BUCKETS = c(
  "pass1b_6m" = "gs://motrpac-portal-transfer-dmaqc/DMAQC_TRANSFER_PASS_1B.6M_1.00/",
  "pass1b_18m" = "gs://motrpac-portal-transfer-dmaqc/DMAQC_TRANSFER_PASS_1B.18M_1.00/"
)
# these two variables determine how the field ids will look like
# e.g., for a prefix a and nchar 3 we will have ids "a001","a002",...
FIELD_ID_PREFIX = "pass1bf"
FIELD_ID_NCHAR = 4
# global counter for the field ids
FIELD_ID_COUNTER = 1
# specify threshold for field security for external sharing
# fields with value >= this threshold are excluded
FIELD_SECURITY_THRESHOLD = 4
# The character to use for adding a version e.g.,
# when specifying "_" we will see ids like "pass1bf001_1",
# "pass1bf001_2",... where the numeric suffix specifies the visit
# number.
FIELD_VISIT_SEP = "_"
# The id-columns in the datasets - these
# should get the first ids
ID_COLS = c(
  "pid","bid","labelid","viallabel"
)

# Load libraries and auxiliary functions
source("~/Desktop/repos/motrpac-bic-norm-qc/tools/gcp_functions.R")
library(reshape2)
# install.packages("tidyverse")
library(tidyverse)
library(stringr)

# Define helper functions to be used below
#' A way to enumerate the visits per id.
#' 
#' @details This function should be used when "spreading" out a table. 
#' For example: if Vo2 is measured twice for pid "1" then the output vector 
#' should contain 1 and 2 for the pid, where 1 marks the first visit in d_visit
#' and 2 marks the later visit. 
extract_visit_vector<-function(ids,d_visit){
  v = rep(-1,length(ids))
  names(v) = ids
  for(id in unique(ids)){
    id_inds = which(ids==id)
    id_inds_ranks = rank(d_visit[id_inds])
    v[id_inds] = id_inds_ranks
  }
  return(v)
}
# Run example:
# > extract_visit_vector(ids,d_visit)[1:4]
# 10023259 10023259 10024735 10024735 
# 2        1        1        2 
# > ids[1:4]
# [1] 10023259 10023259 10024735 10024735
# > d_visit[1:4]
# [1] 84 27 27 84
# > x$d_visit[1:4]
# [1] "28AUG2018" "02JUL2018" "02JUL2018" "28AUG2018"

#' ADD DOC here
spread_data_by_visits<-function(x,id_col = "pid",date_col = "days_visit"){
  visit_number_BIC = extract_visit_vector(x[[id_col]],x[[date_col]])
  xx = x[visit_number_BIC==1,]
  rownames(xx) = xx[[id_col]]
  non_id_inds = colnames(xx) != id_col
  colnames(xx)[non_id_inds] = paste0(colnames(xx)[non_id_inds],FIELD_VISIT_SEP,"1")
  for(i in 2:max(visit_number_BIC)){
    curr_x = x[visit_number_BIC==i,]
    rownames(curr_x) = curr_x[[id_col]]
    colnames(curr_x)[non_id_inds] = paste0(colnames(curr_x)[non_id_inds],FIELD_VISIT_SEP,i)
    curr_x = curr_x[rownames(xx),]
    rownames(curr_x) = rownames(xx)
    curr_x = curr_x[,non_id_inds]
    xx = cbind(xx,curr_x)
  }
  return(xx)
}

#' ADD doc
generate_new_field_id = function(){
  newf = paste0(FIELD_ID_PREFIX,str_pad(FIELD_ID_COUNTER,FIELD_ID_NCHAR,"left","0"))
  FIELD_ID_COUNTER = FIELD_ID_COUNTER + 1
  assign("FIELD_ID_COUNTER", FIELD_ID_COUNTER, envir = .GlobalEnv)
  return(newf)
}
names(ID_COLS) = sapply(ID_COLS,function(x)generate_new_field_id())

# v - vector
# e - entry
get_latest_ind_in_vector<-function(e,v){
  arr = which(v==e)
  return(arr[length(arr)])
}

#########################################################################################
#########################################################################################
# Load data from the from the buckets
print("Loading data from buckets")
local_path = "~/Desktop/MoTrPAC/data/pass_1b/dmaqc_pheno/"
system(paste("mkdir",local_path))
bucket2data = list()
for(bucket_name in names(BUCKETS)){
  bucket = BUCKETS[bucket_name]
  bucket2data[[bucket_name]] = list()
  print(paste("Current bucket:",bucket,", bucket name is:",bucket_name))
  dir.create(paste(local_path,bucket_name,sep="/"))
  dmaqc_data_dir = paste(local_path,bucket_name,"/3-Data_Sets/",sep="")
  dmaqc_dict_dir = paste(local_path,bucket_name,"/1-Data_Dictionary/",sep="")
  download_bucket_files_to_local_dir(
    bucket = paste(bucket,"3-Data_Sets/",sep=""),
    local_path = dmaqc_data_dir)
  download_bucket_files_to_local_dir(
    bucket = paste(bucket,"1-Data_Dictionary/",sep=""),
    local_path = dmaqc_dict_dir)
  
  all_csvs = list.files(dmaqc_data_dir,full.names = T) # get all files in dir
  all_csvs = all_csvs[grepl(".csv$",all_csvs)] # make sure we take csv only
  all_dict_csvs = list.files(dmaqc_dict_dir,full.names = T) # get all files in dir
  all_dict_csvs = all_dict_csvs[grepl(".csv$",all_dict_csvs)] # make sure we take csv only
  # read all files
  csv_data = list()
  for(fname in names(csv_name_main_field)){
    curr_data_csv = all_csvs[grepl(fname,all_csvs)]
    if(length(curr_data_csv)==0){
      print(paste("Warning: no",fname,"csv in data from bucket:",bucket))
      print("skipping")
      next
    }
    if(length(curr_data_csv)>1){
      print(paste("Error:",fname,"is not unique in downloaded data from bucket",bucket))
      print(paste("taking the first csv:",curr_data_csv[1]))
    }
    curr_data_csv = curr_data_csv[1]

    curr_dict_csv = all_dict_csvs[grepl(fname,all_dict_csvs)]
    if(length(curr_dict_csv)==0){
      print(paste("Error: no",fname,"dictionary csv in data from bucket:",bucket))
      print("skipping")
      next
    }
    if(length(curr_dict_csv)>1){
      print(paste("Error:",fname,"is not unique in downloaded data from bucket",bucket))
      print(paste("taking the first dictionary csv:",curr_dict_csv[1]))
    }
    curr_dict_csv = curr_dict_csv[1]
    
    dict = read.csv(curr_dict_csv,stringsAsFactors = F)
    dataset = read.csv(curr_data_csv,stringsAsFactors = F)
    bucket2data[[bucket_name]][[fname]] = list(dataset = dataset, dict = dict)
  }
}
print ("################  Done ##################")

print("Going over dictionaries")
for(bucket_name in names(bucket2data)){
  for(fname in names(bucket2data[[bucket_name]])){
    dataset = bucket2data[[bucket_name]][[fname]][["dataset"]]
    dict = bucket2data[[bucket_name]][[fname]][["dict"]]
    if(length(unique(dict[["Field.Name"]]))!=ncol(dataset)){
      print(paste("Error: dictionary and data do not have the same number of fields -",fname))
    }
  }
}

# Go over dictionaries to create the filtered datasets
#   Take the latest version number record
#   Remove field names using the security field
filtered_dictionary_fields = c()
for(bucket_name in names(bucket2data)){
  for(fname in names(bucket2data[[bucket_name]])){
    dataset = bucket2data[[bucket_name]][[fname]][["dataset"]]
    dict = bucket2data[[bucket_name]][[fname]][["dict"]]
    
    versionNumber_col = which(grepl("^versionNbr",names(dict),perl=T))
    security_col = which(grepl("^External.Data.Release",names(dict),perl=T))
    dict_rows_to_keep = rep(T,nrow(dict))
    if(length(security_col)>0){
      dict_rows_to_keep = dict_rows_to_keep &
        dict[[security_col]] < FIELD_SECURITY_THRESHOLD
    }
    if(length(versionNumber_col)>0){
      dict_rows_to_keep = dict_rows_to_keep &
        dict[[versionNumber_col]] == max(dict[[versionNumber_col]],na.rm=T)
    }
    dict_columns_to_keep = 
      !grepl("^External.Data.Release",names(dict),perl=T) &
      !grepl("^versionNbr",names(dict),perl=T) &
      !grepl("^Data.Set.Variable.Sequence",names(dict),perl=T)
    filtered_dict = dict[dict_rows_to_keep,dict_columns_to_keep]
    filtered_dataset = dataset[,filtered_dict[["Field.Name"]]]
    bucket2data[[bucket_name]][[fname]]$filtered_dataset = filtered_dataset
    bucket2data[[bucket_name]][[fname]]$filtered_dict = filtered_dict
    
    filtered_dictionary_fields = union(filtered_dictionary_fields,colnames(filtered_dict))
  }
}
filtered_dictionary_fields = setdiff(filtered_dictionary_fields,"Data.Set.Name")
filtered_dictionary_fields = setdiff(filtered_dictionary_fields,"fieldGUID..Field.unique.identifier.")

# QA
for(bucket_name in names(bucket2data)){
  for(fname in names(bucket2data[[bucket_name]])){
    dataset = bucket2data[[bucket_name]][[fname]][["filtered_dataset"]]
    dict = bucket2data[[bucket_name]][[fname]][["filtered_dict"]]
    if(nrow(dict)!=ncol(dataset)){
      print(paste("Error: filtered dictionary and data do not have the same number of fields -",fname))
    }
  }
}

#########################################################################################
#########################################################################################

print("Extracting unique field ids")
# Go over dictionaries to create the unique field codes
# this table keeps: the long field name (fname "_" field name),
# the bucket it came from and the latest version
field_unique_ids = cbind(ID_COLS,"","","","")
rownames(field_unique_ids) = ID_COLS
colnames(field_unique_ids) = c("FullName","Field.Name","CsvName","BucketSource","VersionNumber")
for(bucket_name in names(bucket2data)){
  for(fname in names(bucket2data[[bucket_name]])){
    dict = bucket2data[[bucket_name]][[fname]][["filtered_dict"]]
    # check id columns
    curr_id_col_inds = dict[["Field.Name"]] %in% ID_COLS
    if(any(curr_id_col_inds)){
      curr_id_cols = dict[["Field.Name"]][curr_id_col_inds]
      # if we do not have form/csv names use the current fname
      # and move on, this is important as below we will need the 
      # csv name for creating the merged dictionary
      curr_id_cols = curr_id_cols[field_unique_ids[curr_id_cols,"CsvName"]==""]
      if(length(curr_id_cols)>0){
        field_unique_ids[curr_id_cols,"CsvName"] = fname
        field_unique_ids[curr_id_cols,"BucketSource"] = bucket_name
        field_unique_ids[curr_id_cols,"Field.Name"] = curr_id_cols
      }
    }
    
    # we are now analyzing the non-id columns
    # but also take the latest version of a field
    dict = dict[!curr_id_col_inds,]
    shortfnames = dict[["Field.Name"]]
    # get the latest appearance of each field name in the dictionary
    latest_version_inds = sapply(shortfnames,get_latest_ind_in_vector,v=shortfnames)
    latest_version_inds = unique(latest_version_inds)
    # remove previous versions
    dict = dict[latest_version_inds,]
    # obtain the full field names and keep only those not
    # covered already in our field metadata table
    fullnames = paste(fname,dict[["Field.Name"]],sep=".")
    rownames(dict) = fullnames
    to_keep = ! (fullnames %in% field_unique_ids[,1])
    if(sum(to_keep)==0){next}
    dict = dict[to_keep,]
    versionNumber_col = which(grepl("^versionNbr",names(dict),perl=T))
    if(length(versionNumber_col)==0){
      curr_info = cbind(rownames(dict),dict[["Field.Name"]],
                        fname,bucket_name,"")
    }
    else{
      curr_info = cbind(rownames(dict),dict[["Field.Name"]],
                        fname,bucket_name,dict[,versionNumber_col])
    }
    field_unique_ids = rbind(field_unique_ids,curr_info)
  }
}

# add unique identifiers
field_ids = names(ID_COLS)
for(j in (1+length(ID_COLS)):nrow(field_unique_ids)){
  field_ids = c(field_ids,generate_new_field_id())
}
field_unique_ids = cbind(field_ids,field_unique_ids)
colnames(field_unique_ids)[1] = "BICUniqueID"

print("created the metadata table of the fields")
print("Here are the first 10 rows:")
write.table(field_unique_ids[1:10,],sep="\t",row.names = F)

# Go over datasets and merge - once the metadata table is created
# we can now use it to go over forms and fetch the information
# of each field.
merged_dictionary = field_unique_ids[,1:3]
for(dict_f_name in filtered_dictionary_fields){
  merged_dictionary = cbind(merged_dictionary,"")
}
colnames(merged_dictionary)[4:ncol(merged_dictionary)] = filtered_dictionary_fields
for(i in 1:nrow(merged_dictionary)){
  curr_csv = field_unique_ids[i,"CsvName"]
  curr_field = field_unique_ids[i,"Field.Name"]
  curr_bucket = field_unique_ids[i,"BucketSource"]
  curr_version = field_unique_ids[i,"VersionNumber"]
  curr_dict = bucket2data[[curr_bucket]][[curr_csv]]$filtered_dict
  curr_rows = curr_dict[
    curr_dict[["Field.Name"]] == curr_field,
  ]
  if(curr_version == "" &&
     !is.null(dim(curr_rows)) && nrow(curr_rows)>1){
    print(paste("Error: field",curr_field, "in form",curr_csv,
                "has multiple entries in dictionary but no version was parsed"))
  }
  if(!is.null(dim(curr_rows)) && nrow(curr_rows)>1){
    versionNumber_col = which(grepl("^versionNbr",names(dict),perl=T))
    curr_rows = curr_dict[
      curr_dict[["Field.Name"]] == curr_field && curr_dict[[versionNumber_col]] == curr_version,
    ]
  }
  if(length(curr_rows)==0 && curr_field!="viallabel"){
    print(paste("Error: field",curr_field,"has no dictionary entry in csv",curr_csv))
  } 
  if(length(curr_rows)==0){next}
  # fill in by taking the names of the current new entry
  curr_rows = as.matrix(curr_rows)[1,]
  curr_rows = curr_rows[intersect(names(curr_rows),colnames(merged_dictionary))]
  merged_dictionary[i,names(curr_rows)] = curr_rows
}

print("Creating the merged dictionary...")
print("done, first 10 rows:")
write.table(merged_dictionary[1:10,],sep="\t",row.names = F)

# Finally, we go over each bucket and create its merged dataset


