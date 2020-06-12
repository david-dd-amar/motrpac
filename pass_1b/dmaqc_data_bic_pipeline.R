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
#' csv_name_main_field sets the id name in each csv that will be used
#' as the identifier field when the table is merged with the previous 
#' tables in this ordered list. 
#' In addition, this will also be used to check whether the
#' rows in each table are unique. As explained above, if the rows are 
#' not unique then this will trigger a procedure to transform the csv to a "wider"
#' version by looking at the days_visit. 

#' # PASS1B parameters:
#' #-------------------------------------------------
#' csv_name_main_field = c(
#'   "Key" = "pid",
#'   "Registration" = "pid",
#'   "Familiarization" = "pid",
#'   "Training" = "pid",
#'   "VO2.Max.Test" = "pid",
#'   "NMR.Testing" = "pid",
#'   "Terminal.Weight" = "pid",
#'   "Specimen.Collection" = "pid",
#'   "Specimen.Processing" = "pid",
#'   "Calculated.Variables" = "labelid",
#'   # This is the name of the unique csv that provides mapping between
#'   # all ids, especially between label ids and viallabel.
#'   "BICLabelData" = "labelid"
#' )
#' 
#' #' Specify which tables need to be converted into a
#' #' wide format so that each row will correspond to a single
#' #' animal/biospecimen/viallabel
#' longitudinal_csv2id_col = c(
#'   "VO2.Max.Test" = "pid",
#'   "NMR.Testing" = "pid"
#' )
#' 
#' # Specify input buckets and field input
#' BUCKETS = c(
#'   "pass1b_6m" = "gs://motrpac-portal-transfer-dmaqc/DMAQC_TRANSFER_PASS_1B.6M_1.00/",
#'   "pass1b_18m" = "gs://motrpac-portal-transfer-dmaqc/DMAQC_TRANSFER_PASS_1B.18M_1.00/"
#' )
#' 
#' # All output files are written to this bucket:
#' OUTBUCKET = "gs://bic_data_analysis/pass1b/pheno_dmaqc/merged"
#' 
#' # these two variables determine how the field ids will look like
#' # e.g., for a prefix a and nchar 3 we will have ids "a001","a002",...
#' FIELD_ID_PREFIX = "pass1bf"
#' 
#' # local dir for storing the downloaded data
#' local_path = "~/Desktop/MoTrPAC/data/pass_1b/dmaqc_pheno/"

# PASS1A parameters:
#-------------------------------------------------
csv_name_main_field = c(
  "Key" = "pid",
  "Registration" = "pid",
  "Familiarization" = "pid",
  "Acute.Test" = "pid",
  "Specimen.Collection" = "pid",
  "Specimen.Processing" = "pid",
  "Calculated.Variables" = "labelid",
  # This is the name of the unique csv that provides mapping between
  # all ids, especially between label ids and viallabel.
  "BICLabelData" = "labelid"
)

#' Specify which tables need to be converted into a
#' wide format so that each row will correspond to a single
#' animal/biospecimen/viallabel
longitudinal_csv2id_col = c()

# Specify input buckets and field input
BUCKETS = c(
  "pass1a_6m" = "gs://motrpac-portal-transfer-dmaqc/DMAQC_Transfer_PASS_1A.6M_1.03/",
  "pass1a_18m" = "gs://motrpac-portal-transfer-dmaqc/DMAQC_TRANSFER_PASS_1A.18M_1.00/"
)

# All output files are written to this bucket:
OUTBUCKET = "gs://bic_data_analysis/pass1a/pheno_dmaqc/merged"

# these two variables determine how the field ids will look like
# e.g., for a prefix a and nchar 3 we will have ids "a001","a002",...
FIELD_ID_PREFIX = "pass1af"

local_path = "~/Desktop/MoTrPAC/data/pass_1a/dmaqc_pheno/"

###########################
# Other general parameters
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

# Spec for the extra information provided in the
# BICLabelData csv
# Specify the column name of the CAS site
SITE_FIELD_NAME = "shiptositeid"

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

#' Spread a longitudinal table such that each subject/specimen id will have exactly one row 
#' @details for the first visit we take the data as is (i.e., subset of the rows). For subsequent visits
#' we take the non id columns and append them to the table using cbind. 
spread_data_by_visits<-function(x,id_col = "pid",date_col = "days_visit",other_id_inds = c()){
  visit_number_BIC = extract_visit_vector(x[[id_col]],x[[date_col]])
  xx = x[visit_number_BIC==1,]
  rownames(xx) = xx[[id_col]]
  non_id_inds = ! (colnames(xx) %in% union(id_col,other_id_inds))
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

#' Use the global counter FIELD_ID_COUNTER to generate
#' a new field ID and increment it by 1
generate_new_field_id = function(){
  newf = paste0(FIELD_ID_PREFIX,str_pad(FIELD_ID_COUNTER,FIELD_ID_NCHAR,"left","0"))
  FIELD_ID_COUNTER = FIELD_ID_COUNTER + 1
  assign("FIELD_ID_COUNTER", FIELD_ID_COUNTER, envir = .GlobalEnv)
  return(newf)
}
names(ID_COLS) = sapply(ID_COLS,function(x)generate_new_field_id())

# Auxiliary function to get the index of the last appearance of e in v
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
print("--------------------------------------------------------------")
system(paste("mkdir",local_path))
bucket2data = list()
for(bucket_name in names(BUCKETS)){
  bucket = BUCKETS[bucket_name]
  bucket2data[[bucket_name]] = list()
  print(paste("Current bucket:",bucket,", bucket name is:",bucket_name))
  curr_local_dir = paste0(local_path,bucket_name)
  dir.create(curr_local_dir)
  curr_obj = download_bucket_to_local_dir(bucket,curr_local_dir)
  curr_files = curr_obj$downloaded_files
  
  all_csvs = curr_files[grepl("3-Data_Sets",curr_files)] # get all dataset files in dir
  all_csvs = all_csvs[grepl(".csv$",all_csvs)] # make sure we take csv only
  all_dict_csvs = curr_files[grepl("1-Data_Dictionary",curr_files)] # get all dict files in dir
  all_dict_csvs = all_dict_csvs[grepl(".csv$",all_dict_csvs)] # make sure we take csv only
  # read all files
  csv_data = list()
  for(fname in names(csv_name_main_field)){
    curr_data_csv = all_csvs[grepl(fname,all_csvs)]
    if(length(curr_data_csv)==0){
      print(paste("Warning: no",fname,"csv in data from bucket:",bucket))
      print("skipping this csv but keeping the others")
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
    colnames(dataset) = tolower(colnames(dataset))
    dict[["Field.Name"]] = tolower(dict[["Field.Name"]])
    bucket2data[[bucket_name]][[fname]] = list(dataset = dataset, dict = dict)
  }
}
print ("Done")

print("Going over dictionaries")
print("--------------------------------------------------------------")
print("Verifying that all field names appear in the dictionaries")
for(bucket_name in names(bucket2data)){
  for(fname in names(bucket2data[[bucket_name]])){
    dataset = bucket2data[[bucket_name]][[fname]][["dataset"]]
    dict = bucket2data[[bucket_name]][[fname]][["dict"]]
    if(length(unique(dict[["Field.Name"]]))!=ncol(dataset)){
      print(paste("Error: dictionary and data do not have the same number of fields -",fname))
    }
  }
}
print("Done")

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
print("--------------------------------------------------------------")
# Go over dictionaries to create the unique field codes
# this table keeps: the long field name (fname "_" field name),
# the bucket it came from and the latest version
field_unique_ids = cbind(ID_COLS,ID_COLS,"","","")
rownames(field_unique_ids) = ID_COLS
colnames(field_unique_ids) = c("FullName","Field.Name","CsvName","BucketSource","VersionNumber")
for(bucket_name in names(bucket2data)){
  # take all forms but remove the BIC's mapping one
  formnames = names(bucket2data[[bucket_name]])
  formnames = setdiff(formnames,"BICLabelData")
  formnames = intersect(formnames,names(bucket2data[[bucket_name]]))
  for(fname in formnames){
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
write.table(field_unique_ids[1:10,],sep="\t",row.names = F,quote=F)

# Go over datasets and merge - once the metadata table is created
# we can now use it to go over forms and fetch the information
# of each field.
merged_dictionary = field_unique_ids[,1:4]
for(dict_f_name in setdiff(filtered_dictionary_fields,"Field.Name")){
  merged_dictionary = cbind(merged_dictionary,"")
}
colnames(merged_dictionary)[5:ncol(merged_dictionary)] = 
  setdiff(filtered_dictionary_fields,"Field.Name")
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
    break
  } 
  if(length(curr_rows)==0){next}
  # fill in by taking the names of the current new entry
  curr_rows = as.matrix(curr_rows)[1,]
  curr_rows = curr_rows[intersect(names(curr_rows),colnames(merged_dictionary))]
  merged_dictionary[i,names(curr_rows)] = curr_rows
}

# Add the form names for the calculated variables
calcvar_inds = grepl("Calculated",merged_dictionary[,"FullName"])
merged_dictionary[calcvar_inds,"FormName"] = "Calculated.Variables"

print("Creating the merged dictionary...")
print("Done, first 10 rows and 10 columns:")
rownames(merged_dictionary) = merged_dictionary[,1]
write.table(merged_dictionary[1:10,1:6],sep="\t",row.names = F,quote=F)

#########################################################################################
#########################################################################################
# Finally, we go over each bucket and create its merged dataset.
# This code merges the different tables using their respective
# main id. At these point none of them are vial labels. These are
# added in the last step below as these require special handling 
# because their field information is not in a standard csv file.

# A helper object to map long field names to their code/id
longname2fid = merged_dictionary[,1]
names(longname2fid) = merged_dictionary[,2]

# This list will keep the animal level data and the viallabel 
# based data for each bucket.
# Comment: this is a naive, inefficient implementation. However,
# it is extremely simple and is less error prone, so we keep it.
print("Merging data tables")
print("--------------------------------------------------------------")
bucket2merged_data = list()
for(bucket_name in names(bucket2data)){
  print(paste("analyzing data from bucket:",bucket_name))
  m = c()
  bucket2merged_data[[bucket_name]] = list()
  # take all forms but remove the BIC's mapping one
  formnames = names(bucket2data[[bucket_name]])
  formnames = setdiff(formnames,"BICLabelData")
  formnames = intersect(formnames,names(bucket2data[[bucket_name]]))
  for(fname in formnames){
    dataset = bucket2data[[bucket_name]][[fname]][["filtered_dataset"]]
    by_col = csv_name_main_field[fname]
    if(any(is.na(dataset[,by_col]))){
      print(paste("Warning: id column",merged_dictionary[longname2fid[by_col],"Field.Name"],"has NA values in",
                  fname,"removing these rows"))
      dataset = dataset[!is.na(dataset[,by_col]),]
    }
    # convert column names to their longer names
    non_id_inds = which(! colnames(dataset) %in% ID_COLS)
    colnames(dataset)[non_id_inds] = paste0(fname,".",colnames(dataset)[non_id_inds])
    # get the index of the days visit col
    days_visit_col = which(grepl("days_visit",colnames(dataset),ignore.case = T))
    # convert to field id (only for the relevant ones)
    in_fname = colnames(dataset) %in% names(longname2fid)
    colnames(dataset)[in_fname] = longname2fid[colnames(dataset)[in_fname]]
    
    # check if we need to "spread" the table
    if(fname %in% names(longitudinal_csv2id_col)){
      if(length(days_visit_col)!=1){
        print(paste("Error in longitudinal format csv",fname,"days_visit was not detected, skipping"))
        next
      }
      dataset = dataset[,in_fname]
      curr_main_id_col = longitudinal_csv2id_col[fname]
      curr_main_id_col = longname2fid[curr_main_id_col]
      dataset = spread_data_by_visits(dataset,curr_main_id_col,days_visit_col)
    }
    
    if(length(m)==0){
      m = dataset
    }
    else{
      # update the by_col to be in field ids
      by_col = longname2fid[csv_name_main_field[fname]]
      m = merge(m,dataset,by=by_col,suffixes = c("",".yyy"),all=T)
    }
    
    print(paste("After adding data from:",fname,"we now have:",nrow(m),"rows in the merged matrix"))
  }
  # remove undesired colums, merging artifacts, and sort
  to_keep = grepl(FIELD_ID_PREFIX,colnames(m)) &
    (!grepl(".yyy$",colnames(m),perl = T))
  m_cols = sort(colnames(m)[to_keep])
  m = m[,m_cols]
  bucket2merged_data[[bucket_name]][["labelid_data"]] = m
}

# Add the viallabels for the available buckets
for(bucket_name in names(bucket2merged_data)){
  if(!"BICLabelData" %in% names(bucket2data[[bucket_name]])){
    next
  }
  m = bucket2merged_data[[bucket_name]][["labelid_data"]]
  biclabel_data = bucket2data[[bucket_name]][["BICLabelData"]]
  labelid_field = names(ID_COLS)[ID_COLS=="labelid"]
  vialabel_field = names(ID_COLS)[ID_COLS=="viallabel"]
  mapping_info = biclabel_data$filtered_dataset
  mapping_info_dict = biclabel_data$filtered_dict
  colnames(mapping_info) =  tolower(colnames(mapping_info))
  mapping_info = mapping_info[,c("labelid","viallabel","shiptositeid")]
  
  if(!length(unique(mapping_info$viallabel)) == nrow(mapping_info)){
    print(paste("Error in",bucket_name,"vials in bic mapping info matrix are not unique, skipping"))
    next
  }
  
  # add the viallabel description
  merged_dictionary[vialabel_field,"Alias..Field.Name.description."] = 
    mapping_info_dict[
      grepl("viallabel",mapping_info_dict[,"Field.Name"],ignore.case = T),
      grepl("description",colnames(mapping_info_dict),ignore.case = T)]
  
  if(! SITE_FIELD_NAME %in% merged_dictionary[,"Field.Name"]){
    site_field_id = generate_new_field_id()
    names(SITE_FIELD_NAME) = site_field_id
    v = rep("",ncol(merged_dictionary))
    names(v) = colnames(merged_dictionary)
    curr_r_ind = which(grepl(SITE_FIELD_NAME,mapping_info_dict[["Field.Name"]],ignore.case = T))
    if(length(curr_r_ind)!=1){
      print(paste("Error, there is no single",SITE_FIELD_NAME,"field name in the BICLabel csv dictionary"))
    }
    v[1] = site_field_id
    v[2] = paste("BICLabelData",SITE_FIELD_NAME,sep=".")
    v[3] = SITE_FIELD_NAME
    v[4] = "A code for the chemical analysis site (i.e., the omics site)"
    v["Data.Type"] = "varchar"
    v["FormName"] = mapping_info_dict[curr_r_ind,1]
    v["Categorical.Values"] = mapping_info_dict[curr_r_ind,"Categorical.Values"]
    v["Categorical.Definitions"] = mapping_info_dict[curr_r_ind,"Categorical.Definitions"]
    merged_dictionary = rbind(merged_dictionary,v)
    rownames(merged_dictionary) = merged_dictionary[,1]
  }
  
  # create the vial-based table
  newm = merge(m,mapping_info,by.x=labelid_field,
               by.y = "labelid",all=T,suffixes = c("","yyy"))
  newm = newm[,!grepl(".yyy$",colnames(newm),perl = T)]
  newm = newm[!is.na(newm$viallabel),]
  rownames(newm) = mapping_info$viallabel
  colnames(newm)[(ncol(newm)-1):ncol(newm)] = c(vialabel_field,names(SITE_FIELD_NAME))
  # reorganize the field names
  newm = newm[,sort(colnames(newm))]
  
  if(nrow(newm) != length(unique(newm[,vialabel_field]))){
    print(paste("Error in",bucket_name,"vial labels are not unique in the merged data"))
  }
  bucket2merged_data[[bucket_name]][["viallabel_data"]] = newm
}

#########################################################################################
#########################################################################################
print("Running QA to verify the merged datasets, if all are fine no errors are printed below")
print("--------------------------------------------------------------")

# Check the animal-level datasets
print("checking the labelid-level tables; these should contain all animals")
for(bucket_name in names(bucket2merged_data)){
  m = bucket2merged_data[[bucket_name]][[1]]
  for(field in colnames(m)){
    if(field %in% names(ID_COLS)){next}
    fieldid = strsplit(field,split="_")[[1]][1]
    f_visit = NULL
    if(grepl("_",field)){
      f_visit = strsplit(field,split="_")[[1]][2]
    }
    f_info = merged_dictionary[fieldid,]
    f_fullname = f_info["FullName"]
    f_fieldname = f_info["Field.Name"]
    f_fname = gsub(paste0(".",f_fieldname),"",f_fullname)
    if(f_fname==""){next}
    form_main_id = "pid"
    rawdata = bucket2data[[bucket_name]][[f_fname]]$dataset
    if("labelid" %in% colnames(rawdata)){
      form_main_id = "labelid"
    }
    if(is.null(f_visit)){
      rawdata = rawdata[,c(form_main_id,f_fieldname)]
      rawdata = rawdata[!is.na(rawdata[,1]),]
      rawdata = unique(rawdata)
      rownames(rawdata) = as.character(rawdata[,1])
      
      ourdata = unique(m[,longname2fid[c(form_main_id,f_fullname)]])
      ourdata = ourdata[!is.na(ourdata[,1]),]
      rownames(ourdata) = ourdata[,1]
      
      if(!all(rawdata[,1] %in% ourdata[,1])){
        print(paste("Error in bucket:",bucket_name," field:",f_fullname,
                    "using id column:",form_main_id,
                    "not all IDs in rawdata are covered in the merged data"))
        break
      }
      
      # check values
      if(!all(rawdata[,2] == ourdata[rownames(rawdata),2],na.rm = T)){
        print(paste("Error in bucket:",bucket_name," field:",f_fullname,
                    "using id column:",form_main_id,
                    "not all VALUEs in rawdata fit the merged data"))
        break
      }
      # check NAs
      if(!all(is.na(rawdata[,2]) == is.na(ourdata[rownames(rawdata),2]))){
        print(paste("Error in bucket:",bucket_name," field:",f_fullname,
                    "using id column:",form_main_id,
                    "not all NAs indices in rawdata fit the merged data"))
        break
      }
    }
    else{
      rawdata = rawdata[,c(form_main_id,f_fieldname,"days_visit")]
      rawdata = rawdata[!is.na(rawdata[,1]),]
      
      visit_ind = 1
      curr_id = paste(fieldid,"_",visit_ind,sep="")
      ourdata=c()
      while(curr_id %in% colnames(m)){
        currdata = unique(m[,c(longname2fid[form_main_id],curr_id)])
        currdata = currdata[!is.na(currdata[,1]),]
        currdata = cbind(currdata,curr_id)
        colnames(currdata) = c(longname2fid[form_main_id],fieldid,"fullid")
        ourdata = rbind(ourdata,currdata)
        visit_ind = visit_ind+1
        curr_id = paste(fieldid,"_",visit_ind,sep="")
      }
      
      
      if(!all(rawdata[,1] %in% ourdata[,1])){
        print(paste("Error in bucket:",bucket_name," field:",f_fullname,
                    "using id column:",form_main_id,
                    "not all IDs in rawdata are covered in the merged data"))
      }
      
      for(id in unique(rawdata[,1])){
        tmp_our = ourdata[ourdata[,1]==id,]
        tmp_raw = rawdata[rawdata[,1]==id,]
        tmp_raw = tmp_raw[order(tmp_raw[,3]),]
        tmp_our = tmp_our[1:nrow(tmp_raw),]
        if(!all(tmp_our[,1:2]==tmp_raw[,1:2],na.rm = T)){
          print(paste("Error in bucket:",bucket_name," field:",f_fullname,
                      "using id column:",form_main_id,
                      "sample id:",id,
                      "not all VALUEs in rawdata are the same as the merged data"))
          break
        }
        if(!all(is.na(tmp_our[,1:2])==is.na(tmp_raw[,1:2]))){
          print(paste("Error in bucket:",bucket_name," field:",f_fullname,
                      "using id column:",form_main_id,
                      "sample id:",id,
                      "not all NA indices in rawdata are the same as the merged data"))
          break
        }
      }
    }
  }
}

print("checking the vaillabel tables; these are compared to the label-level data")
print("(should be used only if there are no errors printed for the label-level data)")
for(bucket_name in names(bucket2merged_data)){
  m = bucket2merged_data[[bucket_name]][[1]]
  if(!"viallabel_data" %in% names(bucket2merged_data[[bucket_name]])){
    next
  }
  m2 = bucket2merged_data[[bucket_name]][["viallabel_data"]]
  m2 = m2[,intersect(colnames(m),colnames(m2))]
  m = m[,intersect(colnames(m),colnames(m2))]
  labelid_col = longname2fid["labelid"]
  m = m[m[,labelid_col] %in% m2[,labelid_col],]
  viallabel_col = longname2fid["viallabel"]
  m2 = m2[,colnames(m2)!=viallabel_col]
  rownames(m2) = NULL
  m2 = unique(m2)
  if(!all(m==m2,na.rm = T)){
    print(paste("Error in bucket",bucket_name,"vial-level data and labelid-level data do not fit"))
  }
}
print("Done running all tests!")
#########################################################################################
#########################################################################################
print("Saving the results to the bucket:")
print("--------------------------------------------------------------")
currdate = Sys.Date()
OUTBUCKET = paste0(OUTBUCKET,currdate,"/")
print(OUTBUCKET)
system(paste("~/google-cloud-sdk/bin/gsutil mb",OUTBUCKET))

setwd(local_path)
print("writing merged_dictionary.txt")
fwrite(merged_dictionary,file="merged_dictionary.txt",quote = F,sep="\t",row.names = F)
system(paste("~/google-cloud-sdk/bin/gsutil cp", "merged_dictionary.txt",OUTBUCKET))
system(paste("rm","merged_dictionary.txt"))

for(bucket_name in names(bucket2merged_data)){
  for(table_name in names(bucket2merged_data[[bucket_name]])){
    curr_file_name = paste(bucket_name,table_name,sep="_")
    curr_file_name = paste0(curr_file_name,".txt")
    print(paste("writing",curr_file_name))
    curr_dataset = bucket2merged_data[[bucket_name]][[table_name]]
    # Weird bug (fixed on June 2020) - some strings have newlines and it
    # causes issues when printing the tables
    for(j in names(curr_dataset)){
      if(is.character(curr_dataset[[j]])){
        curr_dataset[[j]] = gsub("\n","",curr_dataset[[j]])
      }
    }
    fwrite(curr_dataset,file=curr_file_name,quote = F,sep="\t",row.names = F,col.names=T)
    system(paste("~/google-cloud-sdk/bin/gsutil cp", curr_file_name,OUTBUCKET))
    system(paste("rm",curr_file_name))
  }
}


