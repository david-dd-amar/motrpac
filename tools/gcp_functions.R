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