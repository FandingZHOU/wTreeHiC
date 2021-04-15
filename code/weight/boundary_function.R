ABcompartment_boundaries<-function(data_origin,bin_size,nbin){
  boundary<-c(0)
  for (i in 2:(nbin-1)){
    if (i*bin_size-boundary[length(boundary)]>3*bin_size){
      if (data_origin[i+1,2]>0&&data_origin[i,2]>=0&&data_origin[i-1,2]<=0){
        boundary<-c(boundary,i*bin_size)
      }
      if (data_origin[i+1,2]<0&&data_origin[i,2]<=0&&data_origin[i-1,2]>=0){
        boundary<-c(boundary,i*bin_size)
      }
    }
  }
  return(boundary)
}
conbined_boundaries<-function(data_origin1,data_origin2,bin_size,nbin){
  boundary1<-ABcompartment_boundaries(data_origin1,bin_size,nbin)
  boundary2<-ABcompartment_boundaries(data_origin2,bin_size,nbin)
  boundary_combined<-sort(c(boundary1,boundary2))
  boundary_combined_nodup<-unique(boundary_combined)
  boundary_combined_merged<-c()
  n_nodup<-length(boundary_combined_nodup)
  i=1
  while(i <n_nodup){
    first=boundary_combined_nodup[i]
    second=boundary_combined_nodup[i+1]
    if(second-first<=2*bin_size){
      boundary_combined_merged<-c(boundary_combined_merged,mean(first,second))
      i=i+2
    }
    else{
      boundary_combined_merged<-c(boundary_combined_merged,boundary_combined_nodup[i])
      i=i+1
    }
  }
  return(boundary_combined_merged)
}

