contactMatrix<-function(data_j,bin_size,ch_length){
  n_bin<-floor(ch_length/bin_size)
  counts <- matrix(rep(0,n_bin*n_bin),n_bin, n_bin)#empty matrix
  data_j_index<-cbind(data_j[,2],data_j[,4])
  data_j_index<-(data_j_index+bin_size/2)/bin_size
  data_j_index<-cbind(data_j_index,data_j[,5])
  Ni<-nrow(data_j)
  for (i in 1:Ni){
    counts[data_j_index[i,1],data_j_index[i,2]]<-data_j_index[i,3]
    counts[data_j_index[i,2],data_j_index[i,1]]<-data_j_index[i,3]
  }
  return(counts)
}
