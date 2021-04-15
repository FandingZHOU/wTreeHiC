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
  return(data.frame(counts))
}
multi_ContactMatrices<-function(input_path,out_path,bin_size,ch_length=NULL,n.cores=4){
  sample_list <- list()
  nbin=c()
  for(j in 1:length(inputfile)){
    sample_list[[j]]<-read_tsv(input_path[j])
    nbin[j]=(max(sample_list[[j]][,4])+bin_size/2)/bin_size
  }
  if(is.null(ch_length)){
    ch_length<-max(nbin)*bin_size
  }
  cm<-mclapply(sample_list,contactMatrix,bin_size,ch_length,mc.cores=n.cores)
  if(!dir.exists(paste0(out_path,"/contactmatrix/"))){
    dir.create(paste0(out_path,"/contactmatrix/"))
  }
  nrep<-length(inputfile)/2
  for(j in 1:(nrep*2)){
    if(j<=nrep){
      write.table(cm[[j]], file=paste0(out_path,"/contactmatrix/G1_cm",j,".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
    }else{
      write.table(cm[[j]], file=paste0(out_path,"/contactmatrix/G2_cm",j-nrep,".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
    }
  }
  return(cm)

}
