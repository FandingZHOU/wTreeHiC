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
multi_ContactMatrices<-function(input_path,out_path,bin_size,ch_length=NULL,n.cores=4,ngroup=2){
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
  if(ngroup==2){
    for(j in 1:(nrep*2)){
      if(j<=nrep){
        write.table(cm[[j]], file=paste0(out_path,"/contactmatrix/G1_cm",j,".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
      }else{
        write.table(cm[[j]], file=paste0(out_path,"/contactmatrix/G2_cm",j-nrep,".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
      }
    }
    
  }
  ##just for simulated data if we do not want to recalculate cms for original GM data 
  if(ngroup==1){
    for(j in 1:(nrep*2)){
      write.table(cm[[j]], file=paste0(out_path,"/contactmatrix/G2_cm",j,".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
    }
  }
   return(cm)
}

merge_ContactMatrices<-function(cm,out_path,ngroup){
  cm1=cm[[1]]
  ncm=length(cm)
  if(ngroup==1){
    if(ncm>1){
      for(j in 2:ncm){
        cm1 <- cm1+cm[[j]]
      }
    }
    write.table(cm1, file=paste0(out_path,"/aggregated_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE,col.names=F)
  }
  if(ngroup==2){
    cm2=cm[[ncm/2+1]]
    if(ncm>=4){
      for(j in 2:(ncm/2)){
        cm1 <- cm1+cm[[j]]
      }
      for(j in (ncm/2+1):ncm){
        cm2 <- cm2+cm[[j]]
      }
      write.table(cm1, file=paste0(out_path,"/G1_aggregated_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE,col.names=F)
      write.table(cm2, file=paste0(out_path,"/G2_aggregated_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE,col.names=F)
    }
  }
  return(out_path)
}

