##Calculating O/E Matrix

exp_ma_func<-function(dist,obs_matrix){
  ind<-which(abs(delta)==dist)
  return(mean(obs_matrix[ind]))
}
exp_no0_func<-function(dist,obs_matrix){
  ind<-which(abs(delta)==dist)
  no0_vec<-obs_matrix[ind]
  no0_vec<-no0_vec[no0_vec!=0]
  return(mean(no0_vec))
}

##get expected matrix
get_expected_matrix<-function(obs_matrix,n_bin,bin_size,n.cores){
  distance<-0:(n_bin-1)
  b<-unlist(mclapply(distance,exp_no0_func,obs_matrix,mc.cores=n.cores))
  expect_no0<-data.frame(distance=distance*bin_size,expected=b)
  #a<-unlist(mclapply(distance,exp_ma_func,obs_matrix,mc.cores =mc))
  #expect<-data.frame(distance=distance*bin_size,expected=a)
  return(expect_no0)
}

#calculate OvsE
get_OvsE_matrix<-function(obs_matrix,expect_no0,n_bin){
  delta <- rep(seq_len(n_bin), n_bin) - 
    rep(seq_len(n_bin), each = n_bin)
  distance<-0:(n_bin-1)
  EvsO_no0_matrix<-matrix(rep(0,n_bin*n_bin),n_bin, n_bin)
  
 for(i in distance){
    ind<-which(abs(delta)==i)
    # if(a[i+1]!=0){
    #  EvsO_matrix[ind]<-obs_matrix[ind]/a[i+1]
    # }
    if(!is.na(expect_no0[i+1])){
      EvsO_no0_matrix[ind]<-obs_matrix[ind]/expect_no0[i+1]
    }
  }
  return(EvsO_no0_matrix)
}

get_ABcompartment<-function(EvsO_no0_matrix){
  S2 <- prcomp(EvsO_no0_matrix)
  eigenvec<-scale(as.vector(S2$rotation[,1]))
  return(data.frame(X1=0:(n_bin-1),X2=eigenvec))
}

