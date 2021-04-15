##Calculating O/E Matrix
dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/contactmatrix"
in_path2="/p/keles/fandingzhou/volumeA/TAD/GM12878"
out_path="/p/keles/fandingzhou/volumeA/ABcompartment/pca"
#out_path="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/originaltable/filter"
func_path="/p/keles/fandingzhou/volumeA/DCI/code/diffHic"
options("scipen"=100)

library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(parallel)
library(bigpca)

ch1_length<-249250621
bin_size<-40000
rep<-c(2,3,4,6)
seqdepth<-c(1,3,5)
nrep<-4
scalesize<-c(2,3,5,10)

n_bin<-floor(ch1_length/bin_size)
delta <- rep(seq_len(n_bin), n_bin) - 
  rep(seq_len(n_bin), each = n_bin)


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
distance<-0:(n_bin-1)

seed=14
foldchange=6
obs_matrix<-as.matrix(read_tsv(paste0(out_path,"/seed",seed,"_GM_FC",foldchange,"_aggregated_contactmatrix.tsv"),col_names=F))
n.cores=getOption("mc.cores",detectCores()-10)
expect_no0<-get_expected_matrix(obs_matrix,n_bin,bin_size,n.cores)
write.table(expect_no0,paste0(out_path,"/seed",seed,"_GM_FC",foldchange,"_expect_no0.tsv"))
EvsO_no0_matrix<-get_OvsE_matrix(obs_matrix,expect_no0$expected,n_bin)
write.table(EvsO_no0_matrix,paste0(out_path,"/seed",seed,"_GM_FC",foldchange,"_eovero_no0.tsv"))
compartment_boundaries<-get_ABcompartment(EvsO_no0_matrix)
write_tsv(compartment_boundaries,paste0(out_path,"/seed",seed,"_GM_FC",foldchange,"_eigen_no0"),col_names= F)


#data_j<-read_tsv(paste0(in_path2,"/GM_combined_contactmatrix.tsv"),col_names=TRUE)
foldchange=6
for(seed in 16:20){
  obs_matrix<-as.matrix(read_tsv(paste0(out_path,"/seed",seed,"_GM_FC",foldchange,"_aggregated_contactmatrix.tsv"),col_names=F))
  mc=getOption("mc.cores",detectCores()-10)
  #a<-unlist(mclapply(distance,exp_ma_func,obs_matrix,mc.cores =mc))
  b<-unlist(mclapply(distance,exp_no0_func,obs_matrix,mc.cores=mc))
  #expect<-data.frame(distance=distance*bin_size,expected=a)
  expect_no0<-data.frame(distance=distance*bin_size,expected=b)
  #write.table(expect,paste0(out_path,"/sim3_GM_FC",foldchange,"_expect.tsv"))
  write.table(expect_no0,paste0(out_path,"/seed",seed,"_GM_FC",foldchange,"_expect_no0.tsv"))
}
for(seed in 19:20){
  obs_matrix<-as.matrix(read_tsv(paste0(out_path,"/seed",seed,"_GM_FC",foldchange,"_aggregated_contactmatrix.tsv"),col_names=F))
  EvsO_no0_matrix<-matrix(rep(0,n_bin*n_bin),n_bin, n_bin)
  b<-read.table(paste0(out_path,"/seed",seed,"_GM_FC",foldchange,"_expect_no0.tsv"))$expected
  for(i in distance){
    ind<-which(abs(delta)==i)
    # if(a[i+1]!=0){
    #  EvsO_matrix[ind]<-obs_matrix[ind]/a[i+1]
    # }
    if(!is.na(b[i+1])){
      EvsO_no0_matrix[ind]<-obs_matrix[ind]/b[i+1]
    }
  }
  
  
  
  # write.table(EvsO_matrix,paste0(out_path,"/sim3_GM_FC",foldchange,"_eovero.tsv"))
  write.table(EvsO_no0_matrix,paste0(out_path,"/seed",seed,"_GM_FC",foldchange,"_eovero_no0.tsv"))
  obs_matrix=EvsO_no0_matrix<-data.frame()
  
  
}

  


for(seed in 19:20){
  EvsO_no0_matrix<-read.table(paste0(out_path,"/seed",seed,"_GM_FC",foldchange,"_eovero_no0.tsv"))
  #S1<-prcomp(EvsO_matrix)
  S2 <- prcomp(EvsO_no0_matrix)
  eigenvec<-scale(as.vector(S2$rotation[,1]))
  
  #write.table(as.vector(S1$rotation[,1]),paste0(out_path,"/sim3_GM_FC",foldchange,"_eigen"))
  write.table(data.frame(X1=0:(n_bin-1),X2=eigenvec),paste0(out_path,"/seed",seed,"_GM_FC",foldchange,"_eigen_no0"),row.names = F)
  rm(S2)
  rm(EvsO_no0_matrix)

}
#for(foldchange in c(3,4,6)){
  #}

