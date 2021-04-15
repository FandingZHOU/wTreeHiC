dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/contactmatrix"
in_path2="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/filter"
out_path="/p/keles/fandingzhou/volumeA/DCI/simulation/"
#out_path="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/originaltable"
func_path="/p/keles/fandingzhou/volumeA/DCI/code/diffHic"
options("scipen"=100)

library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(parallel)

ch1_length<-249250621
bin_size<-40000
rep<-c(2,3,4,6)
seqdepth<-c(1,3,5)
nrep<-4
scalesize<-c(2,3,5,10)
#prob_FC<-0.1
prob_FC<-0.05
#dispersion=10^(-4)

data_j<-list()
for (j in 1:nrep){
  data_j[[j]]<-read_tsv(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j],"/chr1/s1_training/binPairs/rep",rep[j],"_chr1.binPairs.chr1"),col_names=FALSE)%>% 
    filter(X1 == "chr1", X3 == "chr1")
  data_j[[j+nrep]]<-read_tsv(paste0("/p/keles/scrna-seq/volumeB/freeHiC/A549/rep",j,"/chr1/s1_training/binPairs/rep",j,"_chr1.binPairs.chr1"),col_names=FALSE)%>% filter(X1 == "chr1", X3 == "chr1")
}
merged_matrix<-Reduce(function(x,y) full_join(x,y,by = c("X1","X3","X2","X4")),data_j);
merged_matrix<-merged_matrix %>%
  mutate_all(~replace(.,is.na(.),0))
colnames(merged_matrix)<-c("chrA","region1","chrB","region2",paste0("IF",1:8))

merged_matrix<-merged_matrix[rowSums(merged_matrix[5:12])>2,]
long_range<-which(merged_matrix$region2-merged_matrix$region1>2*bin_size)
ave_matrix_GM<-rowMeans(merged_matrix[,5:8])
#ave_matrix_A<-rowMeans(merged_matrix[,9:12])
#merged_matrix_GM<-merged_matrix[,c(2,4:8)]


#Strategy3
#larger_IF_index<-which(ave_matrix_GM>16*mean(ave_matrix_GM))
#large_IF_index<-which(ave_matrix_GM>mean(ave_matrix_GM))
med_IF_index<-which(ave_matrix_GM>mean(ave_matrix_GM))
small_IF_index<-which(ave_matrix_GM<=mean(ave_matrix_GM))

#better_insert<-intersect(large_IF_index,long_range)
good_insert<-intersect(med_IF_index,long_range)
other_insert<-intersect(small_IF_index,long_range)
#prob_FC1<-prob_FC*length(long_range)/length(better_insert)*0.1
prob_FC2<-prob_FC*length(long_range)/length(other_insert)*0.2
prob_FC3<-prob_FC*length(long_range)/length(good_insert)*0.8

# get_dispersion<-function(x){
#   if(var(x)<mean(x)){
#     return(0)
#   }else{
#     return((var(x)-mean(x))/mean(x)^2)
#   }
# }
get_dispersion<-function(x){
  if(mean(x)==0){
    return(0)
  }else{
    return((var(x)-mean(x))/mean(x)^2)
  }
}
#dispersion[1:length(long_range)]=mean(apply(as.matrix(merged_matrix[long_range,5:8]),1,get_dispersion))
#dispersion[large_IF_index]=mean(apply(as.matrix(merged_matrix[long_range,5:8]),1,get_dispersion))

  #mean(apply(as.matrix(merged_matrix[large_IF_index,5:8]),1,get_dispersion))
small_disp<-mean(apply(as.matrix(merged_matrix[small_IF_index,5:8]),1,get_dispersion))
med_disp<-mean(apply(as.matrix(merged_matrix[med_IF_index,5:8]),1,get_dispersion))


#short_range<-which(merged_matrix$region2-merged_matrix$region1<=2*bin_size)
#length(intersect(which(ave_matrix_GM>8*mean(ave_matrix_GM)),short_range))
nb_sim2<-function(x,foldchange){
  if(x[1]<0.25){
    x[1]=0.25
  }
  if(x[2]!=0){
    return(rnbinom(n=1,size=1/x[3],mu=x[1]*foldchange^x[2]))
  }else{
    return(rnbinom(n=1,size=1/x[3],mu=x[1]))
  }
}

#set.seed(101)
#for(seed in 5:10){
for(seed in 11:20){
  set.seed(seed)
  
  direction<-rep(0,nrow(merged_matrix))
  #direction[better_insert]<-sample(c(0,1,-1),length(better_insert),replace=T,c(1-prob_FC1,prob_FC1*0.1,prob_FC1*0.9))
  direction[other_insert]<-sample(c(0,1),length(other_insert),replace=T,c(1-prob_FC2,prob_FC2))
  direction[good_insert]<-sample(c(0,1,-1),length(good_insert),replace=T,c(1-prob_FC3,prob_FC3*0.38,prob_FC3*0.62))
  
  #mu_disper<-cbind(ave_matrix_GM,direction)
  True_sig1<-which(direction!=0)
  True_sig<-paste0(merged_matrix[True_sig1,]$region1,"_",merged_matrix[True_sig1,]$region2)
  #write(True_sig,paste0(out_path,"True_signal4"))
  write(True_sig,paste0(out_path,"True_signal",seed))
  dispersion =c()
  dispersion[small_IF_index]=small_disp
  dispersion[unique(c(True_sig1,med_IF_index))]=med_disp
  
  for(foldchange in c(2,3,4,6)){
    for(i in 5:8){
      mu_disper<-cbind(merged_matrix[,i],direction,dispersion)
      detectCores(logical=F)
      mc=getOption("mc.cores",50)
      mu_disper<-mclapply(seq_len(nrow(mu_disper)), function(i) as.numeric(mu_disper[i,]),mc.cores=mc)
      #data_sim<-mclapply(ave_matrix_GM,nb_sim,direction,dispersion,foldchange=4,mc.cores=mc)
      data_sim<-mclapply(mu_disper,nb_sim2,foldchange,mc.cores=mc)
      #stopCluster(mc)
      data_sim<-do.call(rbind,data_sim)
      binpairdata<-cbind(merged_matrix[,1:4],data_sim)
      #write_tsv(binpairdata,paste0(out_path,"GM12878/strategy4_FC",foldchange,"_rep",rep[i-4],"_chr1.binPairs.chr1"))
      write_tsv(binpairdata,paste0(out_path,"GM12878/seed",seed,"_FC",foldchange,"_rep",rep[i-4],"_chr1.binPairs.chr1"))
      
      #data_sim<-data.frame(t(sapply(ave_matrix_GM,nb_sim,direction,dispersion,foldchange=4)))
    }
  }
}
  




test.speed<-function(){
  time1<-Sys.time()
  for(i in 1:37438){
    mu<-runif(0,10)
    sd<-runif(1,10)
    rnbinom(1,10000,mu)
    rnorm(1,mu)
  }
  time2<-Sys.time()
  return(time2-time1)
}


#strategy 1
better_insert<-intersect(which(ave_matrix_GM>2*mean(ave_matrix_GM)),long_range)
other_insert<-intersect(which(ave_matrix_GM<=2*mean(ave_matrix_GM)),long_range)
prob_FC1<-prob_FC*length(long_range)/length(better_insert)*0.8
prob_FC2<-prob_FC*length(long_range)/length(other_insert)*0.2
set.seed(101)
direction<-rep(0,nrow(merged_matrix))
direction[better_insert]<-sample(c(0,1,-1),length(better_insert),replace=T,c(1-prob_FC1,prob_FC1*0.4,prob_FC1*0.6))
direction[other_insert]<-sample(c(0,1),length(other_insert),replace=T,c(1-prob_FC2,prob_FC2))
#mu_disper<-cbind(ave_matrix_GM,direction)
True_sig<-which(direction!=0)
True_sig<-paste0(merged_matrix[True_sig,]$region1,"_",merged_matrix[True_sig,]$region2)
write(True_sig,paste0(out_path,"True_signal2"))


#strategy 2
better_insert<-intersect(which(ave_matrix_GM>3*mean(ave_matrix_GM)),long_range)
good_insert<-intersect(which(ave_matrix_GM<=3*mean(ave_matrix_GM)&ave_matrix_GM>0.5*mean(ave_matrix_GM)),long_range)
other_insert<-intersect(which(ave_matrix_GM<=0.5*mean(ave_matrix_GM)),long_range)
prob_FC1<-prob_FC*length(long_range)/length(better_insert)*0.1
prob_FC2<-prob_FC*length(long_range)/length(other_insert)*0.1
prob_FC3<-prob_FC*length(long_range)/length(good_insert)*0.8

set.seed(101)
direction<-rep(0,nrow(merged_matrix))
direction[better_insert]<-sample(c(0,0.5,-0.5),length(better_insert),replace=T,c(1-prob_FC1,prob_FC1*0.1,prob_FC1*0.9))
direction[other_insert]<-sample(c(0,1),length(other_insert),replace=T,c(1-prob_FC2,prob_FC2))
direction[good_insert]<-sample(c(0,1,-1),length(good_insert),replace=T,c(1-prob_FC3,prob_FC3*0.5,prob_FC3*0.5))

#mu_disper<-cbind(ave_matrix_GM,direction)
True_sig<-which(direction!=0)
True_sig<-paste0(merged_matrix[True_sig,]$region1,"_",merged_matrix[True_sig,]$region2)
write(True_sig,paste0(out_path,"True_signal3"))

#Strategy 1-2
nb_sim<-function(x,dispersion,foldchange){
  if(x[1]<0.25){
    x[1]=0.25
  }
  if(x[2]!=0){
    return(rnbinom(n=1,size=1/dispersion,mu=x[1]*foldchange^x[2]))
  }else{
    return(rnbinom(n=1,size=1/dispersion,mu=x[1]))
  }
}
for(foldchange in c(2,4,6)){
  for(i in 5:8){
    mu_disper<-cbind(merged_matrix[,i],direction)
    detectCores(logical=F)
    mc=getOption("mc.cores",50)
    mu_disper<-mclapply(seq_len(nrow(mu_disper)), function(i) as.numeric(mu_disper[i,]),mc.cores=mc)
    #data_sim<-mclapply(ave_matrix_GM,nb_sim,direction,dispersion,foldchange=4,mc.cores=mc)
    data_sim<-mclapply(mu_disper,nb_sim,dispersion,foldchange,mc.cores=mc)
    #stopCluster(mc)
    data_sim<-do.call(rbind,data_sim)
    binpairdata<-cbind(merged_matrix[,1:4],data_sim)
    write_tsv(binpairdata,paste0(out_path,"GM12878/stratege3_FC",foldchange,"_rep",rep[i-4],"_chr1.binPairs.chr1"))
    
    #data_sim<-data.frame(t(sapply(ave_matrix_GM,nb_sim,direction,dispersion,foldchange=4)))
  }
}



