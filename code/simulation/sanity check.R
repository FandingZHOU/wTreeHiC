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

#MA plot 4reps
par(mfrow=c(2,2))
for(j in 1:nrep){
  data_j<-list()
  data_j[[1]]<-read_tsv(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j],"/chr1/s1_training/binPairs/rep",rep[j],"_chr1.binPairs.chr1"),col_names=FALSE)%>% 
    filter(X1 == "chr1", X3 == "chr1")
  #data_j[[2]]<-read_tsv(paste0(out_path,"GM12878/less_truesig_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"))
  data_j[[2]]<-read_tsv(paste0(out_path,"GM12878/seed",seed,"_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"))
  #data_j[[2]]<-read_tsv(paste0("/p/keles/scrna-seq/volumeB/freeHiC/A549/rep",j,"/chr1/s1_training/binPairs/rep",j,"_chr1.binPairs.chr1"),col_names=FALSE)%>% filter(X1 == "chr1", X3 == "chr1")
  merged_matrix<-Reduce(function(x,y) full_join(x,y,by = c("X1","X3","X2","X4")),data_j);
  merged_matrix<-merged_matrix %>%
    mutate_all(~replace(.,is.na(.),0))
  colnames(merged_matrix)<-c("chrA","region1","chrB","region2",paste0("IF",1:2))
  samp = sample(1:nrow(merged_matrix),nrow(merged_matrix)/8)
  before = log(merged_matrix[[paste0("IF",1)]]+1)
  after = log(merged_matrix[[paste0("IF",2)]]+1)
  plot((before[samp]+after[samp])/2,after[samp]-before[samp],main=paste0("rep",rep[j]),pch = 46,xlab="log(original count+1)/2+log(FC count+1)/2",ylab="log(FC count+1)-log(original count+1)",ylim = c(-6,6))
  #plot((before[!trueindex]+after[!trueindex])/2,after[!trueindex]-before[!trueindex],main=paste0("rep",rep[j],"_FC",foldchange),pch = 46)
  #points((before[trueindex]+after[trueindex])/2,after[trueindex]-before[trueindex],col="red",pch = 46)
}




par(mfrow=c(2,2))
for(j in 1:nrep){
  data_j<-list()
  data_j[[1]]<-read_tsv(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j],"/chr1/s1_training/binPairs/rep",rep[j],"_chr1.binPairs.chr1"),col_names=FALSE)%>% 
    filter(X1 == "chr1", X3 == "chr1")
  #data_j[[2]]<-read_tsv(paste0(out_path,"GM12878/less_truesig_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"))
  #data_j[[2]]<-read_tsv(paste0(out_path,"GM12878/stratege3_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"))
  data_j[[2]]<-read_tsv(paste0("/p/keles/scrna-seq/volumeB/freeHiC/A549/rep",j,"/chr1/s1_training/binPairs/rep",j,"_chr1.binPairs.chr1"),col_names=FALSE)%>% filter(X1 == "chr1", X3 == "chr1")
  merged_matrix<-Reduce(function(x,y) full_join(x,y,by = c("X1","X3","X2","X4")),data_j);
  merged_matrix<-merged_matrix %>%
    mutate_all(~replace(.,is.na(.),0))
  colnames(merged_matrix)<-c("chrA","region1","chrB","region2",paste0("IF",1:2))
  samp = sample(1:nrow(merged_matrix),nrow(merged_matrix)/8)
  before = log(merged_matrix[[paste0("IF",1)]]+1)
  after = log(merged_matrix[[paste0("IF",2)]]+1)
  plot((before[samp]+after[samp])/2,after[samp]-before[samp],main=paste0("rep",rep[j]),pch = 46,xlab="log(original count+1)/2+log(FC count+1)/2",ylab="log(FC count+1)-log(original count+1)",ylim = c(-6,6))
  #plot((before[!trueindex]+after[!trueindex])/2,after[!trueindex]-before[!trueindex],main=paste0("rep",rep[j],"_FC",foldchange),pch = 46)
  #points((before[trueindex]+after[trueindex])/2,after[trueindex]-before[trueindex],col="red",pch = 46)
}

