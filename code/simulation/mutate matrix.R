dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
out_path="/p/keles/fandingzhou/volumeA/DCI/simulation/"
out_path2="/p/keles/fandingzhou/volumeA/DCI/juicebox/"
options("scipen"=100)

library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(parallel)
rep<-c(2,3,4,6)
for(j in 1:4){
  data_j<-read_tsv(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j],"/chr1/s1_training/binPairs/rep",rep[j],"_chr1.binPairs.chr1"),col_names=FALSE)%>% 
    filter(X1 == "chr1", X3 == "chr1")
  dummy1<-rep(0,nrow(data_j))
  dummy2<-rep(1,nrow(data_j))
  data_j<-data.frame(dummy1,data_j$X1,data_j$X2,dummy1,dummy2,data_j$X3,data_j$X4,dummy2,data_j$X5)
  write_delim(data_j,paste0(out_path2,"GM12878_rep",rep[j],"_chr1.binPairs.chr1"),col_names = F)
}
for(foldchange in c(2,3,4,6)){
  for(j in 1:4){
    data_j<-read_tsv(paste0(out_path,"GM12878/strategy5_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"))
    dummy1<-rep(0,nrow(data_j))
    dummy2<-rep(1,nrow(data_j))
    data_j<-data.frame(dummy1,data_j[,1],data_j[,2],dummy1,dummy2,data_j[,3],data_j[,4],dummy2,data_j[,5])
    write_delim(data_j,paste0(out_path2,"GM12878_strategy5_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"),col_names = F)
  }
}
