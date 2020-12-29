dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
out_path="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/contactmatrix"
func_path="/p/keles/fandingzhou/volumeA/DCI/code/diffHic"
library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 

source(paste0(func_path,"/contactmatrix.R"))

ch_length<-249250621
bin_size<-40000
rep<-c(2,3,4,6)
seqdepth<-c(1,3,5)
scalesize<-c(2,3,5,10)

#GMoriginal
for(j in 1:4){
  data_j<-read_tsv(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j],"/chr1/s1_training/binPairs/rep",rep[j],"_chr1.binPairs.chr1"),col_names=FALSE)
  cm<-contactMatrix(data_j,bin_size,ch_length)
  write.table(cm, file=paste0(out_path,"/GM",rep[j],"_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
}

#GMspikein_original
for(k in 1:3){
  for(j in 1:4){
    data_j<-read_tsv(paste0("/p/keles/treehic/volumeD/freeHiC/GM12878/rep",rep[j],"/chr1/spikeIn3_40kb_smooth_seqx",seqdepth[k],"/simuProcess/s3_binPairs/rep",rep[j],"_chr1.binPairs"),col_names=FALSE) %>% filter(X1 == "chr1", X3 == "chr1")
    cm<-contactMatrix(data_j,bin_size,ch_length)
    write.table(cm, file=paste0(out_path,"/GM",rep[j],"_seqx",seqdepth[k],"_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
  }
}

#GMspikein_scale
for(ss in 1:4){
for(k in 1:3){
  for(j in 1:4){
    data_j<-read_tsv(paste0("/p/keles/treehic/volumeD/freeHiC/GM12878/rep",rep[j],"/chr1/spikeIn3_40kb_scale",scalesize[ss],"_seqx",seqdepth[k],"/simuProcess/s3_binPairs/rep",rep[j],"_chr1.binPairs"),col_names=FALSE) %>% filter(X1 == "chr1", X3 == "chr1")
    cm<-contactMatrix(data_j,bin_size,ch_length)
    write.table(cm, file=paste0(out_path,"/GM",rep[j],"_scale",scalesize[ss],"_seqx",seqdepth[k],"_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
  }
}
}

# Aoriginal
for(j in 1:4){
  data_j<-read_tsv(paste0("/p/keles/scrna-seq/volumeB/freeHiC/A549/rep",j,"/chr1/s1_training/binPairs/rep",j,"_chr1.binPairs.chr1"),col_names=FALSE)%>% filter(X1 == "chr1", X3 == "chr1")
  cm<-contactMatrix(data_j,bin_size,ch_length)
  write.table(cm, file=paste0(out_path,"/A",rep[j],"_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
}


