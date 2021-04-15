dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/contactmatrix"
out_path="/p/keles/fandingzhou/volumeA/TAD/GM12878"

library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 


ch1_length<-249250621
bin_size<-40000
rep<-c(2,3,4,6)
seqdepth<-c(1,3,5)
nrep<-4
scalesize<-c(2,3,5,10)

n_bin<-floor(ch1_length/bin_size)
cm <- matrix(rep(0,n_bin*n_bin),n_bin, n_bin)

##GM original
for (j in 1:nrep){
  data_j<-read_tsv(paste0(in_path,"/GM",rep[j],"_contactmatrix.tsv"),col_names=TRUE)
  cm <- cm+data_j
}
write.table(cm/4, file=paste0(out_path,"/GM_combined_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE,col.names=F)
#GM spikein
for(k in 1:3){
  cm <- matrix(rep(0,n_bin*n_bin),n_bin, n_bin)
  for(j in 1:nrep){
    data_j<-read_tsv(paste0(in_path,"/GM",rep[j],"_seqx",seqdepth[k],"_contactmatrix.tsv"),col_names=TRUE)
    cm <- cm +data_j
 }
  write.table(cm/4, file=paste0(out_path,"/GM_spike_in_seqx",seqdepth[k],"_combined_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE,col.names=F)
}

# aggregate
for(k in 1:3){
  cm <- matrix(rep(0,n_bin*n_bin),n_bin, n_bin)
  for (j in 1:nrep){
    data_j<-read_tsv(paste0(in_path,"/GM",rep[j],"_contactmatrix.tsv"),col_names=TRUE)
    cm <- cm+data_j
  }
  for(j in 1:nrep){
    data_j<-read_tsv(paste0(in_path,"/GM",rep[j],"_seqx",seqdepth[k],"_contactmatrix.tsv"),col_names=TRUE)
    cm <- cm +data_j
  }
  write.table(cm, file=paste0(out_path,"/GM_spike_in_seqx",seqdepth[k],"_aggregated_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE,col.names=F)
}

# aggregate for scale
cm1 <- matrix(rep(0,n_bin*n_bin),n_bin, n_bin)
for (j in 1:nrep){
  data_j<-read_tsv(paste0(in_path,"/GM",rep[j],"_contactmatrix.tsv"),col_names=TRUE)
  cm1 <- cm1+data_j
}
for(k in 1:3){
  for(ss in 1:4){
  cm=cm1
  for(j in 1:nrep){
    data_j<-read_tsv(paste0(in_path,"/GM",rep[j],"_scale",scalesize[ss],"_seqx",seqdepth[k],"_contactmatrix.tsv"),col_names=TRUE)
    cm <- cm +data_j
  }
  write.table(cm, file=paste0(out_path,"/GM_spike_in_scale",scalesize[ss],"_seqx",seqdepth[k],"_aggregated_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE,col.names=F)
}}


# aggregate for new sim
cm1=cm2 <- matrix(rep(0,n_bin*n_bin),n_bin, n_bin)
for (j in 1:nrep){
  data_j<-read_tsv(paste0(in_path,"/GM",rep[j],"_contactmatrix.tsv"),col_names=TRUE)
  cm2 <- cm2+data_j
}

out_path2="/p/keles/fandingzhou/volumeA/ABcompartment/pca"

for(seed in 1:10){
  for(foldchange in c(2,3,4,6)){
    cm=cm1
    for(j in 1:nrep){
      data_j<-read_tsv(paste0(in_path,"/seed",seed,"_FC",foldchange,"_GM",rep[j],"_contactmatrix.tsv"),col_names=TRUE)
      cm <- cm +data_j
    }
    write.table(cm, file=paste0(out_path2,"/seed",seed,"_GM_FC",foldchange,"_aggregated_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE,col.names=F)
    if(seed>1){
      write.table(cm+cm2, file=paste0(out_path,"/seed",seed,"_GM_FC",foldchange,"_aggregated_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE,col.names=F)
    }
   
  }
  
}

foldchange=2
cm=read.table(paste0(out_path,"/sim3_GM_FC",foldchange,"_aggregated_contactmatrix.tsv"))
#which(colSums(cm)==0)
write.table(cm[1:3024,1:3024],file=paste0(out_path,"/sub_GM_FC",foldchange,"_aggregated_contactmatrix.tsv"), sep="\t",quote=FALSE, row.names=FALSE,col.names=F)