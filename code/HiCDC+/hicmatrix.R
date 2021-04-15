dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path="/p/keles/fandingzhou/volumeA/DCI/HiCDCPlus/data/hicpromatrix"
out_path="/p/keles/fandingzhou/volumeA/DCI/HiCDCPlus/data/original/test"
options("scipen"=100)

library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 


ch1_length<-249250621
bin_size<-40000
rep<-c(2,3,4,6)
seqdepth<-c(1,3,5)
nrep<-4
scalesize<-c(2,3,5,10)
n_bin<-floor(ch1_length/bin_size)


#bedfile
all.regions<-data.frame(cbind(rep("chr1",n_bin+1),seq(1,n_bin*bin_size+1,bin_size),c(seq(bin_size,n_bin*bin_size,bin_size),ch1_length),1:(n_bin+1)))
data.table::fwrite(all.regions,file=paste0(in_path,"/bedfile"),sep = "\t",col.names = F)

#matrix
for (j in 1:4){
  binpair_data<-read_tsv(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j],"/chr1/s1_training/binPairs/rep",rep[j],"_chr1.binPairs.chr1"),col_names=FALSE)%>% filter(X1 == "chr1", X3 == "chr1")
  binpair_data<-binpair_data[,c(2,4,5)]
  binpair_data[,1]<-(binpair_data[,1]+bin_size/2)/bin_size
  binpair_data[,2]<-(binpair_data[,2]+bin_size/2)/bin_size
  data.table::fwrite(binpair_data,paste0(in_path,"/GM_original_rep",rep[j]),sep = "\t",col.names = F)
}

for(k in 1:3){
  for(j in 1:4){
    binpair_data<-read_tsv(paste0("/p/keles/treehic/volumeD/freeHiC/GM12878/rep",rep[j],"/chr1/spikeIn3_40kb_smooth_seqx",seqdepth[k],"/simuProcess/s3_binPairs/rep",rep[j],"_chr1.binPairs"),col_names=FALSE)%>% filter(X1 == "chr1", X3 == "chr1")
    binpair_data<-binpair_data[,c(2,4,5)]
    binpair_data[,1]<-(binpair_data[,1]+bin_size/2)/bin_size
    binpair_data[,2]<-(binpair_data[,2]+bin_size/2)/bin_size
    data.table::fwrite(binpair_data,paste0(in_path,"/GM_spikein_rep",rep[j],"_seqx",seqdepth[k]),sep = "\t",col.names = F)
  }
}

in_path1="/p/keles/fandingzhou/volumeA/DCI/simulation/"
#new sim
for(seed in 11:20){
  for(foldchange in c(2,3,4,6)){
    for(j in 1:4){
      binpair_data<-read_tsv(paste0(in_path1,"GM12878/seed",seed,"_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"))
      binpair_data<-binpair_data[,c(2,4,5)]
      binpair_data[,1]<-(binpair_data[,1]+bin_size/2)/bin_size
      binpair_data[,2]<-(binpair_data[,2]+bin_size/2)/bin_size
      data.table::fwrite(binpair_data,paste0(in_path,"/simulation/seed",seed,"_FC",foldchange,"_GM",rep[j]),sep = "\t",col.names = F)
    }
  }
}

