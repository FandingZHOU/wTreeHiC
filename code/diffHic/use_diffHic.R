dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/contactmatrix"
in_path2="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/filter"
out_path="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/originaltable/filter"
#out_path="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/originaltable"
func_path="/p/keles/fandingzhou/volumeA/DCI/code/diffHic"
#load("/p/keles/treehic/volumeD/freeHiC/GM12878_A549_dataSign_trueSpikein.RData")
#TSetKey<-paste0(dataSign[,2],"_",dataSign[,4])
options("scipen"=100)
TSetKey<-scan(paste0(in_path2,"/TSetKey"),what = character())


library(BiocGenerics,lib.loc= dict_path)
library(Biobase,lib.loc= dict_path)
library(S4Vectors,lib.loc= dict_path)
library(IRanges,lib.loc= dict_path)
library(BiocParallel,lib.loc = dict_path)#
library(DelayedArray,lib.loc= dict_path)
library(SummarizedExperiment,lib.loc= dict_path)
library(InteractionSet,lib.loc= dict_path)
library(GenomeInfoDb,lib.loc= dict_path)
library(GenomicRanges,lib.loc= dict_path)

library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(limma,lib.loc = dict_path)
library(edgeR,lib.loc = dict_path)# 
library(diffHic,lib.loc = dict_path)
library(csaw,lib.loc = dict_path)
library(statmod,lib.loc = dict_path)
library(ggplot2,lib.loc = dict_path)

#source(paste0(func_path,"/diffHic_table.R"))
#source(paste0(func_path,"/diffHic_table_filter.R"))
#source(paste0(func_path,"/diffHic_MA_plot.R"))

ch1_length<-249250621
bin_size<-40000
rep<-c(2,3,4,6)
seqdepth<-c(1,3,5)
nrep<-4
scalesize<-c(2,3,5,10)

n_bin<-floor(ch1_length/bin_size)
all.regions <- GRanges("chr1", IRanges(0:(n_bin-1)*bin_size+1, 1:n_bin*bin_size))
row.indices <- 1:n_bin
col.indices <- 1:n_bin
row.regions <- all.regions[row.indices]
col.regions <- all.regions[col.indices]

cm<-list()

for (j in 1:nrep){
  data_j<-read_tsv(paste0(paste0(in_path,"/GM",rep[j],"_contactmatrix.tsv")),col_names=TRUE)
  cm[[j]] <- ContactMatrix(as.matrix(data_j), row.regions, col.regions)
}

# GM vs GM spikein filter as multiHiCcompare
for(k in 1:3){
  #transform to contactmatrix dataset
  for(j in 1:nrep){
    data_j<-read_tsv(paste0(in_path,"/GM",rep[j],"_seqx",seqdepth[k],"_contactmatrix.tsv"),col_names=TRUE)
    cm[[j+nrep]] <- ContactMatrix(as.matrix(data_j), row.regions, col.regions)
  }
  #diffHic Pipeline
  data <- mergeCMs(cm[[1]], cm[[2]],cm[[3]],cm[[4]],cm[[5]],cm[[6]],cm[[7]],cm[[8]])
  binpairs<-scan(paste0(in_path2,"/_GMvsvsGMspikein_seqx",seqdepth[k],"_keep"),what = character())
  result<-diffHic_table_filter(data,binpairs,nrep,bin_size)
  write.table(result, file=paste0(out_path,"/binpairs_GMvsGMspikein_seqx",seqdepth[k],"_diffHic.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
}



# GM vs GM spikein
for(k in 1:3){
  #transform to contactmatrix dataset
  for(j in 1:nrep){
    data_j<-read_tsv(paste0(in_path,"/GM",rep[j],"_seqx",seqdepth[k],"_contactmatrix.tsv"),col_names=TRUE)
    cm[[j+nrep]] <- ContactMatrix(as.matrix(data_j), row.regions, col.regions)
  }
  #diffHic Pipeline
  data <- mergeCMs(cm[[1]], cm[[2]],cm[[3]],cm[[4]],cm[[5]],cm[[6]],cm[[7]],cm[[8]])
  result<-diffHic_table(data,nrep,bin_size)
  write.table(result, file=paste0(out_path,"/binpairs_GMvsGMspikein_seqx",seqdepth[k],"_diffHic.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
}







# GM vs GM spikein scaled & filter as multiHiCcompare
for(ss in 1:4){
  for(k in 1:3){
    #transform to contactmatrix dataset
    for(j in 1:nrep){
      data_j<-read_tsv(paste0(in_path,"/GM",rep[j],"_scale",scalesize[ss],"_seqx",seqdepth[k],"_contactmatrix.tsv"),col_names=TRUE)
      cm[[j+nrep]] <- ContactMatrix(as.matrix(data_j), row.regions, col.regions)
    }
    #diffHic Pipeline
    data <- mergeCMs(cm[[1]], cm[[2]],cm[[3]],cm[[4]],cm[[5]],cm[[6]],cm[[7]],cm[[8]])
    binpairs<-scan(paste0(in_path2,"/_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[k],"_keep"),what = character())
    binpairs_tset<-scan(paste0(in_path2,"/_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[k],"_TSet"),what = character())
    result<-diffHic_table_filter(data,binpairs,nrep,bin_size)
    #result1<-diffHic_table(data,binpairs,nrep,bin_size,2)
    write.table(result, file=paste0(out_path,"/binpairs_GMvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[k],"_diffHic.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
  }
}

#check prefiltering
ss=2
k=2
threshold<-c(0.001,0.01,0.05,0.1)
for(j in 1:nrep){
  data_j<-read_tsv(paste0(in_path,"/GM",rep[j],"_scale",scalesize[ss],"_seqx",seqdepth[k],"_contactmatrix.tsv"),col_names=TRUE)
  cm[[j+nrep]] <- ContactMatrix(as.matrix(data_j), row.regions, col.regions)
}
#diffHic Pipeline
data <- mergeCMs(cm[[1]], cm[[2]],cm[[3]],cm[[4]],cm[[5]],cm[[6]],cm[[7]],cm[[8]])
for (i in c(3,4,5,6,7)){
  HiCbinpairs<-diffHic_table(data,nrep,bin_size,i)
  HiCbinpairs_data<-cbind(HiCbinpairs$end2-bin_size/2,HiCbinpairs$end1-bin_size/2,HiCbinpairs$F,HiCbinpairs$PValue,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
  colnames(HiCbinpairs_data)<-c("bin_1","bin_2","F","P","Group","weight")
  HiCbinpairs_data<-data.frame(HiCbinpairs_data)
  fdr<-get_fdr_dci(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="FDR")
  dci<-get_fdr_dci(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="DCI")
  print(fdr)
  print(dci)
}



#sanity check
diffHic_MA_box_plot(data,nrep,bin_size,1,TSetKey,seqdepth[k])

for(ss in 1:4){
  for(j in 1:nrep){
    data_j<-read_tsv(paste0(in_path,"/GM",rep[j],"_scale",scalesize[ss],"_seqx",seqdepth[k],"_contactmatrix.tsv"),col_names=TRUE)
    cm[[j+nrep]] <- ContactMatrix(as.matrix(data_j), row.regions, col.regions)
  }
  #diffHic Pipeline
  data <- mergeCMs(cm[[1]], cm[[2]],cm[[3]],cm[[4]],cm[[5]],cm[[6]],cm[[7]],cm[[8]])
  diffHic_MA_box_plot(data,nrep,bin_size,scalesize[ss],TSetKey,seqdepth[k])
  
  #diffHic_MA_plot(data,nrep,bin_size,scalesize[ss],TSetKey,seqdepth[k])
}

#GM vs A
for(j in 1:nrep){
  data_j<-read_tsv(paste0(in_path,"/A",rep[j],"_contactmatrix.tsv"),col_names=TRUE)
  cm[[j+nrep]] <- ContactMatrix(as.matrix(data_j), row.regions, col.regions)
}
data <- mergeCMs(cm[[1]], cm[[2]],cm[[3]],cm[[4]],cm[[5]],cm[[6]],cm[[7]],cm[[8]])
result<-diffHic_table(data,nrep,bin_size)
write.table(result, file=paste0(out_path,"binpairs_GMvsA_diffHic.tsv"), sep="\t",quote=FALSE, row.names=FALSE)

# GM vs GM
nrep=2
data <- mergeCMs(cm[[1]], cm[[2]],cm[[3]],cm[[4]])
result<-diffHic_table(data,nrep,bin_size)
write.table(result, file=paste0(out_path,"/binpairs_GMvsGM12_diffHic.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
data <- mergeCMs(cm[[1]], cm[[3]],cm[[2]],cm[[4]])
result<-diffHic_table(data,nrep,bin_size)
write.table(result, file=paste0(out_path,"/binpairs_GMvsGM13_diffHic.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
data <- mergeCMs(cm[[1]], cm[[4]],cm[[2]],cm[[3]])
result<-diffHic_table(data,nrep,bin_size)
write.table(result, file=paste0(out_path,"/binpairs_GMvsGM14_diffHic.tsv"), sep="\t",quote=FALSE, row.names=FALSE)

