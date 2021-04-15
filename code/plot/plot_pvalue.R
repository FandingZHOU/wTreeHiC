dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
meth=c("multiHiCcompare","diffHic")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data","/p/keles/fandingzhou/volumeA/DCI/diffHic/data")
doc_name<-c("_Multi_glm","_diffHic.tsv")
seqdepth<-c(1,3,5)

library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 

weight_method<-c("xi","pi0")
for (m in 1:2){
for (n in 1:3){
  #HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/binpairs_GMvsGMspikein_seqx",seqdepth[n],doc_name[m]),col_names=TRUE)
  HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
  hist(HiCbinpairs_data$P,main=paste0("p_",meth[m],"_GMvsGMspikein_seqx",seqdepth[n]))
}}

for (m in 1:2){
  for (n in 1:3){
    for(k in 1:2){
    #HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/binpairs_GMvsGMspikein_seqx",seqdepth[n],doc_name[m]),col_names=TRUE)
    HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
    weighted_p<-HiCbinpairs_data$weighted_P
    weighted_p[weighted_p>1]<-1
    hist(weighted_p,xlim = c(0,1),main=paste0("weighed_p_",meth[m],"_GMvsGMspikein_seqx",seqdepth[n],"_",weight_method[k]))
  }
}
}