options("scipen"=100)

dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path2="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/filter"
meth=c("multiHiCcompare","diffHic")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/weighttable","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/weighttable/filter")
#doc_name<-c("_Multi_glm","_diffHic.tsv")
seqdepth<-c(1,3,5)
weight_method<-c("xi","pi0")
TSetKey<-scan(paste0(in_path2,"/TSetKey"),what = character())

library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(ggplot2)
#accuracy
#detected/real

for(k in 1:2){
  #for (m in 1:2){
    for (n in 1:3){
      #HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/binpairs_GMvsGMspikein_seqx",seqdepth[n],doc_name[m]),col_names=TRUE)
      HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
      # HiCbinpairs_data<-results.r[o.r,]
      #TSet<-scan(paste0(in_path2,"/_GMvsvsGMspikein_seqx",seqdepth[n],"_TSet"),what = character())
      #smoothTSet<-scan(paste0(in_path2,"/_GMvsvsGMspikein_seqx",seqdepth[n],"_smoothTSet"),what = character())
      
      topN<-100000
      interval_thre<-500
      thre<-seq(0,topN,interval_thre)
      accuracy_original<-get_accuracy(HiCbinpairs_data,add_weight=FALSE,topN,interval_thre=500,TSetKey)
      accuracy_weighted<-get_accuracy(HiCbinpairs_data,add_weight=TRUE,topN,interval_thre=500,TSetKey)
      
      supp=c(rep("orignal",length(thre)-1),rep("weighted",length(thre)-1))
      top_N_pvalues=rep(thre[2:(length(thre))],2)
      Accuracy=c(accuracy_original,accuracy_weighted)
      tgg=data.frame(supp,top_N_pvalues,Accuracy)
      #ggplot(tgg, aes(x=threshold, y=Accuracy, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_power_",weight_method[k],"_small"))
      #ggplot(tgg, aes(x=threshold, y=Accuracy, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_power_",weight_method[k],"_large"))
      print(ggplot(tgg, aes(x= top_N_pvalues, y=Accuracy, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_accuracy_",weight_method[k],"_seqx",seqdepth[n],"_TOP",topN)))
      
    }
  #}
}


