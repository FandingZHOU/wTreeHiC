dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
meth=c("multiHiCcompare","diffHic")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data","/p/keles/fandingzhou/volumeA/DCI/diffHic/data")
doc_name<-c("_Multi_glm","_diffHic.tsv")
seqdepth<-c(1,3,5)

library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(ggplot2)
#accuracy
#detected/real
load("/p/keles/treehic/volumeD/freeHiC/GM12878_A549_TSetKey_trueSpikeinSmooth.RData")
weight_method<-c("xi","pi0")
dci_or_fdr<-"FDR"


for(k in 1:2){
  for (m in 1:2){
    for (n in 1:3){
    #HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/binpairs_GMvsGMspikein_seqx",seqdepth[n],doc_name[m]),col_names=TRUE)
    HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
    threshold<-c(0.001,0.005,0.01,0.05,seq(0.1,0.9,0.1))
    FDR_weight_detected<-c()
    DCI_weight_true<-c()
    for (j in 1:length(threshold)){
      maxpi<-max(HiCbinpairs_data[HiCbinpairs_data$adj_w_p<=threshold[j],]$weighted_P)
      detected_sig_data<-HiCbinpairs_data[HiCbinpairs_data$weighted_P<=maxpi,]
      detected_sig_n<-nrow(detected_sig_data)
      datapairs<-c()
      for (i in 1:detected_sig_n){
        datapairs[i]<-paste0(detected_sig_data[i,1],"_",detected_sig_data[i,2])
      }
      DCI_weight_true[j]<-length(intersect(datapairs,TSetKey))
      FDR_weight_detected[j]<-1-DCI_weight_true[j]/detected_sig_n
    }
    adj_p<-c()
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data$P),])
    for (i in 1:nrow(HiCbinpairs_data)){
      adj_p[i]<-nrow(HiCbinpairs_data)*HiCbinpairs_data$P[i]/i
    }
    FDR_orig_detected<-c()
    DCI_orig_true<-c()
    for (j in 1:length(threshold)){
      maxpi<-max(HiCbinpairs_data[adj_p<=threshold[j],]$P)
      detected_sig_data=HiCbinpairs_data[HiCbinpairs_data$P<=maxpi,]
      datapairs<-c()
      detected_sig_n<-nrow(detected_sig_data)
      for (i in 1:detected_sig_n){
        datapairs[i]<-paste0(detected_sig_data[i,1],"_",detected_sig_data[i,2])
      }
      DCI_orig_true[j]<-length(intersect(datapairs,TSetKey))
      FDR_orig_detected[j]<-1-DCI_orig_true[j]/detected_sig_n
    }
    supp=c(rep("orignal",length(threshold)),rep("weighted",length(threshold)))
    theoretical_fdr=rep(threshold,2)
    if (dci_or_fdr=="DCI"){
      power=c(DCI_orig_true,DCI_weight_true)
      tgg=data.frame(supp,theoretical_fdr,power)
      print(ggplot(tgg, aes(x=theoretical_fdr, y=power, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_theoretical_fdrvs_power_",weight_method[k],"_seqx",seqdepth[n])))
    }else{
      detected_fdr=c(FDR_orig_detected,FDR_weight_detected)
      tgg=data.frame(supp,theoretical_fdr,detected_fdr)
      print(ggplot(tgg, aes(x=theoretical_fdr, y=detected_fdr, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_theoretical_vs_detected_fdr_",weight_method[k],"_seqx",seqdepth[n]))+geom_abline(slope=1, intercept=0, linetype=2,color="grey")+expand_limits(y=c(0,1)))
    }
    
  }
}
}