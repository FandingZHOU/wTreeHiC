dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
meth=c("multiHiCcompare","diffHic")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data","/p/keles/fandingzhou/volumeA/DCI/diffHic/data")
func_path="/p/keles/fandingzhou/volumeA/DCI/code/weight/"
group<-c("A","GM12","GM13","GM14")
seqdepth<-c(1,3,5)
bin_size<-40000
weight_method<-c("xi","pi0")

library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(ggplot2)
load("/p/keles/treehic/volumeD/freeHiC/GM12878_A549_TSetKey_trueSpikeinSmooth.RData")

dci_or_fdr="DCI"
for(k in 1:2){
for(m in 1:2){
  for(n in 1:3){
    HiCbinpairs_data<-original_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
    threshold<-c(0.001,0.005,0.01,0.05,0.1)
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
    
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data$P),])
    adj_p<-c()
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
    
    FDR_tree_detected<-c()
    DCI_tree_true<-c()
    for (j in 1:length(threshold)){
      #original group
      HiCbinpairs_data<-original_data
      HiCbinpairs_data$Group2<-HiCbinpairs_data$rowindex*1000+HiCbinpairs_data$colindex
      
      #L1
      HiCbinpairs_data<-tree_filter_group(HiCbinpairs_data,"Group","P",threshold[j])
      
      #L2
      HiCbinpairs_data<-tree_filter_group(HiCbinpairs_data,"Group2","P",threshold[j])
      
      #L3
      HiCbinpairs_data<-HiCbinpairs_data[order(HiCbinpairs_data$P),]
      pL3<-HiCbinpairs_data$P
      adj_pL3<-c()
      for (i in 1:length(pL3)){
        adj_pL3[i]<-length(pL3)*pL3[i]/i
      }
      maxpi<-max(pL3[adj_pL3<=threshold[j]])
      sig_pair=HiCbinpairs_data[pL3<=maxpi,]
      
      #Compare
      datapairs<-c()
      detected_sig_n<-nrow(sig_pair)
      for (i in 1:detected_sig_n){
        datapairs[i]<-paste0(sig_pair[i,1],"_",sig_pair[i,2])
      }
      DCI_tree_true[j]<-length(intersect(datapairs,TSetKey))
      FDR_tree_detected[j]<-1-DCI_tree_true[j]/detected_sig_n
    }
    
    FDR_wtree_detected<-c()
    DCI_wtree_true<-c()
    for (j in 1:length(threshold)){
      #original group
      HiCbinpairs_data<-original_data
      HiCbinpairs_data$Group2<-HiCbinpairs_data$rowindex*1000+HiCbinpairs_data$colindex
      
      #L1
      HiCbinpairs_data<-tree_filter_group(HiCbinpairs_data,"Group","weighted_P",threshold[j])
      
      #L2
      HiCbinpairs_data<-tree_filter_group(HiCbinpairs_data,"Group2","weighted_P",threshold[j])
      
      #L3
      HiCbinpairs_data<-HiCbinpairs_data[order(HiCbinpairs_data$weighted_P),]
      pL3<-HiCbinpairs_data$weighted_P
      adj_pL3<-c()
      for (i in 1:length(pL3)){
        adj_pL3[i]<-length(pL3)*pL3[i]/i
      }
      maxpi<-max(pL3[adj_pL3<=threshold[j]])
      sig_pair=HiCbinpairs_data[pL3<=maxpi,]
      
      #Compare
      datapairs<-c()
      detected_sig_n<-nrow(sig_pair)
      for (i in 1:detected_sig_n){
        datapairs[i]<-paste0(sig_pair[i,1],"_",sig_pair[i,2])
      }
      DCI_wtree_true[j]<-length(intersect(datapairs,TSetKey))
      FDR_wtree_detected[j]<-1-DCI_wtree_true[j]/detected_sig_n
    }
    
    supp=c(rep("orignal",5),rep("weight",5),rep("tree",5),rep("wtree",5))
    theoretical_fdr=rep(threshold,4)
    if (dci_or_fdr=="DCI"){
      true_positive=c(DCI_orig_true,DCI_weight_true,DCI_tree_true,DCI_wtree_true)
      tgg=data.frame(supp,theoretical_fdr,true_positive)
      print(ggplot(tgg, aes(x=theoretical_fdr, y=true_positive, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_theoretical_fdr_vs_tree_tp_seqx",seqdepth[n],"_",weight_method[k])))
    }else{
      detected_fdr=c(FDR_orig_detected,FDR_weight_detected,FDR_tree_detected,FDR_wtree_detected)
      tgg=data.frame(supp,theoretical_fdr,detected_fdr)
      print(ggplot(tgg, aes(x=theoretical_fdr, y=detected_fdr, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_theoretical_vs_detected_fdr_seqx",seqdepth[n],"_",weight_method[k]))+geom_abline(slope=1, intercept=0, linetype=2,color="grey")+expand_limits(y=c(0,0.3)))
    }
  }
}
}

#HiCbinpairs_data$sig<-c()
#HiCbinpairs_data$sig[HiCbinpairs_data$Group %in% not_sig_group1$Group.1]<-0

