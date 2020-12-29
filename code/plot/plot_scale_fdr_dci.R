dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
meth=c("multiHiCcompare","diffHic")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/weighttable","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/weighttable/filter")
in_path2="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/filter"
options("scipen"=100)
TSetKey<-scan(paste0(in_path2,"/TSetKey"),what = character())

seqdepth<-c(1,3,5)
scalesize<-c(2,3,5,10)

library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(ggplot2,lib.loc = dict_path)
#accuracy
#detected/real
#load("/p/keles/treehic/volumeD/freeHiC/GM12878_A549_TSetKey_trueSpikeinSmooth.RData")
weight_method<-c("xi","pi0")
#dci_or_fdr<-"FDR"
k=1
#plot pvalue
for (m in 1:2){
for(ss in 1:4){
  for (n in 1:3){
    HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
    hist(HiCbinpairs_data$P,main=paste0(meth[m],"scale",scalesize[ss],"seqx",seqdepth[n]))
  }
}
}

#Plot original TP
k=1
for(m in 1:2){
  for(n in 1:3){
    threshold<-c(0.001,0.005,0.01,0.05,0.1)
    power<-c()
    HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
    DCI_orig_true<-get_fdr_dci(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="DCI")
    power<-c(power,DCI_orig_true)
    for (ss in 1:4){
      HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
      DCI_orig_true<-get_fdr_dci(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="DCI")
      power<-c(power,DCI_orig_true)
    }
    supp=c(rep("scale1",length(threshold)),rep("scale2",length(threshold)),rep("scale3",length(threshold)),rep("scale5",length(threshold)),rep("scale10",length(threshold)))
    theoretical_fdr=rep(threshold,5)
    tgg=data.frame(supp,theoretical_fdr,power)
    print(ggplot(tgg, aes(x=theoretical_fdr, y=power, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_theoretical_fdrvs_power_seqx",seqdepth[n])))
  }
}

#plot original FDR
k=1
for(m in 1:2){
  for(n in 1:3){
    threshold<-c(0.001,0.005,0.01,0.05,0.1)
    detected_fdr<-c()
    HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
    FDR_orig_detected<-get_fdr_dci(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="FDR")
    detected_fdr<-c(detected_fdr,FDR_orig_detected)
    for (ss in 1:4){
      HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
      FDR_orig_detected<-get_fdr_dci(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="FDR")
      detected_fdr<-c(detected_fdr,FDR_orig_detected)
    }
    supp=c(rep("scale1",length(threshold)),rep("scale2",length(threshold)),rep("scale3",length(threshold)),rep("scale5",length(threshold)),rep("scale10",length(threshold)))
    theoretical_fdr=rep(threshold,5)
    tgg=data.frame(supp,theoretical_fdr,detected_fdr)
    print(ggplot(tgg, aes(x=theoretical_fdr, y=detected_fdr, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_theoretical_vs_detected_fdr_seqx",seqdepth[n]))+geom_abline(slope=1, intercept=0, linetype=2,color="grey")+expand_limits(y=c(0,1)))
  }
}

#plot weighted FDR
for(m in 1:2){
  for(k in 1:2){
    for(n in 1:3){
      threshold<-c(0.001,0.005,0.01,0.05,0.1)
      detected_fdr<-c()
      HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
      FDR_weight_detected<-get_fdr_dci(HiCbinpairs_data,add_weight=TRUE,threshold,TSetKey,gettype="FDR")
      detected_fdr<-c(detected_fdr,FDR_weight_detected)
      for (ss in 1:4){
        HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
        FDR_weight_detected<-get_fdr_dci(HiCbinpairs_data,add_weight=TRUE,threshold,TSetKey,gettype="FDR")
        detected_fdr<-c(detected_fdr,FDR_weight_detected)
      }
      supp=c(rep("scale1",length(threshold)),rep("scale2",length(threshold)),rep("scale3",length(threshold)),rep("scale5",length(threshold)),rep("scale10",length(threshold)))
      theoretical_fdr=rep(threshold,5)
      tgg=data.frame(supp,theoretical_fdr,detected_fdr)
      print(ggplot(tgg, aes(x=theoretical_fdr, y=detected_fdr, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_theoretical_vs_detected_fdr_seqx",seqdepth[n],"_weight_",weight_method[k]))+geom_abline(slope=1, intercept=0, linetype=2,color="grey")+expand_limits(y=c(0,1)))
    }
  }
}