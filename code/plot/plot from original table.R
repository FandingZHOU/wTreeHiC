dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
meth=c("multiHiCcompare","diffHic")
doc_name<-c("_Multi_glm","_diffHic.tsv")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/originaltable/simulation","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/originaltable/simulation")
in_path2="/p/keles/fandingzhou/volumeA/DCI/simulation"
options("scipen"=100)
True_sig=scan(paste0(in_path2,"/True_signal5"),what = character())

seqdepth<-c(1,3,5)
scalesize<-c(2,3,5,10)
bin_size=40000

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
#m=1 p.value
m=2
    for (foldchange in c(2,3,4,6)){
      HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/sim3_binpairs_GMvsGMFC",foldchange,doc_name[m]))
      hist(HiCbinpairs_data$PValue,main=paste0(meth[m],"_FC",foldchange))
    }


#Plot original TP
m=2
threshold<-c(0.001,0.005,0.01,0.05,0.1)
power<-c()
FDR<-c()
for (foldchange in c(2,3,4,6)){
    HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/sim3_binpairs_GMvsGMFC",foldchange,doc_name[m]))
    DCI_ratio<-c()
    FD_ratio<-c()
    for(i in 1:length(threshold)){
      index=which(HiCbinpairs_data$FDR<threshold[i])
      detected_sig<-HiCbinpairs_data
      detected_sig<-paste0(detected_sig$end2-bin_size/2,"_",detected_sig$end1-bin_size/2)
      DCI_ratio[i]<-length(intersect(detected_sig[index],True_sig))/length(intersect(detected_sig,True_sig))
      FD_ratio[i]<-1-length(intersect(detected_sig[index],True_sig))/length(index)
      if(is.na(FD_ratio[i])){FD_ratio[i]=0}
    }
    power<-c(power,DCI_ratio)
    FDR<-c(FDR,FD_ratio)
}   
   supp=c(rep("FC2",length(threshold)),rep("FC3",length(threshold)),rep("FC4",length(threshold)),rep("FC6",length(threshold)))
    theoretical_fdr=rep(threshold,4)
    tgg=data.frame(supp,theoretical_fdr,power)
    print(ggplot(tgg, aes(x=theoretical_fdr, y=power, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_theoretical_fdrvs_power")))
    tgg=data.frame(supp,theoretical_fdr,FDR)
    print(ggplot(tgg, aes(x=theoretical_fdr, y=FDR, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_theoretical_vs_detected_fdr"))+geom_abline(slope=1, intercept=0, linetype=2,color="grey")+expand_limits(y=c(0,1)))
    
    m=1
    threshold<-c(0.001,0.005,0.01,0.05,0.1)
    power<-c()
    FDR<-c()
    for (foldchange in c(2,3,4,6)){
      HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/sim3_binpairs_GMvsGMFC",foldchange,doc_name[m]))
      DCI_ratio<-c()
      FD_ratio<-c()
      for(i in 1:length(threshold)){
        index=which(p.adjust(HiCbinpairs_data$p.value,"BH")<threshold[i])
        detected_sig<-HiCbinpairs_data
        detected_sig<-paste0(detected_sig$region1+bin_size/2,"_",detected_sig$region2+bin_size/2)
        DCI_ratio[i]<-length(intersect(detected_sig[index],True_sig))/length(intersect(detected_sig,True_sig))
        FD_ratio[i]<-1-length(intersect(detected_sig[index],True_sig))/length(index)
        if(is.na(FD_ratio[i])){FD_ratio[i]=0}
      }
      power<-c(power,DCI_ratio)
      FDR<-c(FDR,FD_ratio)
    }   
    supp=c(rep("FC2",length(threshold)),rep("FC3",length(threshold)),rep("FC4",length(threshold)),rep("FC6",length(threshold)))
    theoretical_fdr=rep(threshold,4)
    tgg=data.frame(supp,theoretical_fdr,power)
    print(ggplot(tgg, aes(x=theoretical_fdr, y=power, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_theoretical_fdrvs_power")))
    tgg=data.frame(supp,theoretical_fdr,FDR)
    print(ggplot(tgg, aes(x=theoretical_fdr, y=FDR, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_theoretical_vs_detected_fdr"))+geom_abline(slope=1, intercept=0, linetype=2,color="grey")+expand_limits(y=c(0,1)))
    