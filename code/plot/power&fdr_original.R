##plot
dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
meth=c("multiHiCcompare","diffHic")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/weighttable","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/weighttable")
group<-c("A","GM12","GM13","GM14")
library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(ggplot2)

weight_method<-c("xi","pi0")
#Power
n=1
for (m in 1:2){
  for(k in 1:2){
    #HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_weightedp_",weight_method[k],"_small.tsv"))
    #HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_weightedp_",weight_method[k],"_large.tsv"))
    HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_weightedp_",weight_method[k],".tsv"))
    DCI_weight<-c()
    threshold<-c(0.001,0.005,0.01,0.05,0.1)
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data$weighted_P),])
    for (j in 1:5){
      maxpi<-max(HiCbinpairs_data[HiCbinpairs_data$adj_w_p<=threshold[j],]$weighted_P)
      DCI_weight[j]=nrow(HiCbinpairs_data[HiCbinpairs_data$weighted_P<=maxpi,])
    }
    adj_p<-c()
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data$P),])
    for (i in 1:nrow(HiCbinpairs_data)){
      adj_p[i]<-nrow(HiCbinpairs_data)*HiCbinpairs_data$P[i]/i
    }
    DCI_orig<-c()
    for (j in 1:5){
      maxpi<-max(HiCbinpairs_data[adj_p<=threshold[j],]$P)
      DCI_orig[j]=nrow(HiCbinpairs_data[HiCbinpairs_data$P<=maxpi,])
    }
    DCIOUT<-cbind(DCI_weight,DCI_orig)
    write.table(DCIOUT, file=paste0(out_path[m],"/DCI_GM_A_",weight_method[k],".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
    supp=c(rep("orignal",5),rep("weighted",5))
    threshold=rep(threshold,2)
    DCI=c(DCI_orig,DCI_weight)
    tgg=data.frame(supp,threshold,DCI)
    
   #ggplot(tgg, aes(x=threshold, y=DCI, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_power_",weight_method[k],"_small"))
    #ggplot(tgg, aes(x=threshold, y=DCI, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_power_",weight_method[k],"_large"))
    print(ggplot(tgg, aes(x=threshold, y=DCI, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_power_",weight_method[k])))
  }
}

n=1
for (m in 1:2){
    HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_weightedp_",weight_method[k],".tsv"))
    hist(HiCbinpairs_data$P)
}

#FDR
n=3
for (m in 1:2){
  for(k in 1:2){
    HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_weightedp_",weight_method[k],"_small.tsv"))
    #HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_weightedp_",weight_method[k],".tsv"))
    DCI_weight<-c()
    threshold<-c(0.001,0.005,0.01,0.05,0.1)
    for (j in 1:5){
      maxpi<-max(HiCbinpairs_data[HiCbinpairs_data$adj_w_p<=threshold[j],]$weighted_P)
      DCI_weight[j]=nrow(HiCbinpairs_data[HiCbinpairs_data$weighted_P<=maxpi,])
    }
    adj_p<-c()
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data$P),])
    for (i in 1:nrow(HiCbinpairs_data)){
      adj_p[i]<-nrow(HiCbinpairs_data)*HiCbinpairs_data$P[i]/i
    }
    DCI_orig<-c()
    for (j in 1:5){
      maxpi<-max(HiCbinpairs_data[adj_p<=threshold[j],]$P)
      DCI_orig[j]=nrow(HiCbinpairs_data[HiCbinpairs_data$P<=maxpi,])
    }
    FDR_weight<-DCI_weight/nrow(HiCbinpairs_data)
    FDR_orig<-DCI_orig/nrow(HiCbinpairs_data)
    supp=c(rep("orignal",5),rep("weighted",5))
    threshold=rep(threshold,2)
    FDR=c(FDR_orig,FDR_weight)
    tgg=data.frame(supp,threshold,FDR)
    library(ggplot2)
    print(ggplot(tgg, aes(x=threshold, y=FDR, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_power_",weight_method[k],"_small"))+geom_abline(slope=1, intercept=0, linetype=2,color="grey")+expand_limits(y=c(0,0.1)))
    #ggplot(tgg, aes(x=threshold, y=FDR, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=1.5)+geom_abline(slope=1, intercept=0, linetype=2,color="grey")+expand_limits(y=c(0,0.1))+labs(title = paste0(meth[m],"_fdr_",weight_method[k]))
    FDROUT<-cbind(FDR_weight,FDR_orig)
   
    #write.table(FDROUT, file=paste0(out_path[m],"/FDR_GM13_",weight_method[k],".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
  }
}






#FDR


#GM_GM13
# DCI_orig<-c(155,743,1456,6365,12494)
# DCI_weight<-c(219,1147,2185,9683,19246)
FDR_weight<- c(0.0000244813,0.0003427382,0.0007344391,0.0103984332,0.0304486199)
FDR_orig<-c(0.0000244813,0.0002631740,0.0005691903,0.0060652427,0.0172164759)
FDR_weight<-c(0.0000244813,0.000281535,0.0005753106,0.0062488524,0.0182140890)
FDR_weight<-c(3.060163e-05,4.835057e-04,1.077177e-03,1.085746e-02,2.528307e-02)
DCI_orig<-c(13844,22474,28075,49369,64329)
#DCI_weight<-c(14962,23374,28576,47008,59410)
DCI_weight<-c(14245,22792,28332,49165,63671)
#DCI_weight<-c(15024,24022,29911,52376,68303)
DCI_weight<-c(14462,23043,28565,49317,63577)
DCI_weight<-c(14893,23599, 29177, 50022, 64037)
DCI_weight2<-c(14893,23599, 29177, 50022, 64037)
DCI_weight2<-c(14462,23043,28565,49317,63577)
supp=c(rep("orignal",5),rep("weighted",5),rep("weighted_small",5))
threshold=rep(threshold,3)
DCI=c(DCI_orig,DCI_weight,DCI_weight2)
tgg=data.frame(supp,threshold,DCI)
library(ggplot2)
ggplot(tgg, aes(x=threshold, y=DCI, color=supp,shape=supp)) + geom_line(size=0.5) +geom_point(size=1.5)