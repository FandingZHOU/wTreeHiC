options("scipen"=100)

dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path="/p/keles/fandingzhou/volumeA/TAD/GM12878"
in_path2="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/filter"
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/weighttable","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/weighttable/filter")
func_path="/p/keles/fandingzhou/volumeA/DCI/code/tree/"
GM_spike_in_scale10_seqx5_aggregated_contactmatrix.tsv
meth=c("multiHiCcompare","diffHic")
seqdepth<-c(1,3,5)
group<-c("A","GM12","GM13","GM14")
doc_name<-c("_Multi_glm","_diffHic.tsv")
group_size<-c("regular","large","small")
weight_method<-c("xi","pi0")
scalesize<-c(2,3,5,10)

library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(ggplot2,lib.loc = dict_path)


#source(paste0(func_path,"filter_group_function.R"))
dci_or_fdr="FDR"

bin_size<-40000
ch1_length<-249250621

#dci_or_fdr="DCI"
gettype="power"
#for(m in 1:2){
for(gettype in c("FDR","DCI","power")){

for(m in 1:2){
  tgg<-data.frame()
for(ss in 0:4){
  for(n in 1:3){
    for(k in 1:2){
      threshold<-c(0.001,0.005,0.01,0.05,0.1)
      if(ss==0){
        HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
        #TSetKey<-scan(paste0(in_path2,"/_GMvsvsGMspikein_seqx",seqdepth[n],"_smoothTSet"),what = character())
        TSetKey<-read.csv(paste0(in_path2,"/GMvsvsGMspikein_seqx",seqdepth[n],"TSetmatrix"),stringsAsFactors=F)
        TADs<-read_tsv(paste0(in_path,"/OnTAD_GM_spike_in_seqx",seqdepth[n],"_aggregated.tad"),col_names=F)
      }
      if(ss>0){
        TSetKey<-read.csv(paste0(in_path2,"/GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"TSetmatrix"),stringsAsFactors=F)
        #TSetKey<-scan(paste0(in_path2,"/_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_smoothTSet"),what = character())
        HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
        TADs<-read_tsv(paste0(in_path,"/OnTAD_GM_spike_in_scale",scalesize[ss],"_seqx",seqdepth[n],"_aggregated.tad"),col_names=F)
      }
      HiCbinpairs_data<-tree_layer(HiCbinpairs_data,TADs,bin_size)
      max_layer<-max(TADs$X3)-2
      if(k==1){
        
        if(gettype=="power"){
          power_weight<-get_fdr_dci(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="power")
          power_orig<-get_fdr_dci(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="power")
          power_tree<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="power",max_layer)
          power_wtree<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="power",max_layer)
        }
        if(gettype=="FDR"){
          FDR_weight_detected<-get_fdr_dci(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="FDR")
          FDR_orig_detected<-get_fdr_dci(HiCbinpairs_data,add_weight = FALSE,threshold,TSetKey,gettype="FDR")
          FDR_tree_detected<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="FDR",max_layer)
          FDR_wtree_detected<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="FDR",max_layer)
        }
        if(gettype=="DCI"){
          DCI_weight_true<-get_fdr_dci(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="DCI")
          DCI_orig_true<-get_fdr_dci(HiCbinpairs_data,add_weight =FALSE,threshold,TSetKey,gettype="DCI")
          DCI_tree_true<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="DCI",max_layer)
          DCI_wtree_true<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="DCI",max_layer)
          
        }
        
        }
      if(k==2){
        if(gettype=="power"){
          power_weight2<-get_fdr_dci(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="power")
          power_wtree2<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="power",max_layer)
        }
        if(gettype=="FDR"){
          FDR_weight_detected2<-get_fdr_dci(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="FDR")
          FDR_wtree_detected2<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="FDR",max_layer)
        }
        if(gettype=="DCI"){
          DCI_weight_true2<-get_fdr_dci(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="DCI")
          DCI_wtree_true2<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="DCI",max_layer)
        }
      }
      
      supp=c(rep("orignal",5),rep("weight",5),rep("tree",5),rep("wtree",5),rep("weight2",5),rep("wtree2",5))
      theoretical_fdr=rep(threshold,6)
      if(ss==0){
        scalesizename<-rep("scale1",length(theoretical_fdr))
      }
      if(ss>0){
        scalesizename<-rep(paste0("scale",scalesize[ss]),length(theoretical_fdr))
      }
      seqdepthname<-rep(paste0("seqx",seqdepth[n]),length(theoretical_fdr))
    }
    if (gettype=="DCI"){
      true_positive=c(DCI_orig_true,DCI_weight_true,DCI_tree_true,DCI_wtree_true,DCI_weight_true2,DCI_wtree_true2)
      tgg=rbind(tgg,data.frame(supp,theoretical_fdr,true_positive,scalesizename,seqdepthname))
    }
    if(gettype=="power"){
      true_power=c(power_orig,power_weight,power_tree,power_wtree,power_weight2,power_wtree2)
      tgg=rbind(tgg,data.frame(supp,theoretical_fdr,true_power,scalesizename,seqdepthname))
    }
    if(gettype=="FDR"){
      detected_fdr=c(FDR_orig_detected,FDR_weight_detected,FDR_tree_detected,FDR_wtree_detected,FDR_weight_detected2,FDR_wtree_detected2)
      tgg=rbind(tgg,data.frame(supp,theoretical_fdr,detected_fdr,scalesizename,seqdepthname))
    }
  }
}
if(gettype=="power"){
write.table(tgg,paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_power"))
print(ggplot(tgg, aes(x=theoretical_fdr, y=true_power, color=supp,shape=supp)) + geom_line(size=0.5) +geom_point(size=1)+ facet_grid(scalesizename~seqdepthname)+labs(title = paste0(meth[m],"_theoretical_fdr_vs_tree_power")))
}
  
if (gettype=="DCI"){
  write.table(tgg,paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_DCI"))
  print(ggplot(tgg, aes(x=theoretical_fdr, y=true_positive, color=supp,shape=supp)) 
        + geom_line(size=0.35) +geom_point(size=1)+ facet_grid(scalesizename~seqdepthname)
        +labs(title = paste0(meth[m],"_theoretical_fdr_vs_tree_tp")))
}
if(gettype=="FDR"){
  write.table(tgg,paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_FDR"))
  print(ggplot(tgg, aes(x=theoretical_fdr, y=detected_fdr, color=supp,shape=supp)) 
        + geom_line(size=0.35) +geom_point(size=1)+ geom_abline(slope=1, intercept=0, linetype=2,color="grey")
        +expand_limits(y=c(0,0.3))+ facet_grid(scalesizename~seqdepthname)+ labs(title = paste0(meth[m],"_theoretical_vs_detected_fdr")))

}
}
}
tgg<-read.table(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_power"))


tgg<-data.frame(tgg[451:900,])
for(ss in 1:4){
for(dci_or_fdr in c("DCI","FDR")){
for(k in 1:2){
  for(m in 1:2){
    for(n in 1:3){
      TSetKey<-scan(paste0(in_path2,"/_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_smoothTSet"),what = character())
      HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"))
      threshold<-c(0.001,0.005,0.01,0.05,0.1)
      supp=c(rep("orignal",5),rep("weight",5),rep("tree",5),rep("wtree",5))
      theoretical_fdr=rep(threshold,4)
      TADs<-read_tsv(paste0(in_path,"/OnTAD_GM_spike_in_seqx",seqdepth[n],"_aggregated.tad"),col_names=F)
      
      if (dci_or_fdr=="DCI"){
      DCI_weight_true<-get_fdr_dci(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="DCI")
      DCI_orig_true<-get_fdr_dci(HiCbinpairs_data,add_weight =FALSE,threshold,TSetKey,gettype="DCI")
      
      #get tree structure
      HiCbinpairs_data<-tree_layer(HiCbinpairs_data,TADs,bin_size)
      max_layer<-max(TADs$X3)
      DCI_tree_true<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="DCI",max_layer)
      DCI_wtree_true<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="DCI",max_layer)
      
      true_positive=c(DCI_orig_true,DCI_weight_true,DCI_tree_true,DCI_wtree_true)
      tgg=data.frame(supp,theoretical_fdr,true_positive)
      print(ggplot(tgg, aes(x=theoretical_fdr, y=true_positive, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_theoretical_fdr_vs_tree_tp_seqx",seqdepth[n],"_",weight_method[k],"_scale",scalesize[ss])))
      }else{
        FDR_weight_detected<-get_fdr_dci(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="FDR")
        FDR_orig_detected<-get_fdr_dci(HiCbinpairs_data,add_weight = FALSE,threshold,TSetKey,gettype="FDR")
        
        HiCbinpairs_data<-tree_layer(HiCbinpairs_data,TADs,bin_size)
        max_layer<-max(TADs$X3)
       
        FDR_tree_detected<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="FDR",max_layer)
        FDR_wtree_detected<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,gettype="FDR",max_layer)
        
        detected_fdr=c(FDR_orig_detected,FDR_weight_detected,FDR_tree_detected,FDR_wtree_detected)
        tgg=data.frame(supp,theoretical_fdr,detected_fdr)
        print(ggplot(tgg, aes(x=theoretical_fdr, y=detected_fdr, color=supp,shape=supp)) + geom_line(size=1) +geom_point(size=2)+ labs(title = paste0(meth[m],"_theoretical_vs_detected_fdr_seqx",seqdepth[n],"_",weight_method[k],"_scale",scalesize[ss]))+geom_abline(slope=1, intercept=0, linetype=2,color="grey")+expand_limits(y=c(0,0.3)))
      }
    }
  }
}}}

