options("scipen"=100)

dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path="/p/keles/fandingzhou/volumeA/TAD/GM12878"
in_path2="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/filter"
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/weighttable",
           "/p/keles/fandingzhou/volumeA/DCI/diffHic/data/weighttable/filter",
           "/p/keles/fandingzhou/volumeA/DCI/HiCDCPlus/data/weighttable/filter")
func_path="/p/keles/fandingzhou/volumeA/DCI/code/tree/"
func_path2="/p/keles/fandingzhou/volumeA/DCI/code/tree/plot"
meth=c("multiHiCcompare","diffHic","HiCDCPlus")
seqdepth<-c(1,3,5)
group<-c("A","GM12","GM13","GM14")
doc_name<-c("_Multi_glm","_diffHic.tsv")
group_size<-c("regular","large","small")
weight_method<-c("xi","pi0","IHW")
scalesize<-c(2,3,5,10)

library(cli,lib.loc = dict_path)
library(fansi,lib.loc=dict_path)
library(utf8,lib.loc=dict_path)
library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(ggplot2,lib.loc = dict_path)


#source(paste0(func_path,"filter_group_function.R"))
#source(paste0(func_path,"get_fdr_tp_tree_function.R"))
#source(paste0(func_path,"tree_layer_function.R"))
#source(paste0(func_path2,"get_fdr_tp_function.R"))

bin_size<-40000
ch1_length<-249250621

#dci_or_fdr="DCI"
#new sim
out_path[1]=paste0(out_path[1],"/simulation/")
out_path[2]="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/weighttable/simulation/filter/"
out_path[3]=paste0(out_path[3],"/simulation/")
in_path3="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/filter/simulation"
out_path[4]="/p/keles/fandingzhou/volumeA/DCI/plots/data/"
for(foldchange in c(3,4,6)){
   for(m in 1:3){
     k=1
      HiCbinpairs_data<-read_tsv(paste0(out_path[m],meth[m],"_sim3_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],".tsv"))
      hist(HiCbinpairs_data$P,main=paste0(meth[m],"_FC",foldchange))
   }
}

for(seed in 16:20){
  for(test in c("2vs2","4vs4")){
    #if(seed!=1|test!="4vs4"){
    for(n in 1:2){
      #if(seed!=1|n!=1){
        tgg_f<-data.frame()
        tgg_d<-data.frame()
        tgg_p<-data.frame()
        in_path5=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_",meth[n])
        in_path4=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/",test,"/filter_as_",meth[n])
        for(foldchange in c(3,4,6)){
          TADs<-read_tsv(paste0(in_path,"/OnTAD_seed",seed,"_GM_FC",foldchange,"_aggregated.tad"),col_names=F)
          Truesig=scan(paste0(in_path5,"/Trueset/GMvsGMFC",foldchange,"_TSet"),what = character())
          TSetKey=data.frame(Truesig,Truesig)
          colnames(TSetKey)=c("smooth","origin")
          for(m in 1:3){
            threshold<-c(0.001,0.005,0.01,0.05,0.1)
            len_thre<-length(threshold)
            supp=c(rep("orignal",len_thre),rep("weight",len_thre),rep("tree",len_thre),rep("wtree",len_thre),rep("weight2",len_thre),rep("wtree2",len_thre))
            theoretical_fdr=rep(threshold,6)
            FCname<-rep(paste0("FC",foldchange),length(theoretical_fdr))
            methname<-rep(meth[m],length(theoretical_fdr))
            
            for(k in 1:2){
              #threshold<-c(0.001,0.005,0.01,0.05,0.1)
              HiCbinpairs_data<-read_tsv(paste0(in_path4,"/",meth[m],"/weighttable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],".tsv"))
              #HiCbinpairs_data<-read_tsv(paste0(out_path[m],meth[m],"_sim3_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],".tsv"))
              HiCbinpairs_data<-na.omit(HiCbinpairs_data)
              
              HiCbinpairs_data<-tree_layer(HiCbinpairs_data,TADs,bin_size)
              max_layer<-max(TADs$X3)-3
              if(k==1){
                weight_result<-get_fdr_dci(HiCbinpairs_data,add_weight=T,threshold,TSetKey);
                orig_result<-get_fdr_dci(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey);
                tree_result<-get_fdr_dci_tree(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,max_layer);
                wtree_result<-get_fdr_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,max_layer,weight_method[k]);
                
              }
              if(k==2){
                weight2_result<-get_fdr_dci(HiCbinpairs_data,add_weight=T,threshold,TSetKey)
                wtree2_result<-get_fdr_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,max_layer,weight_method[k])
              }
            }
            for(gettype in c("FDR","DCI","power")){
              if(gettype=="power"){
                power_weight<-weight_result$power
                power_orig<-orig_result$power
                power_tree<-tree_result$power
                power_wtree<-wtree_result$power
                power_weight2<-weight2_result$power
                power_wtree2<-wtree2_result$power
                true_power=c(power_orig,power_weight,power_tree,power_wtree,power_weight2,power_wtree2)
                tgg_p=rbind(tgg_p,data.frame(supp,theoretical_fdr,true_power,FCname,methname))
                
              } 
              if(gettype=="FDR"){
                FDR_weight_detected<-weight_result$FDR
                FDR_orig_detected<-orig_result$FDR
                FDR_tree_detected<-tree_result$FDR
                FDR_wtree_detected<-wtree_result$FDR
                FDR_weight_detected2<-weight2_result$FDR
                FDR_wtree_detected2<-wtree2_result$FDR
                detected_fdr=c(FDR_orig_detected,FDR_weight_detected,FDR_tree_detected,FDR_wtree_detected,FDR_weight_detected2,FDR_wtree_detected2)
                tgg_f=rbind(tgg_f,data.frame(supp,theoretical_fdr,detected_fdr,FCname,methname))
                
              }
              if(gettype=="DCI"){
                DCI_weight_true<-weight_result$DCI
                DCI_orig_true<-orig_result$DCI
                DCI_tree_true<-tree_result$DCI
                DCI_wtree_true<-wtree_result$DCI
                DCI_weight_true2<-wtree_result$DCI
                DCI_wtree_true2<-wtree2_result$DCI
                true_positive=c(DCI_orig_true,DCI_weight_true,DCI_tree_true,DCI_wtree_true,DCI_weight_true2,DCI_wtree_true2)
                tgg_d=rbind(tgg_d,data.frame(supp,theoretical_fdr,true_positive,FCname,methname))
                
              }
              
            }
            
          }
        }
        write.table(tgg_p,paste0(in_path4,"/GMvsGM_power"))
        write.table(tgg_d,paste0(in_path4,"/GMvsGM_DCI"))
        write.table(tgg_f,paste0(in_path4,"/GMvsGM_FDR"))
      }
    }
    
  }
 # }
}

test="4vs4"
n=1
for(seed in 3:10){
        tgg_f<-data.frame()
        tgg_d<-data.frame()
        tgg_p<-data.frame()
        in_path5=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_",meth[n])
        in_path4=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/",test,"/filter_as_",meth[n])
          for(foldchange in c(3,4,6)){
            TADs<-read_tsv(paste0(in_path,"/OnTAD_seed",seed,"_GM_FC",foldchange,"_aggregated.tad"),col_names=F)
            Truesig=scan(paste0(in_path5,"/Trueset/GMvsGMFC",foldchange,"_TSet"),what = character())
            TSetKey=data.frame(Truesig,Truesig)
            colnames(TSetKey)=c("smooth","origin")
            for(m in 1:3){
              threshold<-c(0.001,0.005,0.01,0.05,0.1)
              len_thre<-length(threshold)
              supp=c(rep("wtree3",len_thre),rep("wtree4",len_thre))
              theoretical_fdr=rep(threshold,2)
              FCname<-rep(paste0("FC",foldchange),length(theoretical_fdr))
              methname<-rep(meth[m],length(theoretical_fdr))
              k=1
              HiCbinpairs_data<-read_tsv(paste0(in_path4,"/",meth[m],"/weighttable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],".tsv"))
              #HiCbinpairs_data<-read_tsv(paste0(out_path[m],meth[m],"_sim3_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],".tsv"))
              HiCbinpairs_data<-na.omit(HiCbinpairs_data)
              HiCbinpairs_data<-tree_layer(HiCbinpairs_data,TADs,bin_size)
              max_layer<-max(TADs$X3)-5
              wtree_result<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,max_layer,weight_method[1]);
              wtree2_result<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,max_layer,weight_method[2])
              for(gettype in c("FDR","DCI","power")){
                if(gettype=="power"){
                  power_wtree<-wtree_result$power
                  power_wtree2<-wtree2_result$power
                  true_power=c(power_wtree,power_wtree2)
                  tgg_p=rbind(tgg_p,data.frame(supp,theoretical_fdr,true_power,FCname,methname))
                  
                } 
                if(gettype=="FDR"){
                  FDR_wtree_detected<-wtree_result$FDR
                  FDR_wtree_detected2<-wtree2_result$FDR
                  detected_fdr=c(FDR_wtree_detected,FDR_wtree_detected2)
                  tgg_f=rbind(tgg_f,data.frame(supp,theoretical_fdr,detected_fdr,FCname,methname))
                  
                }
                if(gettype=="DCI"){
                  DCI_wtree_true<-wtree_result$DCI
                  DCI_wtree_true2<-wtree2_result$DCI
                  true_positive=c(DCI_wtree_true,DCI_wtree_true2)
                  tgg_d=rbind(tgg_d,data.frame(supp,theoretical_fdr,true_positive,FCname,methname))
                  
                }
                
              }
              
            }
          }
          write.table(tgg_p,paste0(in_path4,"/GMvsGM_power_new_wtree"))
          write.table(tgg_d,paste0(in_path4,"/GMvsGM_DCI_new_wtree"))
          write.table(tgg_f,paste0(in_path4,"/GMvsGM_FDR_new_wtree"))
}

for(seed in 11){
  for(test in c("2vs2","4vs4")){
    #if(seed!=1|test!="4vs4"){
    for(n in 1:2){
      #if(seed!=1|n!=1){
      tgg_f<-data.frame()
      tgg_d<-data.frame()
      tgg_p<-data.frame()
      in_path5=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_",meth[n])
      in_path4=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/",test,"/filter_as_",meth[n])
      for(foldchange in c(3,4,6)){
        TADs<-read_tsv(paste0(in_path,"/OnTAD_seed",seed,"_GM_FC",foldchange,"_aggregated.tad"),col_names=F)
        Truesig=scan(paste0(in_path5,"/Trueset/GMvsGMFC",foldchange,"_TSet"),what = character())
        TSetKey=data.frame(Truesig,Truesig)
        colnames(TSetKey)=c("smooth","origin")
        for(m in 1:3){
          threshold<-c(0.001,0.005,0.01,0.05,0.1)
          len_thre<-length(threshold)
          supp=c(rep("weightd",len_thre),rep("wtreed",len_thre),rep("weight2d",len_thre),rep("wtree2d",len_thre))
          theoretical_fdr=rep(threshold,4)
          FCname<-rep(paste0("FC",foldchange),length(theoretical_fdr))
          methname<-rep(meth[m],length(theoretical_fdr))
          
          for(k in 1:2){
            #threshold<-c(0.001,0.005,0.01,0.05,0.1)
            HiCbinpairs_data<-read_tsv(paste0(in_path4,"/",meth[m],"/weighttable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],"_distance.tsv"))
            #HiCbinpairs_data<-read_tsv(paste0(out_path[m],meth[m],"_sim3_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],".tsv"))
            HiCbinpairs_data<-na.omit(HiCbinpairs_data)
            
            HiCbinpairs_data<-tree_layer(HiCbinpairs_data,TADs,bin_size)
            max_layer<-max(TADs$X3)-3
            if(k==1){
              weight_result<-get_fdr_dci(HiCbinpairs_data,add_weight=T,threshold,TSetKey);
              wtree_result<-get_fdr_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,max_layer,weight_method[k]);
              
            }
            if(k==2){
              weight2_result<-get_fdr_dci(HiCbinpairs_data,add_weight=T,threshold,TSetKey)
              wtree2_result<-get_fdr_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,max_layer,weight_method[k])
            }
          }
          for(gettype in c("FDR","DCI","power")){
            if(gettype=="power"){
              power_weight<-weight_result$power
              power_wtree<-wtree_result$power
              power_weight2<-weight2_result$power
              power_wtree2<-wtree2_result$power
              true_power=c(power_weight,power_wtree,power_weight2,power_wtree2)
              tgg_p=rbind(tgg_p,data.frame(supp,theoretical_fdr,true_power,FCname,methname))
              
            } 
            if(gettype=="FDR"){
              FDR_weight_detected<-weight_result$FDR
              FDR_wtree_detected<-wtree_result$FDR
              FDR_weight_detected2<-weight2_result$FDR
              FDR_wtree_detected2<-wtree2_result$FDR
              detected_fdr=c(FDR_weight_detected,FDR_wtree_detected,FDR_weight_detected2,FDR_wtree_detected2)
              tgg_f=rbind(tgg_f,data.frame(supp,theoretical_fdr,detected_fdr,FCname,methname))
              
            }
            if(gettype=="DCI"){
              DCI_weight_true<-weight_result$DCI
              DCI_wtree_true<-wtree_result$DCI
              DCI_weight_true2<-wtree_result$DCI
              DCI_wtree_true2<-wtree2_result$DCI
              true_positive=c(DCI_weight_true,DCI_wtree_true,DCI_weight_true2,DCI_wtree_true2)
              tgg_d=rbind(tgg_d,data.frame(supp,theoretical_fdr,true_positive,FCname,methname))
              
            }
            
          }
          
        }
      }
      write.table(tgg_p,paste0(in_path4,"/GMvsGM_power_d"))
      write.table(tgg_d,paste0(in_path4,"/GMvsGM_DCI_d"))
      write.table(tgg_f,paste0(in_path4,"/GMvsGM_FDR_d"))
    }
  }
  
}
k=3
for(seed in 11){
  for(test in c("2vs2","4vs4")){
    #if(seed!=1|test!="4vs4"){
    for(n in 1:2){
      #if(seed!=1|n!=1){
      tgg_f<-data.frame()
      tgg_d<-data.frame()
      tgg_p<-data.frame()
      in_path5=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_",meth[n])
      in_path4=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/",test,"/filter_as_",meth[n])
      for(foldchange in c(3,4,6)){
        TADs<-read_tsv(paste0(in_path,"/OnTAD_seed",seed,"_GM_FC",foldchange,"_aggregated.tad"),col_names=F)
        Truesig=scan(paste0(in_path5,"/Trueset/GMvsGMFC",foldchange,"_TSet"),what = character())
        TSetKey=data.frame(Truesig,Truesig)
        colnames(TSetKey)=c("smooth","origin")
        for(m in 1:3){
          threshold<-c(0.001,0.005,0.01,0.05,0.1)
          len_thre<-length(threshold)
          supp=c(rep("weight3d",len_thre),rep("wtree3d",len_thre))
          theoretical_fdr=rep(threshold,2)
          FCname<-rep(paste0("FC",foldchange),length(theoretical_fdr))
          methname<-rep(meth[m],length(theoretical_fdr))
          
          #for(k in 1:2){
            #threshold<-c(0.001,0.005,0.01,0.05,0.1)
            HiCbinpairs_data<-read_tsv(paste0(in_path4,"/",meth[m],"/weighttable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],"_distance.tsv"))
            #HiCbinpairs_data<-read_tsv(paste0(out_path[m],meth[m],"_sim3_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],".tsv"))
            HiCbinpairs_data<-na.omit(HiCbinpairs_data)
            
            HiCbinpairs_data<-tree_layer(HiCbinpairs_data,TADs,bin_size)
            max_layer<-max(TADs$X3)-3
            weight_result<-get_fdr_dci(HiCbinpairs_data,add_weight=T,threshold,TSetKey);
            wtree_result<-get_fdr_dci_tree(HiCbinpairs_data,add_weight=T,threshold,TSetKey,max_layer,weight_method[k]);
            
          }
          for(gettype in c("FDR","DCI","power")){
            if(gettype=="power"){
              power_weight<-weight_result$power
              power_wtree<-wtree_result$power
              true_power=c(power_weight,power_wtree)
              tgg_p=rbind(tgg_p,data.frame(supp,theoretical_fdr,true_power,FCname,methname))
              
            } 
            if(gettype=="FDR"){
              FDR_weight_detected<-weight_result$FDR
              FDR_wtree_detected<-wtree_result$FDR
              detected_fdr=c(FDR_weight_detected,FDR_wtree_detected)
              tgg_f=rbind(tgg_f,data.frame(supp,theoretical_fdr,detected_fdr,FCname,methname))
              
            }
            if(gettype=="DCI"){
              DCI_weight_true<-weight_result$DCI
              DCI_wtree_true<-wtree_result$DCI
              true_positive=c(DCI_weight_true,DCI_wtree_true)
              tgg_d=rbind(tgg_d,data.frame(supp,theoretical_fdr,true_positive,FCname,methname))
              
            }
            
          }
          
        }
      }
      write.table(tgg_p,paste0(in_path4,"/GMvsGM_power_d3"))
      write.table(tgg_d,paste0(in_path4,"/GMvsGM_DCI_d3"))
      write.table(tgg_f,paste0(in_path4,"/GMvsGM_FDR_d3"))
    }
  }
  
}


tgg_p<-data.frame()
tgg_f<-data.frame()
seed=11
test="4vs4"
for(n in 1:2){
  
}
n=1
in_path4=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/",test,"/filter_as_",meth[n])
tgg_p<-rbind(tgg_p,read.table(paste0(in_path4,"/GMvsGM_power")))
tgg_p<-rbind(tgg_p,read.table(paste0(in_path4,"/GMvsGM_power_d")))
tgg_p<-rbind(tgg_p,read.table(paste0(in_path4,"/GMvsGM_power_d3")))
tgg_f<-rbind(tgg_f,read.table(paste0(in_path4,"/GMvsGM_FDR")))
tgg_f<-rbind(tgg_f,read.table(paste0(in_path4,"/GMvsGM_FDR_d")))
tgg_f<-rbind(tgg_f,read.table(paste0(in_path4,"/GMvsGM_FDR_d3")))
tgg_f1<-tgg_f[tgg_f$supp %in% c("orignal", "tree","weight","wtree", "weightd","wtreed"),]
tgg_f1<-tgg_f[tgg_f$supp %in% c("orignal", "tree","weight2","wtree2", "weight2d","wtree2d"),]
tgg_f3<-tgg_f[tgg_f$supp %in% c("orignal", "weightd","weight2d","weight3d"),]
tgg_f4<-tgg_f[tgg_f$supp %in% c("orignal", "wtreed","wtree2d","wtree3d"),]

print(ggplot(tgg_p, aes(x=theoretical_fdr, y=true_power, color=supp,shape=supp)) + geom_line(size=0.5) +geom_point(size=1)+ facet_grid(FCname~methname)+labs(title = paste0(meth[n],"_theoretical_fdr_vs_tree_power")))
print(ggplot(tgg_f, aes(x=theoretical_fdr, y=detected_fdr, color=supp,shape=supp))
      + geom_line(size=0.35) +geom_point(size=1)+ geom_abline(slope=1, intercept=0, linetype=2,color="grey")
      +expand_limits(y=c(0,0.3))+ facet_grid(FCname~methname)+ labs(title = paste0(meth[n],"_theoretical_vs_detected_fdr")))


if (gettype=="DCI"){
  write.table(tgg,paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_DCI"))
  print(ggplot(tgg, aes(x=theoretical_fdr, y=true_positive, color=supp,shape=supp))
        + geom_line(size=0.35) +geom_point(size=1)+ facet_grid(FCname~methname)
        +labs(title = paste0(meth[m],"_theoretical_fdr_vs_tree_tp")))
}
if(gettype=="FDR"){
  write.table(tgg,paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_FDR"))
  print(ggplot(tgg, aes(x=theoretical_fdr, y=detected_fdr, color=supp,shape=supp))
        + geom_line(size=0.35) +geom_point(size=1)+ geom_abline(slope=1, intercept=0, linetype=2,color="grey")
        +expand_limits(y=c(0,0.3))+ facet_grid(FCname~methname)+ labs(title = paste0(meth[m],"_theoretical_vs_detected_fdr")))
  


library(ggplot2)
meth=c("multiHiCcompare","diffHic","HiCDCPlus")

#for(test in c("2vs2","4vs4")){
 # for(n in 1:2){
test="4vs4"
n=1
    tgg_p<-data.frame()
    tgg_f<-data.frame()
    for(seed in 1:10){
      in_path4=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/",test,"/filter_as_",meth[n])
      tgg_p<-rbind(tgg_p,read.table(paste0(in_path4,"/GMvsGM_power")))
      tgg_p<-rbind(tgg_p,read.table(paste0(in_path4,"/GMvsGM_power_new_wtree")))
      tgg_f<-rbind(tgg_f,read.table(paste0(in_path4,"/GMvsGM_FDR")))
      tgg_f<-rbind(tgg_f,read.table(paste0(in_path4,"/GMvsGM_FDR_new_wtree")))
    }
    tgg_p<-tgg_p[tgg_p$supp %in% c("tree","wtree",paste0("wtree",c(2,3,4))),]
    tgg_f<-tgg_f[tgg_f$supp %in% c("tree","wtree",paste0("wtree",c(2,3,4))),]
    tgg_p$seed=rep(paste0("seed",seed),each=nrow(tgg_p)/10)
    tgg_f$seed=rep(paste0("seed",seed),each=nrow(tgg_f)/10)
    tgg_p_0.1<-tgg_p[tgg_p$theoretical_fdr>=0.01,]
    tgg_f_0.1<-tgg_f[tgg_p$theoretical_fdr>=0.01,]
    tgg_f_0.1[is.na(tgg_f_0.1)] <- 0
    
    print(ggplot(tgg_p_0.1, aes(x=as.factor(theoretical_fdr), y=true_power, color=supp,fill=supp))
          +  geom_boxplot()+ facet_grid(FCname~methname)
          +labs(title = paste0(test,"_filter_by_",meth[n],"_theoretical_fdr_vs_tree_tp"))
    )   
    print(ggplot(tgg_f_0.1, aes(x=as.factor(theoretical_fdr), y=detected_fdr, color=supp,fill=supp))
          +  geom_boxplot()+ facet_grid(FCname~methname)
          +labs(title = paste0(test,"_filter_by_",meth[n],"_theoretical_fdr_vs_detected_fdr"))
    )   
    
  }
}


library(ggplot2)
meth=c("multiHiCcompare","diffHic","HiCDCPlus")

for(test in c("2vs2","4vs4")){
  for(n in 1:2){
    tgg_p<-data.frame()
    tgg_f<-data.frame()
    for(seed in 1:10){
      in_path4=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/",test,"/filter_as_",meth[n])
      tgg_p<-rbind(tgg_p,read.table(paste0(in_path4,"/GMvsGM_power")))
      tgg_f<-rbind(tgg_f,read.table(paste0(in_path4,"/GMvsGM_FDR")))
    }
    tgg_p$seed=rep(paste0("seed",seed),each=nrow(tgg_p)/10)
    tgg_f$seed=rep(paste0("seed",seed),each=nrow(tgg_f)/10)
    tgg_p_0.1<-tgg_p[tgg_p$theoretical_fdr>=0.01,]
    tgg_f_0.1<-tgg_f[tgg_p$theoretical_fdr>=0.01,]
    tgg_f_0.1[is.na(tgg_f_0.1)] <- 0
    
    print(ggplot(tgg_p_0.1, aes(x=as.factor(theoretical_fdr), y=true_power, color=supp,fill=supp))
          +  geom_boxplot()+ facet_grid(FCname~methname)
          +labs(title = paste0(test,"_filter_by_",meth[n],"_theoretical_fdr_vs_tree_tp"))
    )   
    print(ggplot(tgg_f_0.1, aes(x=as.factor(theoretical_fdr), y=detected_fdr, color=supp,fill=supp))
          +  geom_boxplot()+ facet_grid(FCname~methname)
          +labs(title = paste0(test,"_filter_by_",meth[n],"_theoretical_fdr_vs_detected_fdr"))
    )   
    
  }
}
   
#     if(gettype=="power"){
#       #write.table(tgg,paste0(out_path[4],"sim3_GMvsGM_power"))
#        print(ggplot(tgg, aes(x=theoretical_fdr, y=true_power, color=supp,shape=supp)) + geom_line(size=0.5) +geom_point(size=1)+ facet_grid(FCname~methname)+labs(title = paste0("filter_by_",meth[n],"_theoretical_fdr_vs_tree_power")))
#     }
#     
#     if (gettype=="DCI"){
#       #write.table(tgg,paste0(out_path[4],"sim3_GMvsGM_DCI"))
#       print(ggplot(tgg, aes(x=theoretical_fdr, y=true_positive, color=supp,shape=supp))
#             + geom_line(size=0.35) +geom_point(size=1)+ facet_grid(FCname~methname)
#             +labs(title = paste0("filter_by_",meth[n],"_theoretical_fdr_vs_tree_tp")))
#     }
#     if(gettype=="FDR"){
#       #write.table(tgg,paste0(out_path[4],"sim3_GMvsGM_FDR"))
#        print(ggplot(tgg, aes(x=theoretical_fdr, y=detected_fdr, color=supp,shape=supp))
#             + geom_line(size=0.35) +geom_point(size=1)+ geom_abline(slope=1, intercept=0, linetype=2,color="grey")
#             +expand_limits(y=c(0,0.3))+ facet_grid(FCname~methname)+ labs(title = paste0("filter_by_",meth[n],"_theoretical_vs_detected_fdr")))
#       
#     }
#   }
#   
# }
m=1
a<-read.table(paste0(out_path[m],"/",meth[m],"_sim3_GMvsGM_power"))

