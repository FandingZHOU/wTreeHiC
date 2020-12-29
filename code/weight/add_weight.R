#weight methods

dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
meth=c("multiHiCcompare","diffHic")
in_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/grouptable","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/grouptable/filter")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/weighttable","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/weighttable/filter")
func_path="/p/keles/fandingzhou/volumeA/DCI/code/weight/"
group<-c("A","GM12","GM13","GM14")
seqdepth<-c(1,3,5)
bin_size<-40000
weight_method<-c("xi","pi0")
scalesize<-c(2,3,5,10)

library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(splines)
library(qvalue,lib.loc = dict_path)

#source(paste0(func_path,"weight_function.R"))

#regular group for GM vs GM spike-in
for(m in 1:2){
  
  for (n in 1:3){
    for(k in 1:2){
      HiCbinpairs_data<-read_tsv(paste0(in_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_group.tsv"))
      HiCbinpairs_data<-add_weight(HiCbinpairs_data,weight_method[k])
      write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
      }
  }
}

#regular group for GM vs GM spike-in multiple scales

for(ss in 1:4){
  for(m in 1:2){
    for (n in 1:3){
      for(k in 1:2){
        HiCbinpairs_data<-read_tsv(paste0(in_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_group.tsv"))
        HiCbinpairs_data<-add_weight(HiCbinpairs_data,weight_method[k])
        write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
      }
    }
  }
}



#small group for GM vs GM spike-in --xi
k=1
for(m in 1:2){
  for (n in 1:3){
    HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"group_small.tsv"))
    HiCbinpairs_data<-add_weight(HiCbinpairs_data,weight_method[k])
    write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],"_small.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
  }
}

#large group for GM vs GM spike-in --pi0
k=2
for(m in 1:2){
  for (n in 1:3){
    HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_group_large.tsv"))
    HiCbinpairs_data<-add_weight(HiCbinpairs_data,weight_method[k])
    write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],"_large.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
  }
}

#regular GM vs A and GM vs GM
for(m in 1:2){
  for (n in 1:4){
    for (k in 1:2){
      HiCbinpairs_data<-read_tsv(paste0(in_path[m],"/",meth[m],"_GMvs",group[n],"_group.tsv")) 
      HiCbinpairs_data<-add_weight(HiCbinpairs_data,weight_method[k])
      write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_weightedp_",weight_method[k],".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
    }
  }
}


#HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_group_large.tsv"))
#HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_group_small.tsv"))
#write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_weightedp_",weight_method,"_large.tsv"), sep="\t",quote=FALSE, row.names=FALSE)






