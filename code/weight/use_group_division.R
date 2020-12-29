dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path=c("/p/keles/fandingzhou/volumeA/DCI/ABcompartment")
in_path1=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/originaltable","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/originaltable/filter")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/grouptable","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/grouptable/filter")
func_path="/p/keles/fandingzhou/volumeA/DCI/code/weight/"
meth=c("multiHiCcompare","diffHic")
seqdepth<-c(1,3,5)
group<-c("A","GM12","GM13","GM14")
doc_name<-c("_Multi_glm","_diffHic.tsv")
group_size<-c("regular","large","small")
scalesize<-c(2,3,5,10)

library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 


#source(paste0(func_path,"group_division_function.R"))

bin_size<-40000
ch1_length<-249250621


m=2
##partition for GM vs GM spike-in
#for (m in 1:2){
  for (n in 1:3){
    HiCbinpairs<-read_tsv(paste0(in_path1[m],"/binpairs_GMvsGMspikein_seqx",seqdepth[n],doc_name[m]),col_names=TRUE)
    if(m==1){
      HiCbinpairs_data<-cbind(HiCbinpairs$region1+bin_size/2,HiCbinpairs$region2+bin_size/2,HiCbinpairs$F,HiCbinpairs$p.value,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
    }else{
      HiCbinpairs_data<-cbind(HiCbinpairs$end2-bin_size/2,HiCbinpairs$end1-bin_size/2,HiCbinpairs$F,HiCbinpairs$PValue,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
    }
    colnames(HiCbinpairs_data)<-c("bin_1","bin_2","F","P","Group","weight")
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data[,1]),])
    boundary1<-read.table(paste0(in_path,"/GM_spikein_boundaries"))
    boundary1<- boundary1[,1]
    boundary1<-c(boundary1,ch1_length)
    #for(s in 1:3){
    s=1
      HiCbinpairs_data<-group_division(HiCbinpairs_data,boundary1,size=group_size[s])
      if(s==1){
        write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_group.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
      }else{
        write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_group_",group_size[s],".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
      }
    }
  #}
#}


##partition for GM vs GM spike-in scales
m=2
for(ss in 1:4){
#for (m in 1:2){
for (n in 1:3){
  HiCbinpairs<-read_tsv(paste0(in_path1[m],"/binpairs_GMvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],doc_name[m]),col_names=TRUE)
  if(m==1){
    HiCbinpairs_data<-cbind(HiCbinpairs$region1+bin_size/2,HiCbinpairs$region2+bin_size/2,HiCbinpairs$F,HiCbinpairs$p.value,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
  }else{
    HiCbinpairs_data<-cbind(HiCbinpairs$end2-bin_size/2,HiCbinpairs$end1-bin_size/2,HiCbinpairs$F,HiCbinpairs$PValue,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
  }
  colnames(HiCbinpairs_data)<-c("bin_1","bin_2","F","P","Group","weight")
  HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data[,1]),])
  boundary1<-read.table(paste0(in_path,"/GM_spikein_boundaries"))
  boundary1<- boundary1[,1]
  boundary1<-c(boundary1,ch1_length)
  s=1
  HiCbinpairs_data<-group_division(HiCbinpairs_data,boundary1,size=group_size[s])
  write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_group.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
}

#}
}

#partition for GM vs A and GM vs GM
for (m in 1:2){
  for (n in 1:4){
    HiCbinpairs<-read_tsv(paste0(in_path1[m],"/binpairs_GMvs",group[n],doc_name[m]),col_names=TRUE)
    if(m==1){
      HiCbinpairs_data<-cbind(HiCbinpairs$region1+bin_size/2,HiCbinpairs$region2+bin_size/2,HiCbinpairs$F,HiCbinpairs$p.value,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
    }else{
      HiCbinpairs_data<-cbind(HiCbinpairs$end2-bin_size/2,HiCbinpairs$end1-bin_size/2,HiCbinpairs$F,HiCbinpairs$PValue,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
    }
    colnames(HiCbinpairs_data)<-c("bin_1","bin_2","F","P","Group","weight")
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data[,1]),])
    if(n==1){
      boundary1<-read.table(paste0(in_path,"/GM_A_boundaries"))
    }else{
      boundary1<-read.table(paste0(in_path,"/GM_GM_boundaries"))
    }
    
    boundary1<- boundary1[,1]
    boundary1<-c(boundary1,ch1_length)
    #for(s in 1:3){
    s=1
      HiCbinpairs_data<-group_division(HiCbinpairs_data,boundary1,size=group_size[s])
      if(s==1){
        write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_group.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
      }else{
        write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvs",seqdepth[n],"_group_",group_size[s],".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
      }
    #}
  }
}