dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path=c("/p/keles/fandingzhou/volumeA/DCI/ABcompartment")
in_path1=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/originaltable",
          # "/p/keles/fandingzhou/volumeA/DCI/diffHic/data/originaltable/filter",
           "/p/keles/fandingzhou/volumeA/DCI/diffHic/data/originaltable",
           "/p/keles/fandingzhou/volumeA/DCI/HiCDCPlus/data/filter")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/grouptable",
           #"/p/keles/fandingzhou/volumeA/DCI/diffHic/data/grouptable/filter",
           "/p/keles/fandingzhou/volumeA/DCI/diffHic/data/grouptable",
           "/p/keles/fandingzhou/volumeA/DCI/HiCDCPlus/data/grouptable/filter")
func_path="/p/keles/fandingzhou/volumeA/DCI/code/weight/"
meth=c("multiHiCcompare","diffHic","HiCDCPlus")
seqdepth<-c(1,3,5)
group<-c("A","GM12","GM13","GM14")
doc_name<-c("_Multi_glm","_diffHic.tsv")
group_size<-c("regular","large","small")
scalesize<-c(2,3,5,10)

library(utf8,lib.loc = dict_path)
library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(DESeq2,lib.loc = dict_path)# 

binA<-function(x,bin_size,tomid=F){
  if(tomid==F){
    return(as.numeric(x[1])-bin_size/2)
  }else{
    return(as.numeric(x[1])+bin_size/2)
  }
}
binB<-function(x,bin_size,tomid=F){
  if(tomid==F){
    return(as.numeric(x[2])-bin_size/2)
  }else{
    return(as.numeric(x[2])+bin_size/2)
  }
}

#source(paste0(func_path,"group_division_function.R"))

bin_size<-40000
ch1_length<-249250621


m=3
##partition for GM vs GM spike-in
#for (m in 1:2){
  for (n in 1:3){
    if(m!=3){
      HiCbinpairs<-read_tsv(paste0(in_path1[m],"/binpairs_GMvsGMspikein_seqx",seqdepth[n],doc_name[m]),col_names=TRUE)
      if(m==1){
        HiCbinpairs_data<-cbind(HiCbinpairs$region1+bin_size/2,HiCbinpairs$region2+bin_size/2,HiCbinpairs$F,HiCbinpairs$p.value,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
      }
      if(m==2){
        HiCbinpairs_data<-cbind(HiCbinpairs$end2-bin_size/2,HiCbinpairs$end1-bin_size/2,HiCbinpairs$F,HiCbinpairs$PValue,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
      }
      colnames(HiCbinpairs_data)<-c("bin_1","bin_2","F","P","Group","weight")
    }else{
      HiCbinpairs<-results(readRDS(paste0(in_path1[m],"/GM_spikein_seqx",seqdepth[n],"/chr1_DESeq2_obj.rds")))
      #HiCbinpairs<-results(readRDS(paste0(out_path,"GM_spikein_seqx",seqdepth[k],"/chr1_DESeq2_obj.rds")))
      tset<-strsplit(x=rownames(HiCbinpairs), split=":")
      HiCbinpairs_data<-cbind(sapply(tset, binA,bin_size,T),sapply(tset, binB,bin_size,T),HiCbinpairs$stat,HiCbinpairs$pvalue,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
      colnames(HiCbinpairs_data)<-c("bin_1","bin_2","Z","P","Group","weight")
    }
    
    
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
    #}
  }
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

#partition for GM vs new simulation

for (m in 1:2){
  for (foldchange in c(2,3,4,6)){
    HiCbinpairs<-read_tsv(paste0(in_path1[m],"/simulation/sim3_binpairs_GMvsGMFC",foldchange,doc_name[m]),col_names=TRUE)
    if(m==1){
      HiCbinpairs_data<-cbind(HiCbinpairs$region1+bin_size/2,HiCbinpairs$region2+bin_size/2,HiCbinpairs$F,HiCbinpairs$p.value,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
    }else{
      HiCbinpairs_data<-cbind(HiCbinpairs$end2-bin_size/2,HiCbinpairs$end1-bin_size/2,HiCbinpairs$F,HiCbinpairs$PValue,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
    }
    colnames(HiCbinpairs_data)<-c("bin_1","bin_2","F","P","Group","weight")
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data[,1]),])
    boundary1<-read.table(paste0(in_path,"/GM_GM_boundaries"))
    boundary1<- boundary1[,1]
    boundary1<-c(boundary1,ch1_length)
    s=1
    HiCbinpairs_data<-group_division(HiCbinpairs_data,boundary1,size=group_size[s])
    write.table(HiCbinpairs_data, file=paste0(out_path[m],"/simulation/",meth[m],"_sim3_binpairs_GMvsGMFC",foldchange,"_group.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
  }
}

#partition for GM vs new simulation new filtering
m=2
  for (foldchange in c(2,3,4,6)){
     if(m!=3){
       HiCbinpairs<-read_tsv(paste0(in_path1[m],"/simulation/filter/sim3_binpairs_GMvsGMFC",foldchange,doc_name[m]),col_names=TRUE)
      if(m==1){
        HiCbinpairs_data<-cbind(HiCbinpairs$region1+bin_size/2,HiCbinpairs$region2+bin_size/2,HiCbinpairs$F,HiCbinpairs$p.value,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
      }else{
        HiCbinpairs_data<-cbind(HiCbinpairs$end2-bin_size/2,HiCbinpairs$end1-bin_size/2,HiCbinpairs$F,HiCbinpairs$PValue,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
      }
      colnames(HiCbinpairs_data)<-c("bin_1","bin_2","F","P","Group","weight")
    }else{
     
      HiCbinpairs<-results(readRDS(paste0(in_path1[3],"/simulation/GM_filter_FC",foldchange,"/chr1_DESeq2_obj.rds")))

      tset<-strsplit(x=rownames(HiCbinpairs), split=":")
      HiCbinpairs_data<-cbind(sapply(tset, binA,bin_size,T),sapply(tset, binB,bin_size,T),HiCbinpairs$stat,HiCbinpairs$pvalue,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
      colnames(HiCbinpairs_data)<-c("bin_1","bin_2","Z","P","Group","weight")
    }
    
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data[,1]),])
    boundary1<-read.table(paste0(in_path,"/GM_GM_boundaries"))
    boundary1<- boundary1[,1]
    boundary1<-c(boundary1,ch1_length)
    s=1
    HiCbinpairs_data<-group_division(HiCbinpairs_data,boundary1,size=group_size[s])
    #simulation
    #write.table(HiCbinpairs_data, file=paste0(out_path[m],"/simulation/filter/",meth[m],"_sim3_binpairs_GMvsGMFC",foldchange,"_group.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
    write.table(HiCbinpairs_data, file=paste0(out_path[m],"/simulation/",meth[m],"_sim3_binpairs_GMvsGMFC",foldchange,"_group.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
    
  }

seed=12
test="4vs4"
for(seed in 11:20){
  boundary1<-read.table(paste0(in_path,"/seed",seed,"GM_FC_boundaries"))
  boundary1<- boundary1[,1]
  boundary1<-c(boundary1,ch1_length)
  for(test in c("4vs4","2vs2")){
    #if(seed !=8 |test!="2vs2"){
      for(n in 1:2){
        for(m in 1:3){
          for (foldchange in c(2,3,4,6)){
            #in_path2=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/filter_as_",meth[n],"/",meth[m])
            in_path2=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/",test,"/filter_as_",meth[n],"/",meth[m])
             #if(!file.exists(paste0(in_path2,"/grouptable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_group.tsv"))){
          #     print(seed)
          #     print(test)
          #     print(n)
          #     print(m)
          #     print(foldchange)
          #   }
          # }}}}}
            if(m!=3){
              HiCbinpairs<-read_tsv(paste0(in_path2,"/originaltable/binpairs_GMvsGMFC",foldchange,doc_name[m]),col_names=TRUE)
              if(m==1){
                HiCbinpairs_data<-cbind(HiCbinpairs$region1+bin_size/2,HiCbinpairs$region2+bin_size/2,HiCbinpairs$F,HiCbinpairs$p.value,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
              }else{
                HiCbinpairs_data<-cbind(HiCbinpairs$end2-bin_size/2,HiCbinpairs$end1-bin_size/2,HiCbinpairs$F,HiCbinpairs$PValue,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
              }
              colnames(HiCbinpairs_data)<-c("bin_1","bin_2","F","P","Group","weight")
            }else{
              HiCbinpairs<-results(readRDS(paste0(in_path2,"/originaltable/GM_FC",foldchange,"/chr1_DESeq2_obj.rds")))
              tset<-strsplit(x=rownames(HiCbinpairs), split=":")
              HiCbinpairs_data<-cbind(sapply(tset, binA,bin_size,T),sapply(tset, binB,bin_size,T),HiCbinpairs$stat,HiCbinpairs$pvalue,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
              colnames(HiCbinpairs_data)<-c("bin_1","bin_2","Z","P","Group","weight")
            }
            
            HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data[,1]),])
            s=1
            #HiCbinpairs_data<-group_division(HiCbinpairs_data,boundary1,clustern=60,size=group_size[s])
            ##Group by distances
            HiCbinpairs_data<-group_division(HiCbinpairs_data,boundary1,groups.by="distance",size=group_size[s])
            
            #simulation
            #write.table(HiCbinpairs_data, file=paste0(out_path[m],"/simulation/filter/",meth[m],"_sim3_binpairs_GMvsGMFC",foldchange,"_group.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
            #write.table(HiCbinpairs_data, file=paste0(in_path2,"/grouptable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_group.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
            write.table(HiCbinpairs_data, file=paste0(in_path2,"/grouptable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_group_distance.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
            
          }
        }
        
      }
    }

    
  }
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