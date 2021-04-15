#weight methods
in_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/grouptable",
          "/p/keles/fandingzhou/volumeA/DCI/diffHic/data/grouptable/filter",
          "/p/keles/fandingzhou/volumeA/DCI/HiCDCPlus/data/grouptable/filter")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/weighttable",
           "/p/keles/fandingzhou/volumeA/DCI/diffHic/data/weighttable/filter",
           "/p/keles/fandingzhou/volumeA/DCI/HiCDCPlus/data/weighttable/filter")

func_path="/p/keles/fandingzhou/volumeA/DCI/code/weight/"
meth=c("multiHiCcompare","diffHic","HiCDCPlus")

dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
group<-c("A","GM12","GM13","GM14")
seqdepth<-c(1,3,5)
bin_size<-40000
weight_method<-c("xi","pi0")
scalesize<-c(2,3,5,10)
#install.packages('utf8',lib = dict_path)#reshape2,fansi'
library('utf8',lib.loc = dict_path)
library(cli,lib.loc = dict_path)
library(fansi,lib.loc = dict_path)
library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(splines)
library(qvalue,lib.loc = dict_path)
library(IHW,lib.loc = dict_path)

#source(paste0(func_path,"weight_function.R"))

for(seed in 11:20){
  #seed=4
  for(test in c("4vs4","2vs2")){
    #if(seed!=1|test!="4vs4"){
      for(n in 1:2){
        for(m in 1:3){
          in_path2=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/",test,"/filter_as_",meth[n],"/",meth[m])
          for (foldchange in c(2,3,4,6)){
            if(!file.exists(paste0(in_path2,"/weighttable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[2],"_distance.tsv"))){
              
            #in_path2=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/filter_as_",meth[n],"/",meth[m])
             #HiCbinpairs_data<-read_tsv(paste0(in_path2,"/grouptable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_group.tsv"))
              HiCbinpairs_data<-read_tsv(paste0(in_path2,"/grouptable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_group_distance.tsv"))
              for(k in 1:2){
             #k=2
              if(m==3){
                HiCbinpairs_data2<-add_weight(HiCbinpairs_data,weight_method[k],statistic="Z")
              }else{
                HiCbinpairs_data2<-add_weight(HiCbinpairs_data,weight_method[k],statistic="F")
              }
              #write.table(HiCbinpairs_data2, file=paste0(in_path2,"/weighttable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
                write.table(HiCbinpairs_data2, file=paste0(in_path2,"/weighttable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],"_distance.tsv"), sep="\t",quote=FALSE, row.names=FALSE)
                
              }}
          }
        }
      }
    }
    
  }
#}



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

#regular group for GM vs GM spike-in
#for(m in 1:2){
m=2
for (n in 1:3){
  for(k in 1:2){
    HiCbinpairs_data<-read_tsv(paste0(in_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_group.tsv"))
    HiCbinpairs_data<-add_weight(HiCbinpairs_data,weight_method[k],statistic="F")
    write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
  }
}
#}

m=3
for (n in 1:3){
  for(k in 1:2){
    HiCbinpairs_data<-read_tsv(paste0(in_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_group.tsv"))
    HiCbinpairs_data<-add_weight(HiCbinpairs_data,weight_method[k],statistic="Z")
    write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
  }
}

#regular group for GM vs GM spike-in multiple scales
m=2
for(ss in 1:4){
  #for(m in 1:2){
  for (n in 1:3){
    for(k in 1:2){
      HiCbinpairs_data<-read_tsv(paste0(in_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_group.tsv"))
      HiCbinpairs_data<-add_weight(HiCbinpairs_data,weight_method[k])
      write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_weightedp_",weight_method[k],".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
    }
  }
  #}
}

#regular group for GM vs GM new sim
in_path[1]=paste0(in_path[1],"/simulation/")
in_path[2]="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/grouptable/simulation/filter/"
out_path[1]=paste0(out_path[1],"/simulation/")
out_path[2]="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/weighttable/simulation/filter/"
in_path[3]=paste0(in_path[3],"/simulation/")
out_path[3]=paste0(out_path[3],"/simulation/")

for(m in 1:2){
  #m=3
  for(k in 1:2){
    for(foldchange in c(3,4,6)){
      HiCbinpairs_data<-read_tsv(paste0(in_path[m],meth[m],"_sim3_binpairs_GMvsGMFC",foldchange,"_group.tsv"))
      HiCbinpairs_data<-add_weight(HiCbinpairs_data,weight_method[k],statistic="Z")
      write.table(HiCbinpairs_data, file=paste0(out_path[m],meth[m],"_sim3_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],".tsv"), sep="\t",quote=FALSE, row.names=FALSE)
    }
  }
}
n=m=k=1
foldchange=6

#HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_group_large.tsv"))
#HiCbinpairs_data<-read_tsv(paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_group_small.tsv"))
#write.table(HiCbinpairs_data, file=paste0(out_path[m],"/",meth[m],"_GMvs",group[n],"_weightedp_",weight_method,"_large.tsv"), sep="\t",quote=FALSE, row.names=FALSE)






