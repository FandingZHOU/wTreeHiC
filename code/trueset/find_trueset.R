options("scipen"=100)

dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
meth=c("multiHiCcompare","diffHic","HiCDCPlus")
in_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/grouptable","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/grouptable")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/filter")
func_path="/p/keles/fandingzhou/volumeA/DCI/code/weight/"
group<-c("A","GM12","GM13","GM14")
seqdepth<-c(1,3,5)
bin_size<-40000
weight_method<-c("xi","pi0")
scalesize<-c(2,3,5,10)

library(cli,lib.loc = dict_path)
library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
TSetKey<-scan(paste0(out_path[2],"/TSetKey"),what = character())
#load("/p/keles/treehic/volumeD/freeHiC/GM12878_A549_TSetKey_trueSpikeinSmooth.RData")
load("/p/keles/treehic/volumeD/freeHiC/GM12878_A549_dataSign_trueSpikein.RData")
Truesig<-paste0(dataSign$binA,"_",dataSign$binB)

#new TSetKey
# tset<-strsplit(x=TSetKey, split="_")
# linkpair<-function(x){
#   x<-as.numeric(x)
#   paste0(x[1],"_",x[2])
# }
# TSetKey<-sapply(tset, linkpair)
# length(grep(pattern = "e", x = TSetKey))
# write(TSetKey, file=paste0(out_path[2],"/TSetKey"))

#regular group for GM vs GM spike-in
for (n in 1:3){
  binpairs<-list()
  for(m in 1:2){
      HiCbinpairs_data<-read_tsv(paste0(in_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_group.tsv"))
      binpairs[[m]]<-paste0(HiCbinpairs_data$bin_1,"_",HiCbinpairs_data$bin_2) 
  }
  a<-intersect(binpairs[[1]],binpairs[[2]])
  nosmooth<-intersect(Truesig,a)
  smooth<-intersect(TSetKey,a)
  write(nosmooth, file=paste0(out_path[2],"/_GMvsvsGMspikein_seqx",seqdepth[n],"_TSet"))
  write(smooth, file=paste0(out_path[2],"/_GMvsvsGMspikein_seqx",seqdepth[n],"_smoothTSet"))
  write(a, file=paste0(out_path[2],"/_GMvsvsGMspikein_seqx",seqdepth[n],"_keep"))
}

#regular group for GM vs GM spike-in multiple scales
for(ss in 1:4){
  for (n in 1:3){
    binpairs<-list()
    for(m in 1:2){
      HiCbinpairs_data<-read_tsv(paste0(in_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_group.tsv"))
      binpairs[[m]]<-paste0(HiCbinpairs_data$bin_1,"_",HiCbinpairs_data$bin_2)
    }
    a<-intersect(binpairs[[1]],binpairs[[2]])
    nosmooth<-intersect(Truesig,a)
    smooth<-intersect(TSetKey,a)
    write(nosmooth, file=paste0(out_path[2],"/_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_TSet"))
    write(smooth, file=paste0(out_path[2],"/_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_smoothTSet"))
    write(a, file=paste0(out_path[2],"/_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_keep"))
  }
}

in_path[3]="/p/keles/fandingzhou/volumeA/DCI/simulation/"
out_path[3]="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/filter/simulation"
TSetKey<-scan(paste0(in_path[3],"True_signal5"),what = character())
#TSetKey<-scan(paste0(in_path,"True_signal5"),what = character())
j=0
ratiotable<-c()
Keeptable<-c()
for(n in 1:3){
  for(foldchange in c(2,3,4,6)){
    j=j+1
    in_path4=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/filter_as_",meth[n])
    Truesig=scan(paste0(in_path4,"/Trueset/GMvsGMFC",foldchange,"_TSet"),what = character())
    Keepsig=scan(paste0(in_path4,"/Trueset/GMvsGMFC",foldchange,"_keep"),what = character())
    ratiotable[j]=length(Truesig)/length(TSetKey)
    Keeptable[j]=length(Keepsig)
  }
}
ratiotable=matrix(ratiotable,4,3)
colnames(ratiotable)=meth
row.names(ratiotable)=paste0("FC",c(2,3,4,6))
Keeptable=matrix(Keeptable,4,3)
colnames(Keeptable)=meth
row.names(Keeptable)=paste0("FC",c(2,3,4,6))
#regular group for GM vs new GM simulation
for (foldchange in c(2,3,4,6)){
  binpairs<-list()
  for(m in 1:2){
    HiCbinpairs_data<-read_tsv(paste0(in_path[m],"/simulation/",meth[m],"_sim3_binpairs_GMvsGMFC",foldchange,"_group.tsv"))
    binpairs[[m]]<-paste0(HiCbinpairs_data$bin_1,"_",HiCbinpairs_data$bin_2) 
  }
  a<-intersect(binpairs[[1]],binpairs[[2]])
  nosmooth<-intersect(TSetKey,a)
  print(length(a))
  print(length(binpairs[[1]]))
  print(length(binpairs[[2]]))
  print(length(nosmooth)/length(TSetKey))

  write(nosmooth, file=paste0(out_path[3],"/GMvsGMFC",foldchange,"_TSet"))
  write(a, file=paste0(out_path[3],"/GMvsGMFC",foldchange,"_keep"))
}

#filter as diffHic

in_path[3]="/p/keles/fandingzhou/volumeA/DCI/filter/filter_as_diffHic/originaltable"
out_path[3]="/p/keles/fandingzhou/volumeA/DCI/filter/filter_as_diffHic/Trueset"
TSetKey<-scan(paste0(in_path[3],"True_signal5"),what = character())
#regular group for GM vs new GM simulation
for (foldchange in c(2,3,4,6)){
  binpairs<-list()
  HiCbinpairs_data<-read_tsv(paste0(in_path[2],"/simulation/",meth[2],"_sim3_binpairs_GMvsGMFC",foldchange,"_group.tsv"))
  binpairs<-paste0(HiCbinpairs_data$bin_1,"_",HiCbinpairs_data$bin_2) 
  nosmooth<-intersect(TSetKey,binpairs)
  print(length(nosmooth))
  write(nosmooth, file=paste0(out_path[3],"/GMvsGMFC",foldchange,"_TSet"))
  write(binpairs, file=paste0(out_path[3],"/GMvsGMFC",foldchange,"_keep"))
}

#filter as multiHiCcompare
in_path4="/p/keles/fandingzhou/volumeA/DCI/simulation/"
TSetKey<-scan(paste0(in_path4,"True_signal01"),what = character())
in_path2<-paste0("/p/keles/fandingzhou/volumeA/DCI/filter/filter_as_",c("multiHiCcompare","diffHic","HiCDCPlus"))
k=1
kpset=c()
trset=c()
for(foldchange in c(2,3,4,6)){
  HiCbinpairs<-read_tsv(paste0(in_path2[k],"/multiHiCcompare/originaltable/binpairs_GMvsGMFC",foldchange,"_Multi_glm"))
  binpairs<-paste0(HiCbinpairs$region1+bin_size/2,"_",HiCbinpairs$region2+bin_size/2)
  nosmooth<-intersect(TSetKey,binpairs)
  kpset<-c(kpset,length(binpairs))
  trset<-c(trset,length(nosmooth)/length(TSetKey))
  write(nosmooth, file=paste0(in_path2[k],"/Trueset/GMvsGMFC",foldchange,"_TSet"))
  write(binpairs, file=paste0(in_path2[k],"/Trueset/GMvsGMFC",foldchange,"_keep"))
}
print(kpset,trset)


#filter as HiCDCPlus
out_path[4]="/p/keles/fandingzhou/volumeA/DCI/filter/filter_as_HiCDCPlus/Trueset"
for (foldchange in c(2,3,4,6)){
    in_path[3]="/p/keles/fandingzhou/volumeA/DCI/HiCDCPlus/data/hicpromatrix"
    HiCbinpairs_data<-data.table::fread(paste0(in_path[3],"/simulation/sim3_FC",foldchange,'_GM6_indices.txt.gz'))
    binpairs<-paste0(HiCbinpairs_data$startI+bin_size/2,"_",HiCbinpairs_data$startJ+bin_size/2)
    nosmooth<-intersect(TSetKey,binpairs)
    print(length(nosmooth))
    write(nosmooth, file=paste0(out_path[4],"/GMvsGMFC",foldchange,"_TSet"))
    write(binpairs, file=paste0(out_path[4],"/GMvsGMFC",foldchange,"_keep"))
}
for (foldchange in c(2,3,4,6)){
  in_path2=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/filter_as_HiCDCPlus/HiCDCPlus")
  HiCbinpairs<-results(readRDS(paste0(in_path2,"/originaltable/GM_FC",foldchange,"/chr1_DESeq2_obj.rds")))
  tset<-strsplit(x=rownames(HiCbinpairs), split=":")
  binpairs<-paste0(sapply(tset, binA,bin_size,T),"_",sapply(tset, binB,bin_size,T))
  nosmooth<-intersect(TSetKey,binpairs)
  print(length(nosmooth)/length(TSetKey))
}
a<-sapply(tset, binA,bin_size,T)
b<-sapply(tset, binB,bin_size,T)