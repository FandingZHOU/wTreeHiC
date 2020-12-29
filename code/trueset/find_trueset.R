options("scipen"=100)

dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
meth=c("multiHiCcompare","diffHic")
in_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/grouptable","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/grouptable")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/filter")
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
TSetKey<-scan(paste0(out_path[2],"/TSetKey"),what = character())
#load("/p/keles/treehic/volumeD/freeHiC/GM12878_A549_TSetKey_trueSpikeinSmooth.RData")
load("/p/keles/treehic/volumeD/freeHiC/GM12878_A549_dataSign_trueSpikein.RData")
Truesig<-paste0(dataSign$binA,"_",dataSign$binB)
TSetKey2<-unique(smooth_sig)
length(intersect(TSetKey,TSetKey2))
#new TSetKey
# tset<-strsplit(x=TSetKey, split="_")
# linkpair<-function(x){
#   x<-as.numeric(x)
#   paste0(x[1],"_",x[2])
# }
# TSetKey<-sapply(tset, linkpair)
# length(grep(pattern = "e", x = TSetKey))
# write(TSetKey, file=paste0(out_path[2],"/TSetKey"))
length(scan(paste0(out_path[2],"/TSetKey"),what = character()))
# build a matrix for smoothed and non-smooth true set
truesignal<-function(bins,bin_size){
  binA=bins[1]
  binB=bins[2]
  smoothsig = c()
  k=1
  for(i in -2:2){
    for(j in -2:2){
      if(abs(i)+abs(j)<=2){
        smoothsig[k] = paste0(binA+i*bin_size,"_",binB+j*bin_size)
        k=k+1
      }
    }
  }
  return(smoothsig)
}
Truesig2 = cbind(dataSign$binA,dataSign$binB)
smooth_sig<-c(apply(Truesig2,MARGIN = 1,truesignal,bin_size=bin_size))
orig_sig<-c(sapply(Truesig,rep,times=13))
Truesigmatrix<-data.frame(smooth=smooth_sig,origin=orig_sig,stringsAsFactors = FALSE)
write.csv(Truesigmatrix, file=paste0(out_path[2],"/TSetmatix"))
head(smooth_sig)
length(unique(smooth_sig))
#Truesig2=c("540000_740000","540000_860000")
#Truesig2=matrix(c(540000,540000, 740000,860000),2,2)

#regular group for GM vs GM spike-in
Truesigmatrix<-read.csv(paste0(out_path[2],"/TSetmatix"),stringsAsFactors=F)
Trueorig<-unique(Truesigmatrix$origin)

keep_rate<-c()
i=0
for (n in 1:3){
  binpairs<-list()
  for(m in 1:2){
      HiCbinpairs_data<-read_tsv(paste0(in_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_group.tsv"))
      binpairs[[m]]<-paste0(HiCbinpairs_data$bin_1,"_",HiCbinpairs_data$bin_2)
  }
  a<-intersect(binpairs[[1]],binpairs[[2]])
  keep_sig = Truesigmatrix[Truesigmatrix$smooth %in% a,]
  i=i+1
  keep_rate[i]<-length(unique(keep_sig$origin))/length(Trueorig)
  write.csv(keep_sig, file=paste0(out_path[2],"/GMvsvsGMspikein_seqx",seqdepth[n],"TSetmatrix"),row.names =F)
}

# keep_rate<-c()
# i=0
# for (n in 1:3){
#   binpairs<-list()
#   for(m in 1:2){
#       HiCbinpairs_data<-read_tsv(paste0(in_path[m],"/",meth[m],"_GMvsvsGMspikein_seqx",seqdepth[n],"_group.tsv"))
#       binpairs[[m]]<-paste0(HiCbinpairs_data$bin_1,"_",HiCbinpairs_data$bin_2) 
#   }
#   a<-intersect(binpairs[[1]],binpairs[[2]])
#   nosmooth<-intersect(Truesig,a)
#   smooth<-intersect(TSetKey,a)
#   i=i+1
#   keep_rate[i]<-length(smooth)/length(TSetKey)
#   write(nosmooth, file=paste0(out_path[2],"/_GMvsvsGMspikein_seqx",seqdepth[n],"_TSet"))
#   write(smooth, file=paste0(out_path[2],"/_GMvsvsGMspikein_seqx",seqdepth[n],"_smoothTSet"))
#   write(a, file=paste0(out_path[2],"/_GMvsvsGMspikein_seqx",seqdepth[n],"_keep"))
# }


#regular group for GM vs GM spike-in multiple scales
for(ss in 1:4){
for (n in 1:3){
  binpairs<-list()
  for(m in 1:2){
    HiCbinpairs_data<-read_tsv(paste0(in_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_group.tsv"))
    binpairs[[m]]<-paste0(HiCbinpairs_data$bin_1,"_",HiCbinpairs_data$bin_2)
  }
  a<-intersect(binpairs[[1]],binpairs[[2]])
  keep_sig = Truesigmatrix[Truesigmatrix$smooth %in% a,]
  i=i+1
  keep_rate[i]<-length(unique(keep_sig$origin))/length(Trueorig)
  write.csv(keep_sig, file=paste0(out_path[2],"/GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"TSetmatrix"),row.names =F)
}
}
# for(ss in 1:4){
#   for (n in 1:3){
#     binpairs<-list()
#     for(m in 1:2){
#       HiCbinpairs_data<-read_tsv(paste0(in_path[m],"/",meth[m],"_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_group.tsv"))
#       binpairs[[m]]<-paste0(HiCbinpairs_data$bin_1,"_",HiCbinpairs_data$bin_2)
#     }
#     a<-intersect(binpairs[[1]],binpairs[[2]])
#     nosmooth<-intersect(Truesig,a)
#     smooth<-intersect(TSetKey,a)
#     i=i+1
#     keep_rate[i]<-length(smooth)/length(TSetKey)
#     # print(length(binpairs[[1]]))
#     # print(length(binpairs[[2]]))
#     # print(length(a))
#     write(nosmooth, file=paste0(out_path[2],"/_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_TSet"))
#     write(smooth, file=paste0(out_path[2],"/_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_smoothTSet"))
#     write(a, file=paste0(out_path[2],"/_GMvsvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],"_keep"))
#   }
# }
print(keep_rate)
keep_rate<-c( 0.4656710, 0.7551969, 0.8586402, 0.6002461, 0.7639633, 0.8581394, 0.6148333,
              0.7410028, 0.8336012, 0.6162549, 0.7069121, 0.7900630, 0.6120413, 0.6645182,
              0.7181096)
#keep_rate<-c(0.09561097,0.26433753, 0.39779556, 0.14127637, 0.28984463, 0.40834624, 0.15576045, 0.28235451, 0.38367442, 0.16707496, 0.26877395, 0.34694041, 0.17605882, 0.25141179, 0.30164320)
keep_truesig<-matrix(keep_rate,3,5)
rownames(keep_truesig)<-paste0("seqx",c(1,3,5))
colnames(keep_truesig)<-paste0("scale",c(1,2,3,5,10))
keep_truesig
c(2476996, 2949200)