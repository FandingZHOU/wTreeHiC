# install.packages(c('matrixStats','memoise','digest'),lib = dict_path)
# install.packages(c('foreach','HiCcompare'),lib = dict_path)
install.packages(c('checkmate','backports'),lib = dict_path)
# remove.packages('HiCcompare', lib= dict_path)

dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
func_path="/p/keles/fandingzhou/volumeA/DCI/code/multiHiCcompare"
in_path1="/p/keles/fandingzhou/volumeA/DCI/simulation/"
out_path="/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/originaltable"

##NOTE:Sim2-->strategy4;Sim3-->strategy5
library(cli,lib.loc = dict_path)
library(Rcpp,lib.loc = dict_path)
library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(limma,lib.loc = dict_path)
library(edgeR,lib.loc = dict_path)# 
library(BiocParallel,lib.loc = dict_path)# 
library(HiCcompare,lib.loc = dict_path)# 

library(BiocGenerics,lib.loc = dict_path)
library(Biobase,lib.loc = dict_path) 
library(S4Vectors,lib.loc = dict_path)  
library(IRanges,lib.loc = dict_path)
library(AnnotationDbi,lib.loc = dict_path)
library(org.Hs.eg.db,lib.loc = dict_path)
library(hgu133a.db,lib.loc = dict_path)

library(multiHiCcompare,lib.loc = dict_path)# 
library(pander,lib.loc = dict_path)

source(paste0(func_path,"/multiHiCcompare_function.R"))

seqdepth<-c(1,3,5)
rep<-c(2,3,4,6)
res <- 40000 
scalesize<-c(2,3,5,10)
nrep<-4
# Read data
sample_list <- list()
for(j in 1:4){
  data_j<-read_tsv(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j],"/chr1/s1_training/binPairs/rep",rep[j],"_chr1.binPairs.chr1"),col_names=FALSE)
  #original_list<-cbind(rep(1,nrow(data_j[[j]])),data_j[[j]][,2]-res/2,data_j[[j]][,4]-res/2,data_j[[j]][,5])
  original_list<-cbind(rep(1,nrow(data_j)),data_j[,2]-res/2,data_j[,4]-res/2,data_j[,5])
  sample_list[[j]]<-original_list[abs(original_list[,2]-original_list[,3])>res*2,]
  colnames(sample_list[[j]])<-c('chr','region1','region2','IF')
}
options("scipen"=100)
TSetKey<-scan(paste0(in_path1,"True_signal01"),what = character())
#GM vs GM new simulation

k=2
for(seed in 11:20){
  in_path2<-paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_",c("multiHiCcompare","diffHic","HiCDCPlus"))
  
for(foldchange in c(2,3,4,6)){
   for(j in 1:4){
    data_j<-read_tsv(paste0(in_path1,"GM12878/seed",seed,"_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"))
    original_list<-cbind(rep(1,nrow(data_j)),data_j[,2]-res/2,data_j[,4]-res/2,data_j[,5])
    sample_list[[j+4]]<-original_list[abs(original_list[,2]-original_list[,3])>res*2,]
    colnames(sample_list[[j+4]])<-c('chr','region1','region2','IF')
  }
  filter_pairs<-scan(paste0(in_path2[k],"/Trueset/GMvsGMFC",foldchange,"_keep"),what = character())
  HiCbinpairs<-multiHiCcompare_filter_table(sample_list,filter_pairs,bin_size=res,nrep)
  write_tsv(HiCbinpairs,path =  paste0(in_path2[k],"/multiHiCcompare/originaltable/binpairs_GMvsGMFC",foldchange,"_Multi_glm"))
}
}





###2vs2
nrep<-2
in_path4="/p/keles/fandingzhou/volumeA/DCI/simulation/"

# Read data
sample_list <- list()
##Careful about j+1
for(j in 1:nrep){
  data_j<-read_tsv(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j+1],"/chr1/s1_training/binPairs/rep",rep[j+1],"_chr1.binPairs.chr1"),col_names=FALSE)
  #original_list<-cbind(rep(1,nrow(data_j[[j]])),data_j[[j]][,2]-res/2,data_j[[j]][,4]-res/2,data_j[[j]][,5])
  original_list<-cbind(rep(1,nrow(data_j)),data_j[,2]-res/2,data_j[,4]-res/2,data_j[,5])
  sample_list[[j]]<-original_list[abs(original_list[,2]-original_list[,3])>res*2,]
  colnames(sample_list[[j]])<-c('chr','region1','region2','IF')
}

for(seed in 11:20){
  
  TSetKey<-scan(paste0(in_path4,"True_signal",seed),what = character())
  in_path2<-paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_",c("multiHiCcompare","diffHic","HiCDCPlus"))
  k=1
  kpset=c()
  trset=c()
  for(foldchange in c(2,3,4,6)){
    for(j in 1:nrep){
      data_j<-read_tsv(paste0(in_path1,"GM12878/seed",seed,"_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"))
      original_list<-cbind(rep(1,nrow(data_j)),data_j[,2]-res/2,data_j[,4]-res/2,data_j[,5])
      sample_list[[j+nrep]]<-original_list[abs(original_list[,2]-original_list[,3])>res*2,]
      colnames(sample_list[[j+nrep]])<-c('chr','region1','region2','IF')
    }
    HiCbinpairs<-multiHiCcompare_table(sample_list,nrep,1)
    write_tsv(HiCbinpairs,path =  paste0(in_path2[k],"/multiHiCcompare/originaltable/binpairs_GMvsGMFC",foldchange,"_Multi_glm"))
    bin_size=res
    binpairs<-paste0(HiCbinpairs$region1+bin_size/2,"_",HiCbinpairs$region2+bin_size/2)
    nosmooth<-intersect(TSetKey,binpairs)
    kpset<-c(kpset,length(binpairs))
    trset<-c(trset,length(nosmooth)/length(TSetKey))
    write(nosmooth, file=paste0(in_path2[k],"/Trueset/GMvsGMFC",foldchange,"_TSet"))
    write(binpairs, file=paste0(in_path2[k],"/Trueset/GMvsGMFC",foldchange,"_keep"))
  }
  print(kpset)
  print(trset)
}

in_path5=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_",meth[n])
in_path4=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/",test,"/filter_as_",meth[n])
options("scipen"=100)
TSetKey<-scan(paste0(in_path1,"True_signal01"),what = character())
#GM vs GM new simulation
for(seed in 11:20){
#for(seed in c(1)){
  in_path2<-paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_",c("multiHiCcompare","diffHic","HiCDCPlus"))
  in_path3<-paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/2vs2/filter_as_",c("multiHiCcompare","diffHic","HiCDCPlus"))
  
    for(foldchange in c(2,3,4,6)){
      for(j in 1:nrep){
        data_j<-read_tsv(paste0(in_path1,"GM12878/seed",seed,"_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"))
        original_list<-cbind(rep(1,nrow(data_j)),data_j[,2]-res/2,data_j[,4]-res/2,data_j[,5])
        sample_list[[j+nrep]]<-original_list[abs(original_list[,2]-original_list[,3])>res*2,]
        colnames(sample_list[[j+nrep]])<-c('chr','region1','region2','IF')
      }
      for(k in 1:2){
      filter_pairs<-scan(paste0(in_path2[k],"/Trueset/GMvsGMFC",foldchange,"_keep"),what = character())
      HiCbinpairs<-multiHiCcompare_filter_table(sample_list,filter_pairs,bin_size=res,nrep)
      write_tsv(HiCbinpairs,path =  paste0(in_path3[k],"/multiHiCcompare/originaltable/binpairs_GMvsGMFC",foldchange,"_Multi_glm"))
    }
  }
  
}
  
#}




#}


  
  
  
# GM vs GM spike-in
for(k in 1:3){
  for(j in 1:4){
    data_j<-read_tsv(paste0("/p/keles/treehic/volumeD/freeHiC/GM12878/rep",rep[j],"/chr1/spikeIn3_40kb_smooth_seqx",seqdepth[k],"/simuProcess/s3_binPairs/rep",rep[j],"_chr1.binPairs"),col_names=FALSE)%>% filter(X1 == "chr1", X3 == "chr1")
    original_list<-cbind(rep(1,nrow(data_j)),data_j[,2]-res/2,data_j[,4]-res/2,data_j[,5])
    sample_list[[j+4]]<-original_list[abs(original_list[,2]-original_list[,3])>res*2,]
    colnames(sample_list[[j+4]])<-c('chr','region1','region2','IF')
  }
  HiCbinpairs<-multiHiCcompare_table(sample_list,nrep,Amin)
  write_tsv(HiCbinpairs,path =  paste0(out_path,"/binpairs_GMvsGMspikein_seqx",seqdepth[k],"_Multi_glm"))
}

#GMspikein_scale
for(ss in 1:4){
  for(k in 1:3){
    for(j in 1:4){
      data_j<-read_tsv(paste0("/p/keles/treehic/volumeD/freeHiC/GM12878/rep",rep[j],"/chr1/spikeIn3_40kb_scale",scalesize[ss],"_seqx",seqdepth[k],"/simuProcess/s3_binPairs/rep",rep[j],"_chr1.binPairs"),col_names=FALSE) %>% filter(X1 == "chr1", X3 == "chr1")
      original_list<-cbind(rep(1,nrow(data_j)),data_j[,2]-res/2,data_j[,4]-res/2,data_j[,5])
      sample_list[[j+4]]<-original_list[abs(original_list[,2]-original_list[,3])>res*2,]
      colnames(sample_list[[j+4]])<-c('chr','region1','region2','IF')
    }
    GM12878_GMSPIKEIN <- make_hicexp(data_list = sample_list,A.min=2,groups= c(1,1,1,1,2,2,2,2)) #A.min=5,
    ##Nomalization
    # Normalize 
    GM12878_GMSPIKEIN <- fastlo(GM12878_GMSPIKEIN)
    
    #GLM test
    d <- model.matrix(~factor(meta(GM12878_GMSPIKEIN)$group))
    GM12878_GMSPIKEIN_glm <- hic_glm(GM12878_GMSPIKEIN,design=d,coef=2)
    write_tsv(results(GM12878_GMSPIKEIN_glm),path =  paste0(out_path,"/binpairs_GMvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[k],"_Multi_glm"))

  }
}

options("scipen"=100)
TSetKey<-scan(paste0(in_path1,"True_signal5"),what = character())
#GM vs GM new simulation
in_path2<-paste0("/p/keles/fandingzhou/volumeA/DCI/filter/filter_as_",c("multiHiCcompare","diffHic","HiCDCPlus"))
kpset<-c()
trset<-c()
for(foldchange in c(2,3,4,6)){
  for(j in 1:4){
    data_j<-read_tsv(paste0(in_path1,"GM12878/seed1_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"))
    original_list<-cbind(rep(1,nrow(data_j)),data_j[,2]-res/2,data_j[,4]-res/2,data_j[,5])
    sample_list[[j+4]]<-original_list[abs(original_list[,2]-original_list[,3])>res*2,]
    colnames(sample_list[[j+4]])<-c('chr','region1','region2','IF')
  }
  HiCbinpairs<-multiHiCcompare_table(sample_list,nrep,1)
  binpairs<-paste0(HiCbinpairs$region1+res/2,"_",HiCbinpairs$region2+res/2)
  nosmooth<-intersect(TSetKey,binpairs)
  kpset<-c(kpset,length(binpairs))
  trset<-c(trset,length(nosmooth)/length(TSetKey))
  write(nosmooth, file=paste0(in_path2[1],"/Trueset/GMvsGMFC",foldchange,"_TSet"))
  write(binpairs, file=paste0(in_path2[1],"/Trueset/GMvsGMFC",foldchange,"_keep"))
  write_tsv(HiCbinpairs,path =  paste0(in_path2[1],"/multiHiCcompare/originaltable/binpairs_GMvsGMFC",foldchange,"_Multi_glm"))
}
write(kpset,file=paste0(in_path2[1],"/Trueset/num_kp"))
write(trset,file=paste0(in_path2[1],"/Trueset/rate_tko"))
      


#GM vs GM new simulation & filter as diffHic/HiCDCPlus
#for(k in 2:3){
options("scipen"=100)
TSetKey<-scan(paste0(in_path1,"True_signal5"),what = character())
#GM vs GM new simulation
in_path2<-paste0("/p/keles/fandingzhou/volumeA/DCI/filter/filter_as_",c("multiHiCcompare","diffHic","HiCDCPlus"))

k=2
  for(foldchange in c(2,3,4,6)){
    for(j in 1:4){
      data_j<-read_tsv(paste0(in_path1,"GM12878/strategy5_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"))
      original_list<-cbind(rep(1,nrow(data_j)),data_j[,2]-res/2,data_j[,4]-res/2,data_j[,5])
      sample_list[[j+4]]<-original_list[abs(original_list[,2]-original_list[,3])>res*2,]
      colnames(sample_list[[j+4]])<-c('chr','region1','region2','IF')
    }
    filter_pairs<-scan(paste0(in_path2[k],"/Trueset/GMvsGMFC",foldchange,"_keep"),what = character())
    HiCbinpairs<-multiHiCcompare_filter_table(sample_list,filter_pairs,bin_size=res,nrep)
    write_tsv(HiCbinpairs,path =  paste0(in_path2[k],"/multiHiCcompare/originaltable/binpairs_GMvsGMFC",foldchange,"_Multi_glm"))
  }
#}


multihic_result<-results(GM12878_GMSPIKEIN_glm)
multihic_result$p.adj=p.adjust(multihic_result$p.value,"BH")
#p.adjust(method = )
length(True_sig)
length(which(multihic_result$p.adj<0.1))
detected_sig<-multihic_result
detected_sig<-paste0(detected_sig$region1+bin_size/2,"_",detected_sig$region2+bin_size/2)
length(intersect(detected_sig,True_sig))
length(intersect(detected_sig[which(multihic_result$p.adj<0.1)],True_sig))
(54396-43263)/54396
head(which(multihic_result$p.adj<0.1&(!(detected_sig %in%True_sig))))
head(hic_table(GM12878_GMSPIKEIN))
a=which(multihic_result$p.adj<0.1&(!(detected_sig %in%True_sig)))
head(multihic_result[multihic_result$logFC<0,])

##GM vs GM
group_matrix<-matrix(c(1,1,2,2,1,2,1,2,1,2,2,1),4,3)
group_title<-c("12","13","14")
for(k in 2:3){
  GM12878 <- make_hicexp(data_list = sample_list,groups= group_matrix[,k],zero.p=0.4,A.min=5,filter=TRUE)
  GM12878 <- fastlo(GM12878)
  #GLM test
  d <- model.matrix(~factor(meta(GM12878)$group))
  GM12878_glm <- hic_glm(GM12878,design=d,coef=2)
  write_tsv(results(GM12878_glm),path =  paste0(out_path,"/binpairs_GMvsGM",group_title[k],"_Multi_glm"),col_names = TRUE)
}
# 
## GM vs A
for(j in 1:4){
  data_j<-read_tsv(paste0("/p/keles/scrna-seq/volumeB/freeHiC/A549/rep",j,"/chr1/s1_training/binPairs/rep",j,"_chr1.binPairs.chr1"),col_names=FALSE)
  original_list<-cbind(rep(1,nrow(data_j)),data_j[,2]-res/2,data_j[,4]-res/2,data_j[,5])
  sample_list[[j+4]]<-original_list[abs(original_list[,2]-original_list[,3])>res*2,]
  colnames(sample_list[[j+4]])<-c('chr','region1','region2','IF')
}
GM12878_A549<- make_hicexp(data_list = sample_list,groups= c(1,1,1,1,2,2,2,2),zero.p=0.4,filter=TRUE)

##Nomalization
# Normalize
GM12878_A549<- fastlo(GM12878_A549)

#GLM test
d <- model.matrix(~factor(meta(GM12878_A549)$group))
GM12878_A549_glm <- hic_glm(GM12878_A549,design=d,coef=2)
write_tsv(results(GM12878_A549_glm),path =  paste0(out_path,'/binpairs_GMvsA_Multi_glm'))

for(i in 1:8){
  a[[i]]<-paste0(sample_list[[i]]$region1,'_',sample_list[[i]]$region2)
}
b1<-union(a[[1]],a[[2]],a[[3]],a[[4]])
b2<-union(a[[5]],a[[6]],a[[7]],a[[8]])
b3<-intersect(b1,b2)
gm<-sample_list[[1]]%>%
  full_join(sample_list[[2]],by=c("chr","region1","region2"))%>%
  full_join(sample_list[[3]],by=c("chr","region1","region2"))%>%
  full_join(sample_list[[4]],by=c("chr","region1","region2"))
a<-sample_list[[5]]%>%
  full_join(sample_list[[6]],by=c("chr","region1","region2"))%>%
  full_join(sample_list[[7]],by=c("chr","region1","region2"))%>%
  full_join(sample_list[[8]],by=c("chr","region1","region2"))
a[is.na(a)]<-0
an<-a[,5]+a[,4]+a[,6]+a[,7]
a<-data.frame(a,an)
gm_a<-gm%>%
  full_join(a,by=c("chr","region1","region2"))
  
length(which(gm_a$gmn[is.na(gm_a$an)]>5)) #16166
length(which(gm_a$an[is.na(gm_a$gmn)]>5)) #71592->17287 if 5 is 10! but must be filtered out
head(gm_a[is.na(gm_a$gmn)&gm_a$an>5])

