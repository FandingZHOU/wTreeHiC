dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"

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

#change a function to get f-value
.glm_reformat <- function(result, hic_table, p.method) {
  # create table of location info and p-value results
  result <- cbind(hic_table, result$table)
  colnames(result)[ncol(result)] <-"p.value"
  # adjust p-values
  result$p.adj <- p.adjust(result$p.value, method = p.method)
  # convert logFC from natural log to log2
  result$logFC <- log2(exp(result$logFC))
  return(result)
}

out_path="/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/originaltable"
seqdepth<-c(1,3,5)
rep<-c(2,3,4,6)
res <- 40000 
scalesize<-c(2,3,5,10)
# Read data

sample_list <- list()
for(j in 1:4){
  data_j<-read_tsv(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j],"/chr1/s1_training/binPairs/rep",rep[j],"_chr1.binPairs.chr1"),col_names=FALSE)
  original_list<-cbind(rep(1,nrow(data_j)),data_j[,2]-res/2,data_j[,4]-res/2,data_j[,5])
  sample_list[[j]]<-original_list[abs(original_list[,2]-original_list[,3])>res*2,]
  colnames(sample_list[[j]])<-c('chr','region1','region2','IF')
}

# GM vs GM spike-in
for(k in 1:3){
  for(j in 1:4){
    data_j<-read_tsv(paste0("/p/keles/treehic/volumeD/freeHiC/GM12878/rep",rep[j],"/chr1/spikeIn3_40kb_smooth_seqx",seqdepth[k],"/simuProcess/s3_binPairs/rep",rep[j],"_chr1.binPairs"),col_names=FALSE)%>% filter(X1 == "chr1", X3 == "chr1")
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
  write_tsv(results(GM12878_GMSPIKEIN_glm),path =  paste0(out_path,"/binpairs_GMvsGMspikein_seqx",seqdepth[k],"_Multi_glm"))
  print(nrow(results(GM12878_GMSPIKEIN_glm)))
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
    GM12878_GMSPIKEIN <- make_hicexp(data_list = sample_list,A.min=2,groups= c(1,1,1,1,2,2,2,2)) #,zero.p=0.4,
    ##Nomalization
    # Normalize 
    GM12878_GMSPIKEIN <- fastlo(GM12878_GMSPIKEIN)
    
    #GLM test
    d <- model.matrix(~factor(meta(GM12878_GMSPIKEIN)$group))
    GM12878_GMSPIKEIN_glm <- hic_glm(GM12878_GMSPIKEIN,design=d,coef=2)
    write_tsv(results(GM12878_GMSPIKEIN_glm),path =  paste0(out_path,"/binpairs_GMvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[k],"_Multi_glm"))
    print(nrow(results(GM12878_GMSPIKEIN_glm)))
  }
}


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

