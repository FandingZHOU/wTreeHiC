dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path="/p/keles/fandingzhou/volumeA/DCI/HiCDCPlus/data/hicpromatrix"
in_path2="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/filter"
out_path="/p/keles/fandingzhou/volumeA/DCI/HiCDCPlus/data/original/"
out_path2="/p/keles/fandingzhou/volumeA/DCI/HiCDCPlus/data/filter/"
options("scipen"=100)
#install.packages('xtable',lib = dict_path)#htmlwidgets,stringr
library(progress,lib.loc=dict_path)
library(igraph,lib.loc=dict_path)
library(backports,lib.loc = dict_path)
library(checkmate,lib.loc = dict_path)
library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(DESeq2,lib.loc = dict_path)# 
library(HiCDCPlus,lib.loc= dict_path)
library("BSgenome.Hsapiens.UCSC.hg19",lib.loc = dict_path)
library("BSgenome.Hsapiens.UCSC.hg38",lib.loc = dict_path)
library(dplyr,lib.loc = dict_path)# 

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

ch1_length<-249250621
bin_size<-40000
rep<-c(2,3,4,6)
seqdepth<-c(1,3,5)
nrep<-4
scalesize<-c(2,3,5,10)
n_bin<-floor(ch1_length/bin_size)



indexfile<-data.frame()
for(j in 1:4){
  gi_list<-generate_binned_gi_list(bin_size,chrs=paste0('chr',1),Dthreshold = ch1_length)
  
  gi_list<-add_hicpro_matrix_counts(
    gi_list,
    paste0(in_path,"/bedfile"),
    paste0(in_path,"/GM_original_rep",rep[j])
  )
  
  gi_list<-expand_1D_features(gi_list)
  #run HiC-DC+ on 2 cores
  set.seed(1010) #HiC-DC downsamples rows for modeling
  gi_list<-HiCDCPlus_parallel(gi_list,Dmin = 3*bin_size,Dmax = ch1_length,ncore=2)
  gi_list_qval<-gi_list[[1]]$qvalue
  gi_list[[1]]<-gi_list[[1]][(!is.na(gi_list_qval))]
  indexfile<-unique(rbind(indexfile,
                          as.data.frame(gi_list[[1]][gi_list[[1]]$qvalue<=0.05])[c('seqnames1','start1','start2')]))
  #gi_list_write(gi_list,fname=paste0(in_path,"/GM_original_rep",rep[j],"sig_interact.txt.gz"))
  saveRDS(gi_list,paste0(in_path,"/GM_original_rep",rep[j],"sig_interact.rds"))
}



##new sim

indexfile<-data.table::fread(paste0(in_path,'/GM_origin_filter.txt'))

for(seed in 14:20){
  for(foldchange in c(2,3,4,6)){
    if(seed !=14|!foldchange %in% c(2,3,4)){
    indexfile_sim<-indexfile
    for(j in 1:4){
      gi_list<-generate_binned_gi_list(bin_size,chrs=paste0('chr',1),Dthreshold = ch1_length)
      
      gi_list<-add_hicpro_matrix_counts(
        gi_list,
        paste0(in_path,"/bedfile"),
        paste0(in_path,"/simulation/seed",seed,"_FC",foldchange,"_GM",rep[j])
      )
      
      gi_list<-expand_1D_features(gi_list)
      #run HiC-DC+ on 2 cores
      set.seed(1010) #HiC-DC downsamples rows for modeling
      gi_list<-HiCDCPlus_parallel(gi_list,Dmin = 3*bin_size,Dmax = ch1_length,ncore=10)
      gi_list_qval<-gi_list[[1]]$qvalue
      gi_list[[1]]<-gi_list[[1]][(!is.na(gi_list_qval))]
      indexfile_sim<-unique(rbind(indexfile_sim,
                                  as.data.frame(gi_list[[1]][gi_list[[1]]$qvalue<=0.05])[c('seqnames1',
                                                                                           'start1','start2')]))
      #gi_list_write(gi_list,fname=paste0(in_path,"/GM_spikein_rep",rep[j],"_seqx",seqdepth[k],"sig_interact.txt.gz"))
      saveRDS(gi_list,paste0(in_path,"/simulation/seed",seed,"_FC",foldchange,"_GM",rep[j],"sig_interact.rds"))
    }
    colnames(indexfile_sim)<-c('chr','startI','startJ')
    data.table::fwrite(indexfile_sim,
                       paste0(in_path,"/simulation/seed",seed,"_FC",foldchange,"_GM",rep[j],'_indices.txt.gz'),
                       sep='\t',row.names=FALSE,quote=FALSE)
  }
  }
}


for(foldchange in c(2,3,4,6)){
  out_path1<-paste0(out_path2,"simulation/GM_filter_FC",foldchange,"/")
  hicdcdiff(input_paths=list(origin1=paste0(in_path,"/GM_original_rep",rep,"sig_interact.rds"),
                             spikein=paste0(in_path,"/simulation/seed1_FC",foldchange,"_GM",rep,"sig_interact.rds")),
            output_path=paste0(out_path1),
            filter_file=paste0(out_path2,"simulation/GM_filter_FC",foldchange,'_indices.txt.gz'),
            Dmin = 3*bin_size,
            Dmax = ch1_length,
            binsize=40000,
            DESeq.save=TRUE,
            diagnostics=T
  )
}

#filter by hicdcplus
in_path3=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/filter_as_",c("multiHiCcompare","diffHic","HiCDCPlus"))

for(foldchange in c(2,3,4,6)){
  out_path1<-paste0(in_path3[3],"/HiCDCPlus/GM_filter_FC",foldchange,"/")
  hicdcdiff(input_paths=list(origin=paste0(in_path,"/GM_original_rep",rep[1:2],"sig_interact.rds"),
                             spikein=paste0(in_path,"/simulation/sim3_FC",foldchange,"_GM",rep[1:2],"sig_interact.rds")),
            output_path=paste0(out_path1),
            filter_file=paste0(in_path,"/simulation/seed1_FC",foldchange,'_indices.txt.gz'),
            Dmin = 3*bin_size,
            Dmax = ch1_length,
            binsize=40000,
            granularity=40000,
            DESeq.save=TRUE,
            diagnostics=F
  )
}

#filter by diffHic
for(seed in 12:20){
#seed=1
  in_path4=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/2vs2/filter_as_",c("multiHiCcompare","diffHic","HiCDCPlus"))
  in_path3=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_",c("multiHiCcompare","diffHic","HiCDCPlus"))
  for(hh in 1:2){
    if(seed!=12|hh!=1){
    for(foldchange in c(2,3,4,6)){
      #
        filter1<-scan(paste0(in_path3[hh],"/Trueset/GMvsGMFC",foldchange,"_keep"),what = character())
        filterset<-strsplit(x=filter1, split="_")
        indexfile<- data.frame()
        indexfile<-cbind(rep('chr1',length(filter1)),sapply(filterset, binA,bin_size),sapply(filterset, binB,bin_size))
        colnames(indexfile)<-c('chr','startI','startJ')
        data.table::fwrite(indexfile,
                           paste0(in_path3[hh],"/HiCDCPlus/GM_filter_FC",foldchange,'_indices.txt.gz'),
                           sep='\t',row.names=FALSE,quote=FALSE)
        hicdcdiff(input_paths=list(origin=paste0(in_path,"/GM_original_rep",rep[2:3],"sig_interact.rds"),
                                   spikein=paste0(in_path,"/simulation/seed",seed,"_FC",foldchange,"_GM",rep[1:2],"sig_interact.rds")),
                  output_path=paste0(in_path4[hh],"/HiCDCPlus/originaltable/GM_FC",foldchange,"/"),
                  filter_file=paste0(in_path3[hh],"/HiCDCPlus/GM_filter_FC",foldchange,'_indices.txt.gz'),
                  Dmin = 3*bin_size,
                  Dmax = ch1_length,
                  binsize=40000,
                  granularity=40000,
                  DESeq.save=TRUE,
                  diagnostics=F
        )
        hicdcdiff(input_paths=list(origin=paste0(in_path,"/GM_original_rep",rep,"sig_interact.rds"),
                                   spikein=paste0(in_path,"/simulation/seed",seed,"_FC",foldchange,"_GM",rep,"sig_interact.rds")),
                  output_path=paste0(in_path3[hh],"/HiCDCPlus/originaltable/GM_FC",foldchange,"/"),
                  filter_file=paste0(in_path3[hh],"/HiCDCPlus/GM_filter_FC",foldchange,'_indices.txt.gz'),
                  Dmin = 3*bin_size,
                  Dmax = ch1_length,
                  binsize=40000,
                  granularity=40000,
                  DESeq.save=TRUE,
                  diagnostics=F
        )
      }
    }
     }
      
  
  
}



###original spike-in
data.table::fwrite(indexfile,
                   paste0(in_path,'/GM_origin_filter.txt'),
                   sep='\t',row.names=FALSE,quote=FALSE)

indexfile<-data.table::fread(paste0(in_path,'/GM_origin_filter.txt'))
for(k in 1:3){
  indexfile_spikein<-indexfile
  for(j in 1:4){
    gi_list<-generate_binned_gi_list(bin_size,chrs=paste0('chr',1),Dthreshold = ch1_length)
    
    gi_list<-add_hicpro_matrix_counts(
      gi_list,
      paste0(in_path,"/bedfile"),
      paste0(in_path,"/GM_spikein_rep",rep[j],"_seqx",seqdepth[k])
    )
    
    gi_list<-expand_1D_features(gi_list)
    #run HiC-DC+ on 2 cores
    set.seed(1010) #HiC-DC downsamples rows for modeling
    gi_list<-HiCDCPlus_parallel(gi_list,Dmin = 3*bin_size,Dmax = ch1_length,ncore=2)
    gi_list_qval<-gi_list[[1]]$qvalue
    gi_list[[1]]<-gi_list[[1]][(!is.na(gi_list_qval))]
    indexfile_spikein<-unique(rbind(indexfile_spikein,
                                    as.data.frame(gi_list[[1]][gi_list[[1]]$qvalue<=0.05])[c('seqnames1',
                                                                                             'start1','start2')]))
    #gi_list_write(gi_list,fname=paste0(in_path,"/GM_spikein_rep",rep[j],"_seqx",seqdepth[k],"sig_interact.txt.gz"))
    saveRDS(gi_list,paste0(in_path,"/GM_spikein_rep",rep[j],"_seqx",seqdepth[k],"sig_interact.rds"))
  }
  colnames(indexfile_spikein)<-c('chr','startI','startJ')
  data.table::fwrite(indexfile_spikein,
                     paste0(in_path,"/GM_spikein_filter_seqx",seqdepth[k],'_indices.txt.gz'),
                     sep='\t',row.names=FALSE,quote=FALSE)
}


#check filering
for(k in 1:3){
  filter1<-scan(paste0(in_path2,"/_GMvsvsGMspikein_seqx",seqdepth[k],"_keep"),what = character())
  filter2<-data.table::fread(paste0(in_path,"/GM_spikein_filter_seqx",seqdepth[k],'_indices.txt.gz'))
  filter3<-paste0(filter2[,2]+bin_size/2,"_",filter2[,3]+bin_size/2)
  length(interaction(filter1,filter3))
  length(filter1)
  length(filter3)
}

#Check result
#a<-readRDS("/p/keles/fandingzhou/volumeA/DCI/HiCDCPlus/data/originalchr1_DESeq2_obj.rds")
TSetKey<-scan(paste0(in_path2,"/TSetKey"),what = character())

for(k in 1:3){
  result1<-results(readRDS(paste0(out_path2,"GM_spikein_seqx",seqdepth[k],"/chr1_DESeq2_obj.rds")))
  tset<-strsplit(x=rownames(result1), split=":")
  
  result1$Bin1<-sapply(tset, binA,bin_size,T)
  result1$Bin2<-sapply(tset, binB,bin_size,T)
  
  padj<-p.adjust(result1$pvalue,"BH")
  #padj2<-na.omit(result1$padj)
  
  #tru<-paste0(result1$Bin1[result1$pvalue<0.1],"_",result1$Bin2[result1$pvalue<0.1])
  
  tru<-paste0(result1$Bin1[padj<0.1],"_",result1$Bin2[padj<0.1])
  print(length(intersect(tru,TSetKey)))
  print(length(intersect(tru,TSetKey))/length(tru))
}

# [1] 727
# [1] 0.9405
# [1] 594
# [1] 0.8959
# [1] 372
# [1] 0.767



#filter by intersection of multi/diff/hicdc+
for(k in 1:3){
  filter1<-scan(paste0(in_path2,"/_GMvsvsGMspikein_seqx",seqdepth[k],"_keep"),what = character())
  filterset<-strsplit(x=filter1, split="_")
  indexfile<- data.frame()
  indexfile<-cbind(rep('chr1',length(filter1)),sapply(filterset, binA,bin_size),sapply(filterset, binB,bin_size))
  colnames(indexfile)<-c('chr','startI','startJ')
  data.table::fwrite(indexfile,
                     paste0(out_path2,"GM_spikein_filter_seqx",seqdepth[k],'_indices.txt.gz'),
                     sep='\t',row.names=FALSE,quote=FALSE)
}

for(k in 1:3){
  out_path1<-paste0(out_path2,"GM_spikein_seqx",seqdepth[k],"/")
  hicdcdiff(input_paths=list(origin1=paste0(in_path,"/GM_original_rep",rep,"sig_interact.rds"),
                             spikein=paste0(in_path,"/GM_spikein_rep",rep,"_seqx",seqdepth[k],"sig_interact.rds")),
            output_path=paste0(out_path1),
            filter_file=paste0(out_path2,"/GM_spikein_filter_seqx",seqdepth[k],'_indices.txt.gz'),
            Dmin = 3*bin_size,
            Dmax = ch1_length,
            binsize=40000,
            DESeq.save=TRUE,
            diagnostics=T
  )
}
HiCbinpairs<-results(readRDS(paste0("/p/keles/fandingzhou/volumeA/DCI/HiCDCPlus/data/originalchr1_DESeq2_obj.rds")))

k=3
a<-data.table::fread(
  paste0(out_path2,"GM_spikein_filter_seqx",seqdepth[k],'_indices.txt.gz'))

#filter by hicdcplus
for(k in 1:3){
  out_path1<-paste0(out_path,"GM_spikein_seqx",seqdepth[k],"/")
  hicdcdiff(input_paths=list(origin=paste0(in_path,"/GM_original_rep",rep,"sig_interact.rds"),
                             spikein=paste0(in_path,"/GM_spikein_rep",rep,"_seqx",seqdepth[k],"sig_interact.rds")),
            output_path=paste0(out_path1),
            filter_file=paste0(in_path,"/GM_spikein_filter_seqx",seqdepth[k],'_indices.txt.gz'),
            Dmin = 3*bin_size,
            Dmax = ch1_length,
            binsize=40000,
            fitType='parametric',
            DESeq.save=TRUE,
            diagnostics=T
  )
}



res<-readRDS(paste0(out_path2,"GM_spikein_seqx",seqdepth[k],"/chr1_DESeq2_obj.rds"))
HiCbinpairs<-results(res, cooksCutoff=F,independentFiltering=FALSE)

HiCbinpairs<-results(res, cooksCutoff=F,independentFiltering=F)

HiCbinpairs<-results(res)


for(hh in 1:2){
  hh=2
  for(foldchange in c(2,3,4,6)){
    #if(hh!=1 | !(foldchange %in% c(2,3,4))){
      
 # }
}

input_paths=list(origin=paste0(in_path,"/GM_original_rep",rep,"sig_interact.rds"),
                 spikein=paste0(in_path,"/simulation/sim3_FC",foldchange,"_GM",rep,"sig_interact.rds"))
output_path=paste0(in_path3[1],"/HiCDCPlus/originaltable/GM_FC",foldchange,"/")
filter_file=paste0(in_path3[1],"/HiCDCPlus/GM_filter_FC",foldchange,'_indices.txt.gz')
Dmin = 3*bin_size
Dmax = ch1_length
binsize=40000
DESeq.save=TRUE
diagnostics=T
bin_type = "Bins-uniform"
binsize = 40000
granularity = 40000
diagnostics = FALSE
DESeq.save = FALSE
fitType = "local"
chr="chr1"

filter_file=paste0(in_path,"/simulation/sim3_FC",foldchange,"_GM",rep[4],'_indices.txt.gz')




res<-readRDS(paste0(out_path2,"simulation/GM_filter_FC",foldchange,"/chr1_DESeq2_obj.rds"))
HiCbinpairs<-results(res, cooksCutoff=F,independentFiltering=FALSE)

HiCbinpairs<-results(res, cooksCutoff=F,independentFiltering=F)

HiCbinpairs<-results(res)
nrow(HiCbinpairs[HiCbinpairs$padj<0.1,])
