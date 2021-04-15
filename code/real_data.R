
dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
func_path="/p/keles/fandingzhou/volumeA/DCI/code/multiHiCcompare"
in_path="/p/keles/fandingzhou/volumeA/real_data"
#out_path="/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/originaltable"

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

dat <- cooler2bedpe(path = paste0(in_path,"/GSM2644946_Untreated-R2.20000.cool"))
                    
for(i in 1:2){
  dat <- cooler2bedpe(path = paste0(in_path,"/GSM",2644946+i,"_Auxin2days-R",i,".20000.cool")
  )
  #dat <- cooler2bedpe(path = paste0(in_path,"/GSM2644946_Untreated-R2.20000.cool"))
  for(chr in c(1:19,"X","Y")){
    out_path=paste0(in_path,"/data/Auxin/chr",chr)
    #out_path=paste0(in_path,"/data/Untreated/chr",chr)
    if(!file.exists(out_path)){
      dir.create(out_path)
    }
    HiCbinpairs<-dat$cis[[paste0("chr",chr)]]
    HiCbinpairs<-HiCbinpairs %>%
      mutate(
        binA=(start1+end1)/2,
        binB=(start2+end2)/2
      )
    HiCbinpairs_data<-HiCbinpairs[,c(1,8,4,9,7)]
    write_tsv(HiCbinpairs_data,paste0(out_path,"/rep",i,"_chr",chr,".binPairs.chr",chr))
  }
}



