dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path3="/p/keles/fandingzhou/volumeA/ABcompartment/A549"
in_path4="/p/keles/fandingzhou/volumeA/ABcompartment/GM12878"
in_path5="/p/keles/fandingzhou/volumeA/ABcompartment/GM12878_spikein"
in_path6="/p/keles/fandingzhou/volumeA/ABcompartment/pca"
out_path=c("/p/keles/fandingzhou/volumeA/DCI/ABcompartment")
func_path="/p/keles/fandingzhou/volumeA/DCI/code/weight/"
library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
bin_size<-40000
ch_length<-249250621
nbin=floor(ch_length/bin_size)
source(paste0(func_path,"boundary_function.R"))
#boundary of ABcompartment for A549 and GM12878
##large compartments A549
data_origin1<-read_tsv(paste0(in_path3,"/chr1_cscore.txt"),col_names=FALSE)
#GM12878
data_origin2<-read_tsv(paste0(in_path4,"/chr1_cscore.txt"),col_names=FALSE)

#GM12878_spikein
data_origin3<-read_tsv(paste0(in_path5,"/chr1_cscore.txt"),col_names=FALSE)

#GM12878_foldchange
for(seed in 11:20){
  #length(ABcompartment_boundaries(data_origin2,bin_size,nbin))
  data_origin4<-read.table(paste0(in_path6,"/seed",seed,"_GM_FC",foldchange,"_eigen_no0"),header=T)
  boundary_GM_FC<-conbined_boundaries(data_origin4,data_origin2,bin_size,nbin)
  print(length(boundary_GM_FC))
  write.table(boundary_GM_FC,file=paste0(out_path,"/seed",seed,"GM_FC_boundaries"))
}

##combining boundaries
boundary_GM_A<-conbined_boundaries(data_origin1,data_origin2,bin_size,nbin)
boundary_GM_spikein<-conbined_boundaries(data_origin3,data_origin2,bin_size,nbin)

boundary_GM_GM<-ABcompartment_boundaries(data_origin2,bin_size,nbin)
write.table(boundary_GM_A,file=paste0(out_path,"/GM_A_boundaries"))
write.table(boundary_GM_A,file=paste0(out_path,"/GM_spikein_boundaries"))
write.table(boundary_GM_GM,file=paste0(out_path,"/GM_GM_boundaries"))

