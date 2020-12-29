dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path2="/p/keles/fandingzhou/volumeA/DCI/diffHic/data/filter"
options("scipen"=100)
library(HiCDCPlus)
library("BSgenome.Hsapiens.UCSC.hg19",lib.loc = dict_path)


ch1_length<-249250621
bin_size<-40000
rep<-c(2,3,4,6)
seqdepth<-c(1,3,5)
nrep<-4
scalesize<-c(2,3,5,10)

gi_list<-generate.binned.gi.list(bin_size,Dthreshold = 40000,chrs=paste0('chr',1))

n_bin<-floor(ch1_length/bin_size)
all.regions<-data.frame(cbind(rep("chr1",n_bin+1),seq(1,n_bin*bin_size+1,bin_size),c(seq(bin_size,n_bin*bin_size,bin_size),ch1_length),1:(n_bin+1)))
data.table::fwrite(all.regions,file=paste0(in_path2,"/bedfile"),sep = "\t",col.names = F)
head(read.table(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j],"/chr1/s1_training/validPairs/rep",rep[j],"_chr1.validPairs.chr1")))
a<-read.table(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j],"/chr1/s1_training/binPairs/rep",rep[j],"_chr1.binPairs.chr1"))[,c(2,4,5)]
data.table::fwrite(a,paste0(in_path2,"/matrix"),sep = "\t",col.names = F)
gi_list<-add.hicpro.allvalidpairs.counts(
  gi_list,
  paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j],"/chr1/s1_training/validPairs/rep",rep[j],"_chr1.validPairs.chr1"),
  chrs = NULL,
  binned = TRUE,
  add_inter = FALSE
)
gi_list<-add.hicpro.matrix.counts(
  gi_list,
  paste0(in_path2,"/bedfile"),
  paste0(in_path2,"/matrix"),
)
data_j<-read_tsv