#package directory
dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"

#function directory
func_path1="/p/keles/fandingzhou/volumeA/DCI/code/multiHiCcompare"
func_path2="/p/keles/fandingzhou/volumeA/DCI/code/diffHic"
func_path3="/p/keles/fandingzhou/volumeA/DCI/code/HiCDC+"

options("scipen"=100)

#packages for multiHiCcompare
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

#additional pacakges for diffHic
library(curl)
library(BiocParallel,lib.loc = dict_path)#
library(DelayedArray,lib.loc= dict_path)
library(SummarizedExperiment,lib.loc= dict_path)
library(InteractionSet,lib.loc= dict_path)
library(GenomeInfoDb,lib.loc= dict_path)
library(GenomicRanges,lib.loc= dict_path)
library(diffHic,lib.loc = dict_path)
library(csaw,lib.loc = dict_path)
library(statmod,lib.loc = dict_path)
library(ggplot2,lib.loc = dict_path)
library(parallel)

#additional pacakges for HiCDCPlus
library(progress,lib.loc=dict_path)
library(igraph,lib.loc=dict_path)
library(backports,lib.loc = dict_path)
library(checkmate,lib.loc = dict_path)
library(DESeq2,lib.loc = dict_path)# 
library(HiCDCPlus,lib.loc= dict_path)
library("BSgenome.Mmusculus.UCSC.mm9",verbose=T)


source(paste0(func_path1,"/multiHiCcompare_function.R"))
source(paste0(func_path2,"/diffHic_table.R"))


#real data

##multiHiCcompare
in_path="/p/keles/fandingzhou/volumeA/real_data"
chr=1
bin_size=20000
nrep=2
Amin=1
other_filtering=F
filter_file=NULL
inputfile=c(paste0(in_path,"/data/Auxin/chr",chr,"/rep",1:nrep,"_chr",chr,".binPairs.chr",chr),
             paste0(in_path,"/data/Untreated/chr",chr,"/rep",1:nrep,"_chr",chr,".binPairs.chr",chr))
reformat=T
long_range=T
TSetKey=NULL
ch_length=197195432
out_path1="/p/keles/fandingzhou/volumeA/real_data/filter_as_multiHiCcompare"

#inputfile: 4 or 8 hic tables (format the same with original GM12878)
#outpath: some outputs to outpath, including the table with pvalue and fvalue and kept binpairs after filtering
#res: bin_size
#nrep: number of replicates in one group
#Amin: filtering parameter - average IF for all replicates
#other_filtering: if true, use filtered binpairs from diffHic
#filter_file: if other_filtering is true, put in the filter_file (format: "binA_binB")
#reformat: if true, output file only contains 4 columns - bin_1 bin_2 F/Zstat Pvalue, otherwise, we can get the full table
#long_range: if true, filter out short-range interactions(<=2 bins)
#TSetkey: for simulated data, we can get the true sets after filtering

hic_table<-use_multiHiCcompare(inputfile,out_path1,res=bin_size,nrep,Amin=1,other_filtering=F,filter_file=NULL,reformat=T,long_range=T,TSetKey=NULL)
head(sort(p.adjust(hic_table$P)))
##diffHic
source(paste0(func_path2,"/contactmatrix.R"))
out_path2="/p/keles/fandingzhou/volumeA/real_data/filter_as_diffHic/diffHic"
out_path2_2="/p/keles/fandingzhou/volumeA/real_data/filter_as_diffHic"

#convert into symmetric matrices
cm<-multi_ContactMatrices(inputfile,out_path2,bin_size,ch_length,n.cores=4)
#use diffHic
hic_table2<-use_diffHic(cm,out_path2_2,bin_size,filter_par=1.5,reformat=T,other_filtering=F,filter_file=NULL,long_range=T,TSetKey=NULL,chr)
#filter as multiHiCcompare
filter_file_multi<-paste0(out_path1,"/Trueset/pairs_keep")
hic_table3<-use_diffHic(cm,out_path1,bin_size,filter_par=1.5,reformat=T,other_filtering=T,filter_file=filter_file_multi,long_range=T,TSetKey=NULL,chr)

#multiHiCcompare filtered as diffHic
filter_file_diffHic<-paste0(out_path2_2,"/Trueset/pairs_keep")
hic_table4<-use_multiHiCcompare(inputfile,out_path1,res=bin_size,nrep,Amin=1,other_filtering=T,filter_file=filter_file_diffHic,reformat=T,long_range=T,TSetKey=NULL)



##hic-dc+
source(paste0(func_path3,"/HiCDCPlus_function.R"))

out_path3="/p/keles/fandingzhou/volumeA/real_data/filter_as_HiCDCPlus/HiCDCPlus/hicmatrix"
if(!dir.exists(out_path3)){
  dir.create(out_path3)
}
n_bin<-nrow(cm[[1]])
##convert matrix to hic-pro format
a<-hicmatrix(inputfile,out_path3,bin_size=bin_size,ch1_length=ch_length,chr=chr)
b<-hicdcmatrix(out_path3,out_path3,nrep=2,chr=chr,bin_size=bin_size,ch1_length=ch_length,ncore=2,long_range=T)
out_path3.2="/p/keles/fandingzhou/volumeA/real_data/filter_as_diffHic/HiCDCPlus"
filter_file_multi<-paste0(out_path1,"/Trueset/pairs_keep")
filter_file_diffHic<-paste0(out_path2_2,"/Trueset/pairs_keep")

hicdcdiff2(out_path3,out_path3.2,filter_path,bin_size,chr)
  







#simulated data 4vs4 example
nrep<-4
rep<-c(2,3,4,6)
seed=11
foldchange=4
bin_size=40000
ch_length<-249250621

TSetKey<-scan(paste0("/p/keles/fandingzhou/volumeA/DCI/simulation/True_signal",seed),what = character())
inputfile=c(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep,"/chr1/s1_training/binPairs/rep",rep,"_chr1.binPairs.chr1"),
            paste0("/p/keles/fandingzhou/volumeA/DCI/simulation/GM12878/seed",seed,"_FC",foldchange,"_rep",rep,"_chr1.binPairs.chr1"))
out_path<-paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_multiHiCcompare/FC",foldchange)
#filter_pairs<-scan(paste0(in_path2[k],"/Trueset/GMvsGMFC",foldchange,"_keep"),what = character())
hic_table<-use_multiHiCcompare(inputfile,out_path,res=bin_size,nrep,Amin=1,other_filtering=F,filter_file=NULL,reformat=T,long_range=T,TSetKey)

out_path2<-paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_diffHic/FC",foldchange)
cm<-multi_ContactMatrices(inputfile,out_path2,bin_size,ch_length,n.core=4)
hic_table2<-use_diffHic(cm,out_path2,bin_size,filter_par=1,reformat=T,other_filtering=F,filter_file=NULL,long_range=T,TSetKey=NULL)





###2vs2



  