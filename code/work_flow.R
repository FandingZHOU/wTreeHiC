options(echo=TRUE)
options("scipen"=100)
args <- commandArgs(trailingOnly = TRUE)
print(args)

#Chr length
ch1_length=args[7]

#weight method
weight_method = args[1]

#bin_size
bin_size = args[2]

#output directory
outputdir = args[3]

#input file name
inputfile = args[4]

#methods
meth = args[5]

#packages dictionary
dict_path = args[6]

#reformat
reformat=args[8]

#A/B compartment result
compartment_input1=arg[9]
compartment_input2=arg[10]

#TAD result
tadfile=arg[12]
#threshold (can be vector or a singal value)
threshold=arg[13]

# Function path
func_path=arg[14]

#n clusters
clustern=arg[11]
###


#load packages
library(utf8,lib.loc = dict_path)
library(cli,lib.loc = dict_path)
library(fansi,lib.loc = dict_path)
library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(DESeq2,lib.loc = dict_path)# 
library(splines)
library(qvalue,lib.loc = dict_path)

#load functions
source(paste0(func_path,"boundary_function.R"))
source(paste0(func_path,"group_division_function.R"))
source(paste0(func_path,"reformat_function.R"))
source(paste0(func_path,"weight_function.R"))
source(paste0(func_path2,"filer_group_function"))
source(paste0(func_path2,"get_fdr_tp_tree_function.R"))
source(paste0(func_path2,"get_fdr_tp_function.R"))
source(paste0(func_path2,"tree_layer_function.R"))




#create dictionary
if(!dir.exists(outputdir)){
  dir.create(outputdir)
}

#reformat
if(reformat){
  HiCbinpairs_data<-glm.reformat(inputfile,meth)
}else{
  HiCbinpairs_data<-read_tsv(inputfile,col_names=T)
}

##combining boundaries
nbin=floor(ch_length/bin_size)
#GM12878
compartment1<-read_tsv(compartment_input1,col_names=FALSE)
#GM12878 spike-in
compartment2<-read_tsv(compartment_input2,col_names=FALSE)
boundary1<-conbined_boundaries(data_origin4,data_origin2,bin_size,nbin)
write.table(boundary1,file=paste0(outputdir,"/group_boundaries"))

## group partition
boundary1<- boundary1[,1]
boundary1<-c(boundary1,ch1_length)
HiCbinpairs_data<-group_division(HiCbinpairs_data,boundary1,clustern)
write.table(HiCbinpairs_data, file=paste0(outputdir,"/grouptable"), sep="\t",quote=FALSE, row.names=FALSE)

## weight function
if(meth!="HiCDCPlus"){
  statistic="F"
}else{
  statistic="Z"
}
HiCbinpairs_data<-add_weight(HiCbinpairs_data,weight_method,statistic)
write.table(HiCbinpairs_data, file=paste0(outputdir,"/weighttable"), sep="\t",quote=FALSE, row.names=FALSE)



## Tree
TADs<-read_tsv(tadfile,col_names=F)

HiCbinpairs_data<-na.omit(HiCbinpairs_data)
HiCbinpairs_data<-tree_layer(HiCbinpairs_data,TADs,bin_size)
max_layer<-max(TADs$X3)-3
tree_dci<-get_dci_tree(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,weight_method,gamma=0.05)
original_dci<-get_dci(HiCbinpairs_data,add_weight=FALSE,threshold)
weight_dci<-get_dci(HiCbinpairs_data,add_weight=T,threshold)
wtree_dci<-get_dci_tree(HiCbinpairs_data,add_weight=T,threshold)
saveRDS(tree_dci,paste0(outputdir,"/tree_dci"))
saveRDS(original_dci,paste0(outputdir,"/original_dci"))
saveRDS(weight_dci,paste0(outputdir,"/weight_dci"))
saveRDS(wtree_dci,paste0(outputdir,"/wtree_dci"))
# write(tree_dci,paste0(outputdir,"/tree_dci"))
# write(original_dci,paste0(outputdir,"/original_dci"))
# write(weight_dci,paste0(outputdir,"/weight_dci"))
# write(wtree_dci,paste0(outputdir,"/wtree_dci"))
fp_tree<-obtain_power_fdr(tree_dci,TSetKey)
fp_weight<-obtain_power_fdr(tree_dci,TSetKey)
fp_wtree<-obtain_power_fdr(tree_dci,TSetKey)
fp_original<-obtain_power_fdr(tree_dci,TSetKey)

