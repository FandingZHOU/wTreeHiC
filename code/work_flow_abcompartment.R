##Calculating O/E Matrix
options(echo=TRUE)
options("scipen"=100)

args <- commandArgs(trailingOnly = TRUE)
print(args)


#input matrix
combined_matrix_dir = args[1]

#bin_size
bin_size = args[2]

#output directory
out_path = args[3]

#packages dictionary
dict_path = args[4]

#Chr length
ch1_length=args[5]

#n cores
n.cores=arg[6]

#Function path
func_path=arg[7]

###


library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(parallel)


n_bin<-floor(ch1_length/bin_size)
delta <- rep(seq_len(n_bin), n_bin) - 
  rep(seq_len(n_bin), each = n_bin)

distance<-0:(n_bin-1)

obs_matrix<-as.matrix(read_tsv(combined_matrix_dir,col_names=F))
expect_no0<-get_expected_matrix(obs_matrix,n_bin,bin_size,n.cores)
write.table(expect_no0,paste0(out_path,"/expect_no0.tsv"))
EvsO_no0_matrix<-get_OvsE_matrix(obs_matrix,expect_no0$expected,n_bin)
write.table(EvsO_no0_matrix,paste0(out_path,"/eovero_no0.tsv"))
compartment_boundaries<-get_ABcompartment(EvsO_no0_matrix)
write_tsv(compartment_boundaries,paste0(out_path,"/eigen_no0"),col_names= F)
