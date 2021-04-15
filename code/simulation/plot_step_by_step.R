dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path1="/p/keles/fandingzhou/volumeA/DCI/simulation/"
options("scipen"=100)

library(utf8,lib.loc = dict_path )
library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 
library(parallel)

ch1_length<-249250621
bin_size<-40000
rep<-c(2,3,4,6)
nrep<-4
meth=c("multiHiCcompare","diffHic","HiCDCPlus")
#prob_FC<-0.1
prob_FC<-0.05

#Truesig=scan(paste0(out_path,"True_signal2"),what = character())
Truesig=scan(paste0(in_path1,"True_signal12"),what = character())
in_path1="/p/keles/fandingzhou/volumeA/DCI/simulation/"
seed=12


for(foldchange in c(2,3,4,6)){
  
  foldchange=4
data_j<-list()
#nrep=2
  for(j in 1:nrep){
     data_j[[j]]<-read_tsv(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j],"/chr1/s1_training/binPairs/rep",rep[j],"_chr1.binPairs.chr1"),col_names=FALSE)%>% 
       filter(X1 == "chr1", X3 == "chr1")
    #data_j[[j]]<-read_tsv(paste0("/p/keles/scrna-seq/volumeC/freeHiC/GM12878/rep",rep[j+1],"/chr1/s1_training/binPairs/rep",rep[j+1],"_chr1.binPairs.chr1"),col_names=FALSE)%>%
     # filter(X1 == "chr1", X3 == "chr1")

    #data_j[[2]]<-read_tsv(paste0(out_path,"GM12878/less_truesig_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"))
    data_j[[j+nrep]]<-read_tsv(paste0(in_path1,"GM12878/seed",seed,"_FC",foldchange,"_rep",rep[j],"_chr1.binPairs.chr1"))
    colnames(data_j[[j+nrep]])<-paste0("X",1:5)
  }
  merged_matrix<-Reduce(function(x,y) full_join(x,y,by = c("X1","X3","X2","X4")),data_j);
  merged_matrix<-merged_matrix %>%
      mutate_all(~replace(.,is.na(.),0))
  colnames(merged_matrix)<-c("chrA","region1","chrB","region2",paste0("IF",1:2*nrep))
  allpairs=paste0(merged_matrix$region1,"_",merged_matrix$region2)
  trueindex = which(allpairs %in% Truesig)
  falseindex=which(!allpairs %in% Truesig)
  set.seed(1)
  samptrue = sample(trueindex,length(trueindex)/6)
  sampfalse = sample(falseindex,length(falseindex)/6)
  before = log2(rowMeans(merged_matrix[,5:(4+nrep)])+1)
  after = log2(rowMeans(merged_matrix[,(5+nrep):(4+2*nrep)])+1)
  plot((before[samptrue]+after[samptrue])/2,after[samptrue]-before[samptrue],col="red",main=paste0("seed",seed,"_FC",foldchange),pch = 46,xlab="log2(avg(original count)+1)/2+log2(avg(FC count)+1)/2",ylab="log2(avg(FC count)+1)-log2(avg(original count)+1)",ylim = c(-5,5),xlim=c(0,9))
  points((before[sampfalse]+after[sampfalse])/2,after[sampfalse]-before[sampfalse],pch = 46)
    #plot((before[!trueindex]+after[!trueindex])/2,after[!trueindex]-before[!trueindex],main=paste0("rep",rep[j],"_FC",foldchange),pch = 46)
    #points((before[trueindex]+after[trueindex])/2,after[trueindex]-before[trueindex],col="red",pch = 46)
  #}
}
#q()
  #plot filter
  n=2
  in_path2=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_",meth[n])
  in_path3=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/2vs2/filter_as_",meth[n])
  
  kept_pairs<-scan(paste0(in_path2,"/Trueset/GMvsGMFC",foldchange,"_keep"),what = character())
  kept_pairs<-which(allpairs %in% kept_pairs)
  filter_pairs<-which(!allpairs %in% kept_pairs)
  filter_pairs<-intersect(c(samptrue,sampfalse),filter_pairs)
  kept_pairs_true<-intersect(c(samptrue),kept_pairs)
  kept_pairs_false<-intersect(c(sampfalse),kept_pairs)
  plot((before[kept_pairs_true]+after[kept_pairs_true])/2,after[kept_pairs_true]-before[kept_pairs_true],col="red",main=paste0("seed",seed,"_FC",foldchange,"_keptsig"),pch = 46,xlab="log2(avg(original count)+1)/2+log2(avg(FC count)+1)/2",ylab="log2(avg(FC count)+1)-log2(avg(original count)+1)",ylim = c(-5,5),xlim=c(0,9))
  points((before[kept_pairs_false]+after[kept_pairs_false])/2,after[kept_pairs_false]-before[kept_pairs_false],pch = 46)
  #plot((before[filter_pairs]+after[filter_pairs])/2,after[filter_pairs]-before[filter_pairs],col="darkgrey",pch = 46,main=paste0("seed",seed,"_FC",foldchange,"filteredsig"),xlab="log2(avg(original count)+1)/2+log2(avg(FC count)+1)/2",ylab="log2(avg(FC count)+1)-log2(avg(original count)+1)",ylim = c(-5,5),xlim=c(0,9))
  
  
  Truesig2=scan(paste0(in_path2,"/Trueset/GMvsGMFC",foldchange,"_TSet"),what = character())
  TSetKey=data.frame(Truesig2,Truesig2)
  for(m in c(2:3,1)){
  #m=3
    #plot original
    in_path4=paste0(in_path2,"/",meth[m])
    #in_path4=paste0(in_path3,"/",meth[m])
    HiCbinpairs_data<-read_tsv(paste0(in_path4,"/weighttable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_weightedp_xi.tsv"))
    HiCbinpairs_data<-na.omit(HiCbinpairs_data)
    detected_orig_p<-get_fdr_dci(HiCbinpairs_data,add_weight=FALSE,0.1,TSetKey)$detected_sig
    detected_orig<-which(allpairs %in% detected_orig_p)
    not_detected_orig<-which(!allpairs %in% detected_orig_p)
    kept_true_detected<-intersect(kept_pairs_true,detected_orig)
    kept_false_detected<-intersect(kept_pairs_false,detected_orig)
    kept_true_not_detected<-intersect(kept_pairs_true,not_detected_orig)
    kept_false_not_detected<-intersect(kept_pairs_false,not_detected_orig)
    plot((before[kept_false_not_detected]+after[kept_false_not_detected])/2,after[kept_false_not_detected]-before[kept_false_not_detected],main=paste0("seed",seed,"_FC",foldchange,"_nonsig_o_",meth[m]),pch = 46,xlab="log2(avg(original count)+1)/2+log2(avg(FC count)+1)/2",ylab="log2(avg(FC count)+1)-log2(avg(original count)+1)",ylim = c(-5,5),xlim=c(0,9))
    points((before[kept_false_detected]+after[kept_false_detected])/2,after[kept_false_detected]-before[kept_false_detected],col="blue",pch = 46)
    legend("topright", legend =c("TN","FP"),col=c("black","blue"),pch=19,cex=0.6)
    plot((before[kept_true_not_detected]+after[kept_true_not_detected])/2,after[kept_true_not_detected]-before[kept_true_not_detected],col="orange",main=paste0("seed",seed,"_FC",foldchange,"_truesig_o_",meth[m]),pch = 46,xlab="log2(avg(original count)+1)/2+log2(avg(FC count)+1)/2",ylab="log2(avg(FC count)+1)-log2(avg(original count)+1)",ylim = c(-5,5),xlim=c(0,9))
    points((before[kept_true_detected]+after[kept_true_detected])/2,after[kept_true_detected]-before[kept_true_detected],col="darkorange3",pch = 46)
    legend("topright", legend =c("TP","FN"),col=c("darkorange3","orange"),pch=19,cex=0.6)
    
  #}
  
    #for(m in 2:3){
    #plot weight
    detected_weighted_p<-get_fdr_dci(HiCbinpairs_data,add_weight=T,0.1,TSetKey)$detected_sig
    detected_weighted<-which(allpairs %in% detected_weighted_p)
    not_detected_weighted<-which(!allpairs %in% detected_weighted_p)
    kept_true_detected_w<-intersect(kept_pairs_true,detected_weighted)
    kept_true_new_detected_w<-setdiff(kept_true_detected_w,kept_true_detected)
    kept_false_detected_w<-intersect(kept_pairs_false,detected_weighted)
    kept_false_new_detected_w<-setdiff(kept_false_detected_w,kept_false_detected)
    kept_true_not_detected_w<-intersect(kept_pairs_true,not_detected_weighted)
    kept_false_not_detected_w<-intersect(kept_pairs_false,not_detected_weighted)
    plot((before[kept_false_not_detected_w]+after[kept_false_not_detected_w])/2,after[kept_false_not_detected_w]-before[kept_false_not_detected_w],main=paste0("seed",seed,"_FC",foldchange,"_nonsig_w_",meth[m]),pch = 46,xlab="log2(avg(original count)+1)/2+log2(avg(FC count)+1)/2",ylab="log2(avg(FC count)+1)-log2(avg(original count)+1)",ylim = c(-5,5),xlim=c(0,9))
    points((before[kept_false_detected_w]+after[kept_false_detected_w])/2,after[kept_false_detected_w]-before[kept_false_detected_w],col="blue",pch = 46)
    points((before[kept_false_new_detected_w]+after[kept_false_new_detected_w])/2,after[kept_false_new_detected_w]-before[kept_false_new_detected_w],col="cyan3",pch = 46)
    legend("topright", legend =c("TN","FP","FP_new"),col=c("black","blue","cyan3"),pch=19,cex=0.6)
    plot((before[kept_true_not_detected_w]+after[kept_true_not_detected_w])/2,after[kept_true_not_detected_w]-before[kept_true_not_detected_w],col="orange",main=paste0("seed",seed,"_FC",foldchange,"_truesig_w_",meth[m]),pch = 46,xlab="log2(avg(original count)+1)/2+log2(avg(FC count)+1)/2",ylab="log2(avg(FC count)+1)-log2(avg(original count)+1)",ylim = c(-5,5),xlim=c(0,9))
    points((before[kept_true_detected_w]+after[kept_true_detected_w])/2,after[kept_true_detected_w]-before[kept_true_detected_w],col="darkorange3",pch = 46)
    points((before[kept_true_new_detected_w]+after[kept_true_new_detected_w])/2,after[kept_true_new_detected_w]-before[kept_true_new_detected_w],col="darkred",pch = 46)
    legend("topright", legend =c("TP","FN","TP_new"),col=c("darkorange3","orange","darkred"),pch=19,cex=0.6)
    
    #}
    #plot tree
    in_path="/p/keles/fandingzhou/volumeA/TAD/GM12878"
    
    TADs<-read_tsv(paste0(in_path,"/OnTAD_seed",seed,"_GM_FC",foldchange,"_aggregated.tad"),col_names=F)
    HiCbinpairs_data<-tree_layer(HiCbinpairs_data,TADs,bin_size)
    max_layer<-max(TADs$X3)-3
    detected_tree_p<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=FALSE,0.1,TSetKey,max_layer)$detected_sig
    detected_tree<-which(allpairs %in% detected_tree_p)
    not_detected_tree<-which(!allpairs %in% detected_tree_p)
    kept_true_detected_t<-intersect(kept_pairs_true,detected_tree)
    kept_true_new_detected_t<-setdiff(kept_true_detected_t,kept_true_detected)
    kept_false_detected_t<-intersect(kept_pairs_false,detected_tree)
    kept_false_new_detected_t<-setdiff(kept_false_detected_t,kept_false_detected)
    kept_true_not_detected_t<-intersect(kept_pairs_true,not_detected_tree)
    kept_false_not_detected_t<-intersect(kept_pairs_false,not_detected_tree)
    plot((before[kept_false_not_detected_t]+after[kept_false_not_detected_t])/2,after[kept_false_not_detected_t]-before[kept_false_not_detected_t],main=paste0("seed",seed,"_FC",foldchange,"_nonsig_t_",meth[m]),pch = 46,xlab="log2(avg(original count)+1)/2+log2(avg(FC count)+1)/2",ylab="log2(avg(FC count)+1)-log2(avg(original count)+1)",ylim = c(-5,5),xlim=c(0,9))
    points((before[kept_false_detected_t]+after[kept_false_detected_t])/2,after[kept_false_detected_t]-before[kept_false_detected_t],col="blue",pch = 46)
    points((before[kept_false_new_detected_t]+after[kept_false_new_detected_t])/2,after[kept_false_new_detected_t]-before[kept_false_new_detected_t],col="cyan3",pch = 46)
    legend("topright", legend =c("TN","FP","FP_new"),col=c("black","blue","cyan3"),pch=19,cex=0.6)
    plot((before[kept_true_not_detected_t]+after[kept_true_not_detected_t])/2,after[kept_true_not_detected_t]-before[kept_true_not_detected_t],col="orange",main=paste0("seed",seed,"_FC",foldchange,"_truesig_t_",meth[m]),pch = 46,xlab="log2(avg(original count)+1)/2+log2(avg(FC count)+1)/2",ylab="log2(avg(FC count)+1)-log2(avg(original count)+1)",ylim = c(-5,5),xlim=c(0,9))
    points((before[kept_true_detected_t]+after[kept_true_detected_t])/2,after[kept_true_detected_t]-before[kept_true_detected_t],col="darkorange3",pch = 46)
    points((before[kept_true_new_detected_t]+after[kept_true_new_detected_t])/2,after[kept_true_new_detected_t]-before[kept_true_new_detected_t],col="darkred",pch = 46)
    legend("topright", legend =c("TP","FN","TP_new"),col=c("darkorange3","orange","darkred"),pch=19,cex=0.6)
    
    #plot wtree
    detected_wtree_p<-get_fdf_dci_tree(HiCbinpairs_data,add_weight=T,0.1,TSetKey,max_layer)$detected_sig
    detected_wtree<-which(allpairs %in% detected_wtree_p)
    not_detected_wtree<-which(!allpairs %in% detected_wtree_p)
    kept_true_detected_wt<-intersect(kept_pairs_true,detected_wtree)
    kept_true_new_detected_wt<-setdiff(kept_true_detected_wt,kept_true_detected_w)
    kept_false_detected_wt<-intersect(kept_pairs_false,detected_wtree)
    kept_false_new_detected_wt<-setdiff(kept_false_detected_wt,kept_false_detected_w)
    kept_true_not_detected_wt<-intersect(kept_pairs_true,not_detected_wtree)
    kept_false_not_detected_wt<-intersect(kept_pairs_false,not_detected_wtree)
    plot((before[kept_false_not_detected_wt]+after[kept_false_not_detected_wt])/2,after[kept_false_not_detected_wt]-before[kept_false_not_detected_wt],main=paste0("seed",seed,"_FC",foldchange,"_nonsig_wt_",meth[m]),pch = 46,xlab="log2(avg(original count)+1)/2+log2(avg(FC count)+1)/2",ylab="log2(avg(FC count)+1)-log2(avg(original count)+1)",ylim = c(-5,5),xlim=c(0,9))
    points((before[kept_false_detected_wt]+after[kept_false_detected_wt])/2,after[kept_false_detected_wt]-before[kept_false_detected_wt],col="blue",pch = 46)
    points((before[kept_false_new_detected_wt]+after[kept_false_new_detected_wt])/2,after[kept_false_new_detected_wt]-before[kept_false_new_detected_wt],col="cyan3",pch = 46)
    legend("topright", legend =c("TN","FP","FP_new"),col=c("black","blue","cyan3"),pch=19,cex=0.6)
    plot((before[kept_true_not_detected_wt]+after[kept_true_not_detected_wt])/2,after[kept_true_not_detected_wt]-before[kept_true_not_detected_wt],col="orange",main=paste0("seed",seed,"_FC",foldchange,"_truesig_wt_",meth[m]),pch = 46,xlab="log2(avg(original count)+1)/2+log2(avg(FC count)+1)/2",ylab="log2(avg(FC count)+1)-log2(avg(original count)+1)",ylim = c(-5,5),xlim=c(0,9))
    points((before[kept_true_detected_wt]+after[kept_true_detected_wt])/2,after[kept_true_detected_wt]-before[kept_true_detected_wt],col="darkorange3",pch = 46)
    points((before[kept_true_new_detected_wt]+after[kept_true_new_detected_wt])/2,after[kept_true_new_detected_wt]-before[kept_true_new_detected_wt],col="darkred",pch = 46)
    legend("topright", legend =c("TP","FN","TP_new"),col=c("darkorange3","orange","darkred"),pch=19,cex=0.6)
    
  }    
  
  
n=1  
seed=12
in_path="/p/keles/fandingzhou/volumeA/TAD/GM12878"

for(foldchange in c(4,6)){
  TADs<-read_tsv(paste0(in_path,"/OnTAD_seed",seed,"_GM_FC",foldchange,"_aggregated.tad"),col_names=F)
  
in_path2=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_",meth[n])
#in_path3=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/2vs2/filter_as_",meth[n])
for(wei in c("xi","pi0")){
for(m in 2:3){
  #plot weight
  in_path4=paste0(in_path2,"/",meth[m])
  HiCbinpairs_data<-read_tsv(paste0(in_path4,"/weighttable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_weightedp_",wei,".tsv"))
  HiCbinpairs_data<-na.omit(HiCbinpairs_data)
  Truesig2=scan(paste0(in_path2,"/Trueset/GMvsGMFC",foldchange,"_TSet"),what = character())
  TSetKey=data.frame(Truesig2,Truesig2)
   
  allpairs<-paste0(HiCbinpairs_data$bin_1,"_",HiCbinpairs_data$bin_2)
  #allpairs2<-sample(allpairs2,length(allpairs)/4)
  orig_p_neg_log10=-log10(HiCbinpairs_data$P)
  weighted_P<-ifelse(HiCbinpairs_data$weighted_P<=1,HiCbinpairs_data$weighted_P,1)
  weight_p_neg_log10=-log10(weighted_P)
  kept_pairs_true<-which(allpairs %in% Truesig2)
  kept_pairs_false<-which(!allpairs %in% Truesig2)
  HiCbinpairs_data<-tree_layer(HiCbinpairs_data,TADs,bin_size)
  max_layer<-max(TADs$X3)-3
  
  #p
  #detected_weighted_p<-get_fdr_dci(HiCbinpairs_data,add_weight=T,0.1,TSetKey)$detected_sig
   
  #wtree_p
  for(types in c("original","wtree")){
    if(types == "original"){
      detected_weighted_p<-get_fdr_dci(HiCbinpairs_data,add_weight=F,0.1,TSetKey,max_layer)$detected_sig
    }else{
      detected_weighted_p<-get_fdr_dci_tree(HiCbinpairs_data,add_weight=T,0.1,TSetKey,max_layer)$detected_sig
    }
    detected_weighted<-which(allpairs %in% detected_weighted_p)
    not_detected_weighted<-which(!allpairs %in% detected_weighted_p)
    kept_true_detected_w<-intersect(kept_pairs_true,detected_weighted)
    kept_false_detected_w<-intersect(kept_pairs_false,detected_weighted)
    kept_true_not_detected_w<-intersect(kept_pairs_true,not_detected_weighted)
    kept_false_not_detected_w<-intersect(kept_pairs_false,not_detected_weighted)
    plot(orig_p_neg_log10[kept_false_not_detected_w],weight_p_neg_log10[kept_false_not_detected_w],main=paste0("seed",seed,"_FC",foldchange,"_",meth[m],"_nonsig_",types),pch = 46,xlab="-log10(original_p)",ylab="-log10(weighted_p)",xlim = c(0,10),ylim=c(-2,10))
    points(orig_p_neg_log10[kept_false_detected_w],weight_p_neg_log10[kept_false_detected_w],col="blue",pch = 46)
    legend("topright", legend =c("TN","FP"),col=c("black","blue"),pch=19,cex=0.6)
    plot(orig_p_neg_log10[kept_true_not_detected_w],weight_p_neg_log10[kept_true_not_detected_w],main=paste0("seed",seed,"_FC",foldchange,"_",meth[m],"_truesig_",types),xlab="-log10(original_p)",ylab="-log10(weighted_p)",col="orange",pch = 46,xlim = c(0,10),ylim=c(-2,10))
    points(orig_p_neg_log10[kept_true_detected_w],weight_p_neg_log10[kept_true_detected_w],col="darkred",pch = 46)
    legend("topright", legend =c("FN","TP"),col=c("orange","darkred"),pch=19,cex=0.6)
    
  }
   
}  
}
}
q()
n
  trueindex = which(allpairs %in% Truesig) samptrue
  plot((before[sampfalse]+after[sampfalse])/2,after[sampfalse]-before[sampfalse],main=paste0("seed",seed,"_FC",foldchange),pch = 46,xlab="log2(avg(weightedinal count)+1)/2+log2(avg(FC count)+1)/2",ylab="log2(avg(FC count)+1)-log2(avg(original count)+1)",ylim = c(-5,5))
  

