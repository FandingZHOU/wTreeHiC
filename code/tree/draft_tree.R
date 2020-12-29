#give tads
Group2<-HiCbinpairs_data$rowindex*1000+HiCbinpairs_data$colindex+0.5
layer<-rep(0,nrow(HiCbinpairs_data))
Lgroup<-rep(0,nrow(HiCbinpairs_data))
TADs<-read_tsv(paste0(in_path,"/OnTAD_GM_spike_in_seqx",seqdepth[n],"_aggregated.tad"),col_names=F)
TADs$X1<-(TADs$X1-1)*bin_size
TADs$X2<-TADs$X2*bin_size
max_layer<-max(TADs$X3)
for(la in 1:max_layer){
  TADs_layer<-TADs[TADs$X3==la,]
  for(Lg in 1:nrow(TADs_layer)){
    chosen_index<-which(HiCbinpairs_data$bin_1>TADs_layer$X1[Lg]&HiCbinpairs_data$bin_2<TADs_layer$X2[Lg])
    layer[chosen_index]<-la
    Lgroup[chosen_index]<-Lg
  }
  Group<-layer*1000+Lgroup
  if(la==1){
    Group[Group==0]<-HiCbinpairs_data$Group[Group==0]+0.5
  }else{
    Group[Group==0]<-Group2[Group==0]
  }
  HiCbinpairs_data<-data.frame(cbind(HiCbinpairs_data,Group))
}
original_data<-HiCbinpairs_data
#HiCbinpairs_data<-data.frame(cbind(HiCbinpairs_data,layer, Lgroup))
#table(HiCbinpairs_data$Group.1)[1:20]

FDR_tree_detected<-c()
DCI_tree_true<-c()
HiCbinpairs_data<-original_data
for (j in 1:length(threshold)){
  for(lay in 1:max_layer){
    HiCbinpairs_data<-tree_filter_group(HiCbinpairs_data,paste0("Group.",lay),"P",threshold[j])
  }
  #final L
  HiCbinpairs_data<-HiCbinpairs_data[order(HiCbinpairs_data$P),]
  pL<-HiCbinpairs_data$P
  adj_pL<-c()
  for (i in 1:length(pL)){
    adj_pL[i]<-length(pL)*pL[i]/i
  }
  maxpi<-max(pL[adj_pL<=threshold[j]])
  sig_pair=HiCbinpairs_data[pL<=maxpi,]
  
  #Compare
  datapairs<-c()
  detected_sig_n<-nrow(sig_pair)
  for (i in 1:detected_sig_n){
    datapairs[i]<-paste0(sig_pair[i,1],"_",sig_pair[i,2])
  }
  DCI_tree_true[j]<-length(intersect(datapairs,TSetKey))
  FDR_tree_detected[j]<-1-DCI_tree_true[j]/detected_sig_n
}


FDR_wtree_detected<-c()
DCI_wtree_true<-c()
HiCbinpairs_data<-original_data
for (j in 1:length(threshold)){
  for(lay in 1:max_layer){
    HiCbinpairs_data<-tree_filter_group(HiCbinpairs_data,paste0("Group.",lay),"weighted_P",threshold[j])
  }
  #final L
  HiCbinpairs_data<-HiCbinpairs_data[order(HiCbinpairs_data$weighted_P),]
  pL<-HiCbinpairs_data$weighted_P
  adj_pL<-c()
  for (i in 1:length(pL)){
    adj_pL[i]<-length(pL)*pL[i]/i
  }
  maxpi<-max(pL[adj_pL<=threshold[j]])
  sig_pair=HiCbinpairs_data[pL<=maxpi,]
  
  #Compare
  datapairs<-c()
  detected_sig_n<-nrow(sig_pair)
  for (i in 1:detected_sig_n){
    datapairs[i]<-paste0(sig_pair[i,1],"_",sig_pair[i,2])
  }
  DCI_wtree_true[j]<-length(intersect(datapairs,TSetKey))
  FDR_wtree_detected[j]<-1-DCI_wtree_true[j]/detected_sig_n
}

