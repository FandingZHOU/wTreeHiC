get_fdf_dci_tree<-function(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="DCI",max_layer){
  if(add_weight==T){
    pvalue="weighted_P"
  }else{
    pvalue="P"
  }
  FDR_tree_detected<-c()
  DCI_tree_true<-c()
  power_tree<-c()
  for (j in 1:length(threshold)){
    for(lay in 1:max_layer){
      HiCbinpairs_data<-tree_filter_group(HiCbinpairs_data,paste0("Group.",lay),pvalue,threshold[j])
    }
    #final L
    HiCbinpairs_data<-HiCbinpairs_data[order(HiCbinpairs_data[[pvalue]]),]
    pL<-HiCbinpairs_data[[pvalue]]
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
    DCI_tree_true[j]<-length(unique(TSetKey[TSetKey$smooth %in% datapairs,]$origin))
    #DCI_true[j]<-length(intersect(datapairs,TSetKey$smooth))
    FDR_tree_detected[j]<-1-length(intersect(TSetKey$smooth, datapairs))/detected_sig_n
    power_tree[j]<-DCI_tree_true[j]/length(unique(TSetKey$origin))
  }
  if(gettype=="DCI"){
    return(DCI_tree_true)
  }else{
    if(gettype=="FDR"){
      return(FDR_tree_detected)
    }else{
      return(power_tree)
    }
    
  }
}


