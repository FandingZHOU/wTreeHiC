get_dci_tree<-function(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,weight_method="xi",gamma=0.05){
  if(add_weight==T){
    pvalue="weighted_P"
  }else{
    pvalue="P"
  }
  datapairs<-list()
  for (j in 1:length(threshold)){
  for(lay in 1:max_layer){
      HiCbinpairs_data<-tryCatch(
        tree_filter_group(HiCbinpairs_data,paste0("Group.",lay),pvalue,threshold,weight_method,gamma),
        error=function(e){
          return(HiCbinpairs_data[1,])
        }
      )
    }
    if(nrow(HiCbinpairs_data)==1){
      datapairs1=NULL
    }else{
    #final L
    HiCbinpairs_data<-HiCbinpairs_data[order(HiCbinpairs_data[[pvalue]]),]
    pL<-HiCbinpairs_data[[pvalue]]
    adj_pL<-c()
    for (i in 1:length(pL)){
      adj_pL[i]<-length(pL)*pL[i]/i
    }
    maxpi<-max(pL[adj_pL<=threshold])
    sig_pair=HiCbinpairs_data[pL<=maxpi,]
    
    #Compare
    datapairs1<-c()
    detected_sig_n<-nrow(sig_pair)
    for (i in 1:detected_sig_n){
      datapairs1[i]<-paste0(sig_pair[i,1],"_",sig_pair[i,2])
    }
    }
    datapairs[[j]]=datapairs1
    names(datapairs)=threshold
  }
  
  return(datapairs)
}
obtain_power_fdr<-function(datapairs,TSetKey){
  result<-list()
  if(is.null(ncol(TSetKey))){
    TSetKey1=data.frame()
    TSetKey1$smooth=TSetKey1$origin=TSetKey
    TSetKey=TSetKey1
  }
  if(!is.list(datapairs)&is.vector(datapairs)){
    datapairs1=list()
    datapairs1[[1]]=datapairs
    datapairs=datapairs1
  }else{
    result[["threshold"]]<-as.numeric(names(datapairs))
  }
  FDR_tree_detected<-c()
  DCI_tree_true<-c()
  power_tree<-c()
  for(j in 1:length(datapairs)){
    DCI_tree_true[j]<-length(unique(TSetKey[TSetKey$smooth %in% datapairs[[j]],]$origin))
    FDR_tree_detected[j]<-1-length(intersect(TSetKey$smooth, datapairs[[j]]))/detected_sig_n
    power_tree[j]<-DCI_tree_true[j]/length(unique(TSetKey$origin))
  }
  # result[["detected_sig"]]<-datapairs
  result[["FDR"]]<-FDR_tree_detected
  result[["DCI"]]<-DCI_tree_true
  result[["power"]]<-power_tree
  
  return(result)
}

get_fdr_dci_tree<-function(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,max_layer,weight_method="xi",gamma=0.05){
  if(add_weight==T){
    pvalue="weighted_P"
  }else{
    pvalue="P"
  }
  result<-list()
  FDR_tree_detected<-c()
  DCI_tree_true<-c()
  power_tree<-c()
  for (j in 1:length(threshold)){
    for(lay in 1:max_layer){
      
      HiCbinpairs_data<-tryCatch(
        tree_filter_group(HiCbinpairs_data,paste0("Group.",lay),pvalue,threshold[j],weight_method,gamma),
        error=function(e){
          return(HiCbinpairs_data[1,])
        }
      )
      }
    if(nrow(HiCbinpairs_data)==1){
      DCI_tree_true[j]<-0
      #DCI_true[j]<-length(intersect(datapairs,TSetKey$smooth))
      FDR_tree_detected[j]<-0
      power_tree[j]<-0
      next
      #len_thre=length(threshold)
      #result=list(FDR=rep(0,len_thre),DCI=rep(0,len_thre),power=rep(0,len_thre))
      #return(result)
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
  #result[["detected_sig"]]<-datapairs
  result[["FDR"]]<-FDR_tree_detected
  result[["DCI"]]<-DCI_tree_true
  result[["power"]]<-power_tree
  return(result)
  # if(gettype=="DCI"){
  #   return(DCI_tree_true)
  # }else{
  #   if(gettype=="FDR"){
  #     return(FDR_tree_detected)
  #   }else{
  #     return(power_tree)
  #   }
  #   
  # }
}


