
##without having True signals
get_dci<-function(HiCbinpairs_data,add_weight=FALSE,threshold){
  if(add_weight==T){
    pvalue="weighted_P"
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data[[pvalue]]),])
    adj_p<-HiCbinpairs_data$adj_w_p
  }else{
    pvalue="P"
    adj_p<-c()
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data[[pvalue]]),])
    for (i in 1:nrow(HiCbinpairs_data)){
      adj_p[i]<-nrow(HiCbinpairs_data)*HiCbinpairs_data$P[i]/i
    }
  }
  maxpi<-max(HiCbinpairs_data[adj_p<=threshold[j],][[pvalue]])
  detected_sig_data=HiCbinpairs_data[HiCbinpairs_data[[pvalue]]<=maxpi,]
  datapairs<-c()
  detected_sig_n<-nrow(detected_sig_data)
  for (i in 1:detected_sig_n){
      datapairs[i]<-paste0(detected_sig_data[i,1],"_",detected_sig_data[i,2])
  }
  return(datapairs)
}

#with true signals
get_fdr_dci<-function(HiCbinpairs_data,add_weight=FALSE,threshold,TSetKey,gettype="DCI"){
  if(add_weight==T){
    pvalue="weighted_P"
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data[[pvalue]]),])
    adj_p<-HiCbinpairs_data$adj_w_p
  }else{
    pvalue="P"
    adj_p<-c()
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data[[pvalue]]),])
    for (i in 1:nrow(HiCbinpairs_data)){
      adj_p[i]<-nrow(HiCbinpairs_data)*HiCbinpairs_data$P[i]/i
    }
  }
  result<-list()
  FDR_detected<-c()
  DCI_true<-c()
  power_detected<-c()
  for (j in 1:length(threshold)){
    maxpi<-max(HiCbinpairs_data[adj_p<=threshold[j],][[pvalue]])
    detected_sig_data=HiCbinpairs_data[HiCbinpairs_data[[pvalue]]<=maxpi,]
    datapairs<-c()
    detected_sig_n<-nrow(detected_sig_data)
    for (i in 1:detected_sig_n){
      datapairs[i]<-paste0(detected_sig_data[i,1],"_",detected_sig_data[i,2])
    }
    DCI_true[j]<-length(unique(TSetKey[TSetKey$smooth %in% datapairs,]$origin))
    #DCI_true[j]<-length(intersect(datapairs,TSetKey$smooth))
    FDR_detected[j]<-1-length(intersect(TSetKey$smooth, datapairs))/detected_sig_n
    power_detected[j]<-DCI_true[j]/length(unique(TSetKey$origin))
  }
  #result[["detected_sig"]]<-datapairs
  result[["FDR"]]<-FDR_detected
  result[["DCI"]]<-DCI_true
  result[["power"]]<-power_detected
  return(result)
  # if(gettype=="DCI"){
  #   return(DCI_true)
  # }
  # if(gettype=="FDR"){
  #   return(FDR_detected)
  # }
  # if(gettype=="power"){
  #   return(power_detected)
  # }
}
