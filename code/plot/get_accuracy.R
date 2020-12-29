get_accuracy<-function(HiCbinpairs_data,add_weight=FALSE,topN,interval_thre=500,TSetKey){
  thre<-seq(0,topN,interval_thre)
  datapairs<-c()
  
  for (i in 1:nrow(HiCbinpairs_data)){
    datapairs[i]<-paste0(HiCbinpairs_data[i,1],"_",HiCbinpairs_data[i,2])
  }
  accuracy<-c()
  if(add_weight==TRUE){
    datapairs<-datapairs[order(HiCbinpairs_data$weighted_P)]
  }else{
    datapairs<-datapairs[order(HiCbinpairs_data$P)]
  }
  sum_sig<-0
  for(i in 1:(length(thre)-1)){
    sig_pairs<-datapairs[(thre[i]+1):thre[i+1]]
    dupl<-intersect(sig_pairs,TSetKey)
    sum_sig<-sum_sig+length(dupl)
    accuracy[i]<-sum_sig/thre[i+1]
  }
  return(accuracy)
}
