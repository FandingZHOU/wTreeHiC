tree_filter_group<-function(HiCbinpairs_data,Group_name,pvalue,thres){
  pL1=aggregate(x = HiCbinpairs_data[[pvalue]], by= list(HiCbinpairs_data[[Group_name]]), FUN = min)
  pL1<-pL1[order(pL1$x),]
  adj_pL1<-c()
  for (i in 1:length(pL1$x)){
    adj_pL1[i]<-length(pL1$x)*pL1$x[i]/i
  }
  maxpi<-max(pL1$x[adj_pL1<=thres])
  sig_group1=pL1[pL1$x<=maxpi,]
  HiCbinpairs_data<-HiCbinpairs_data[HiCbinpairs_data[[Group_name]] %in% sig_group1$Group.1,]
  return(HiCbinpairs_data)
}