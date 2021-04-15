tree_layer<-function(HiCbinpairs_data,TADs,bin_size){
  Group2<-HiCbinpairs_data$rowindex*1000+HiCbinpairs_data$colindex+0.5
  layer<-rep(0,nrow(HiCbinpairs_data))
  Lgroup<-rep(0,nrow(HiCbinpairs_data))
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
    # if(la==1){
    # Group[Group==0]<-HiCbinpairs_data$Group[Group==0]+0.5 #first layer: clustered group
    # }else{
      Group[Group==0]<-Group2[Group==0] #other layers except final layer: original group
    #}
    HiCbinpairs_data<-data.frame(cbind(HiCbinpairs_data,Group))
  }
  return(HiCbinpairs_data)
}