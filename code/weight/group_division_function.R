group_division<-function(HiCbinpairs_data,boundary1,size="regular",groups.by="compartment",clustern=100,large_size=50,small_size=20,alpha=600){
  HiCbinpairs_data<-data.frame(HiCbinpairs_data)
  if(is.null(HiCbinpairs_data$weight)|is.null(HiCbinpairs_data$Group)){
    HiCbinpairs_data$weight<-rep(0,nrow(HiCbinpairs))
    HiCbinpairs_data$Group<-rep(0,nrow(HiCbinpairs))
  }
  
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data[,1]),])
    rowindex<-as.numeric(cut(HiCbinpairs_data[,1],boundary1))
    colindex<-as.numeric(cut(HiCbinpairs_data[,2],boundary1))
    HiCbinpairs_data<-data.frame(cbind(HiCbinpairs_data,rowindex,colindex))
  if(groups.by=="compartment"){
    HiCbinpairs_data$Group<-rowindex*1000+colindex
    
    #clustering
    #table(HiCbinpairs_data$Group)[1:100]
    HiCbinpairs_data_group =
      HiCbinpairs_data %>% 
      group_by(Group) %>% 
      dplyr::count(Group) %>% 
      ungroup
    
    ##large/regular group
    if(size=="large"){
      small_group<-HiCbinpairs_data_group[HiCbinpairs_data_group$n<large_size,]
    }else{
      small_group<-HiCbinpairs_data_group[HiCbinpairs_data_group$n<small_size,]
      large_group<-HiCbinpairs_data_group[HiCbinpairs_data_group$n>large_size,]
    }
    rowind<-small_group$Group%/%1000
    colind<-small_group$Group%%1000 
    kmeandata<-data.frame(cbind(rowind,colind))
    result<-kmeans(kmeandata,clustern) 
    small_group<-data.frame(cbind(small_group,result$cluster)) 
    for (i in 1:clustern){
      for (j in small_group$Group[small_group$result.cluster==i]){
        HiCbinpairs_data$Group[HiCbinpairs_data$Group==j]=i
      }
    }
    
    #small group
    if(size=="small"){
      n_large_group=nrow(large_group)
      for(j in 1:n_large_group){
        tim<-floor(large_group$n[j]/20)-1
        a<-c()
        for(t in 1:tim){
          a<-c(a,rep(0.001*t,20))
        }
        a<-c(a,rep(0,large_group$n[j]-length(a)))
        HiCbinpairs_data$Group[HiCbinpairs_data$Group==large_group$Group[j]]=a+large_group$Group[j]
      }
    }
    
  }
  
  if(groups.by=="distance"){
    D=(HiCbinpairs_data$bin_2-HiCbinpairs_data$bin_1)/bin_size
    D.table=table(D)
    merge_begin=as.numeric(names(D.table[min(which(D.table<large_size))]))
    merge_D_seq=1:floor((4*alpha*(max(D)))^(1/3))
    merge_D_seq=merge_begin+merge_D_seq*(merge_D_seq+1)*(merge_D_seq*2+1)/alpha
    Group<-as.numeric(cut(D,c(1:merge_begin,merge_D_seq)))
    HiCbinpairs_data$Group=Group
    return(HiCbinpairs_data)
  }
}
