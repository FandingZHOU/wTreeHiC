group_division<-function(HiCbinpairs_data,boundary1,size="regular",clustern=100,large_size=50,small_size=20){
  
  rowindex<-as.numeric(cut(HiCbinpairs_data[,1],boundary1))
  colindex<-as.numeric(cut(HiCbinpairs_data[,2],boundary1))
  HiCbinpairs_data<-data.frame(cbind(HiCbinpairs_data,rowindex,colindex))
  HiCbinpairs_data$Group<-rowindex*1000+colindex
  
  #clustering
  #table(HiCbinpairs_data$Group)[1:100]
  HiCbinpairs_data_group =
    HiCbinpairs_data %>% 
    group_by(Group) %>% 
    count %>% 
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
  return(HiCbinpairs_data)
}
