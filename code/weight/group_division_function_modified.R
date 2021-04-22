group_division<-function(HiCbinpairs_data,boundary1,groups.by="compartment",clustern=NULL,small_size=20,alpha=600,iter=20){
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
    original_group=HiCbinpairs_data$Group=rowindex*1000+colindex
    
    #clustering
    #table(HiCbinpairs_data$Group)[1:100]
    for(k in 1:iter){
      HiCbinpairs_data_group =
        HiCbinpairs_data %>% 
        group_by(Group) %>% 
        dplyr::count(Group) %>% 
        ungroup
      
      small_group<-as.vector(HiCbinpairs_data_group$Group[HiCbinpairs_data_group$n<small_size])
      large_group<-unique(original_group[!original_group %in% small_group])
      #large_group2=as.vector(HiCbinpairs_data_group$Group[HiCbinpairs_data_group$n>=small_size])
      n=0
      for(i in small_group){
        surroundings<-c(i+1000,i-1000,i+1,i-1)
        nearby_group<-surroundings[surroundings %in% large_group]
        if(length(nearby_group)==0){
          surroundings<-c(i+1001,i-1001,i+999,i-999)
          nearby_group<-surroundings[surroundings %in% large_group]
        }
        if(length(nearby_group)>0){
          n=n+1
          #mapping original groups to previously merged groups
          nearby_group2=unique(HiCbinpairs_data$Group[original_group %in% nearby_group])
          #find the smallest adjacent large group
          smallest_nearby_Lgroup_n<-min(HiCbinpairs_data_group$n[HiCbinpairs_data_group$Group %in% nearby_group2])
          smallest_nearby_Lgroup<-HiCbinpairs_data_group$Group[HiCbinpairs_data_group$n==smallest_nearby_Lgroup_n][1]
          #set small group index to be larger group index
          HiCbinpairs_data$Group[HiCbinpairs_data$Group==i]=smallest_nearby_Lgroup
        }
      }
      if(n<10){
        break
      }
    }
    if(is.null(clustern)){
      avg_small<-mean(HiCbinpairs_data_group$n[HiCbinpairs_data_group$Group %in% small_group]) 
      clustern=floor(length(small_group)*avg_small/(small_size*4))
    }
    rowind<-small_group%/%1000
    colind<-small_group%%1000 
    kmeandata<-data.frame(cbind(rowind,colind))
    result<-kmeans(kmeandata,clustern) 
    small_group<-data.frame(Group=small_group,result.cluster=result$cluster)
    for (i in 1:clustern){
      for (j in small_group$Group[small_group$result.cluster==i]){
        HiCbinpairs_data$Group[HiCbinpairs_data$Group==j]=i
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
    second_large_dist<-sort(as.numeric(names(table(HiCbinpairs_data$Group))),decreasing = T)[2]
    HiCbinpairs_data$Group[HiCbinpairs_data$Group==max(HiCbinpairs_data$Group)]=second_large_dist
    
  }
  return(HiCbinpairs_data)
}
