#additional function for pi0
pi0est1<-function(x){
  i=0
  x<-x[!is.na(x)]
  if (length(x)<20){
    0.9
  }else{
    if(max(x)<=0.95){
      0.8
    }else{
      a<-min(table(cut(x,seq(0,1,0.05))))
      if(a==0){
        0.95
      }else{
        tryCatch({
          pi0est(x,lambda=seq(0.05, 0.95, 0.05))$pi0
        },error = function(e) {
          0.9
        })
        
      }
    }
  }
}

# for(n in seq(0.05,0.95,0.05)){
#   if(length(x[n-0.05<=x&x<n])==0){
#     i=88
#     break
#   }
# }
# if(i==88){
#   0.95
# }

add_weight<-function(HiCbinpairs_data,weight_method,gamma=0.05,statistic="F",type="table"){
  if(weight_method=="xi"){
    if(statistic=="F"){
      HiCstat<-HiCbinpairs_data$F
    }else{
      #HiCstat<-HiCbinpairs_data$Z
      HiCstat<-HiCbinpairs_data$Z^2
    }
    Y=aggregate(x = HiCstat, by= list(HiCbinpairs_data$Group), FUN = mean)
    var=aggregate(x = HiCstat, by= list(HiCbinpairs_data$Group), FUN = var)
    rk=aggregate(x = HiCstat, by= list(HiCbinpairs_data$Group), FUN = length)
    n_group<-nrow(Y)
    xisq=c()
    pik<-c()
    xi=c()
    #if(statistic=="F"){
      for (i in 1:n_group){
        if(Y$x[i]<=1 | rk$x[i]<=1){
          xisq[i]<-0
          pik[i]<-0
        }else{
          xisq[i]<-(var$x[i]+Y$x[i]^2+3)/(Y$x[i]-1)
          pik[i]<-(Y$x[i]-1)/xisq[i]
          if(1/rk$x[i] >= pik[i] | pik[i]>= (rk$x[i]-1)/rk$x[i]){ #
            xisq[i]=0
         
          }
        }
      }
      xi<-sqrt(xisq)
    #}else{
      
      # for (i in 1:n_group){
      #   if(rk$x[i]<=1){
      #     xi[i]<-0
      #     pik[i]<-0
      #   }else{
      #     pik[i]<-(Y$x[i]^2)/(Y$x[i]^2+var$x[i]-1)
      #     xi[i]<-abs(Y$x[i]/pik[i])
      #     if(1/rk$x[i] >= pik[i]){
      #       xi[i]=0
      #     }
      #   }
      # }
    #}
    
    weight<-(1-gamma)*xi+gamma*mean(xi)
    weight<-weight/sum(weight*rk$x)*nrow(HiCbinpairs_data)
    head(HiCbinpairs_data[HiCbinpairs_data$weight<0.05,])
    #weight<-xi/sum(xi*rk$x)*nrow(HiCbinpairs_data)
    #weight<-(1-gamma)*weight+gamma
    weight<-data.frame(cbind(rk,weight))
    for (i in 1:n_group){
      HiCbinpairs_data$weight[HiCbinpairs_data$Group==weight$Group.1[i]]<-rep(weight$weight[i],weight$x[i])
    }
    weighted_P<-HiCbinpairs_data$P/HiCbinpairs_data$weight
  }

  
  
  if(weight_method=="pi0"){
    rk=aggregate(x = HiCbinpairs_data$P, by= list(HiCbinpairs_data$Group), FUN = length)
    nullratio<-aggregate(x = HiCbinpairs_data$P, by= list(HiCbinpairs_data$Group), FUN = pi0est1)
    pi0_pre<-nullratio$x
    weight<-(1-pi0_pre)/pi0_pre
    weight<-(1-gamma)*weight+gamma*mean(weight)
    weight<-weight/sum(weight*rk$x)*nrow(HiCbinpairs_data)
    weight<-data.frame(cbind(rk$Group.1,rk$x,weight))
    for (i in 1:nrow(rk)){
      HiCbinpairs_data$weight[HiCbinpairs_data$Group==weight$V1[i]]<-rep(weight$weight[i],weight$V2[i])
    }
    weighted_P<-HiCbinpairs_data$P/HiCbinpairs_data$weight
    
  }
  
  if(weight_method=="IHW"){
    HiCbinpairs_data$Group<-as.numeric(HiCbinpairs_data$Group)
    second_large_dist<-sort(as.numeric(names(table(HiCbinpairs_data$Group))),decreasing = T)[2]
    HiCbinpairs_data$Group[HiCbinpairs_data$Group==max(HiCbinpairs_data$Group)]=second_large_dist
    HiCbinpairs_data$Group<-as.factor(HiCbinpairs_data$Group)
    ihw_object<-ihw(P~Group, data=HiCbinpairs_data,alpha=0.1)#, covariate_type = "nominal")
    HiCbinpairs_data$weight<-weights(ihw_object)
    weighted_P<-weighted_pvalues(ihw_object)
  }
  
  weighted_P[weighted_P>1]<-1
  if(type=="table"){
    HiCbinpairs_data<-data.frame(cbind(HiCbinpairs_data,weighted_P))
    HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data$weighted_P),])
    adj_w_p<-p.adjust(HiCbinpairs_data$weighted_P,"BH")
    HiCbinpairs_data<-data.frame(cbind(HiCbinpairs_data,adj_w_p))
    return(HiCbinpairs_data)
  }else{
    return(weighted_P)
  }
  
}
