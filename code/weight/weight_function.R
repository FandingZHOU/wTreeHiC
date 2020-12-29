#additional function for pi0
pi0est1<-function(x){
  i=0
  if (length(x)<20){
    0.95
  }else{
    if(max(x)<=0.95){
      0.9
    }else{
      a<-min(table(cut(x,seq(0,1,0.05))))
      if(a==0){
        0.99
      }else{
        pi0est(x,lambda=seq(0.05, 0.95, 0.05))$pi0
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

add_weight<-function(HiCbinpairs_data,weight_method,gamma=0.05){
  if(weight_method=="xi"){
    Y=aggregate(x = HiCbinpairs_data$F, by= list(HiCbinpairs_data$Group), FUN = mean)
    var=aggregate(x = HiCbinpairs_data$F, by= list(HiCbinpairs_data$Group), FUN = var)
    rk=aggregate(x = HiCbinpairs_data$F, by= list(HiCbinpairs_data$Group), FUN = length)
    n_group<-nrow(Y)
    xisq=c()
    pik<-c()
    for (i in 1:n_group){
      if(Y$x[i]<=1 | rk$x[i]<=1){
        xisq[i]<-0
        pik[i]<-0
      }
      else{
        xisq[i]<-(var$x[i]+Y$x[i]^2+3)/(Y$x[i]-1)
        pik[i]<-(Y$x[i]-1)/xisq[i]
        if(1/rk$x[i] > pik[i] | pik[i]> (rk$x[i]-1)/rk$x[i]){
          xisq[i]=0
        }
      }
    }
    xi<-sqrt(xisq)
    weight<-xi/sum(xi*rk$x)*nrow(HiCbinpairs_data)
    weight<-(1-gamma)*weight+gamma
    weight<-data.frame(cbind(rk,weight))
    for (i in 1:n_group){
      HiCbinpairs_data$weight[HiCbinpairs_data$Group==weight$Group.1[i]]<-rep(weight$weight[i],weight$x[i])
    }
    weighted_P<-HiCbinpairs_data$P/HiCbinpairs_data$weight
    
    
    ##plot density
    # chosen_group<-rk[rk$x>600,]$Group.1
    # chosen_xisq<-xisq[rk$x>600]
    # chosen_pi<-pik[rk$x>600]
    # 
    # chosen_group<-rk[rk$x<30&rk$x>=20,]$Group.1
    # chosen_xisq<-xisq[rk$x<30&rk$x>=20]
    # chosen_pi<-pik[rk$x<30&rk$x>=20]
    # 
    # 
    # plot_F_density<-function(chosen_group,chosen_xisq,chosen_pi){
    #   for (b in 1:length(chosen_group)){
    #     n=1000
    #     temp = cbind(rchisq(n,1),(rnorm(n,sqrt(chosen_xisq[b]),1))^2)
    #     id = sample(1:2,n,rep = T,prob = c((1-chosen_pi[b]),chosen_pi[b]))  
    #     id = cbind(1:n,id)
    #     theory_distribution<-temp[id]
    #     Fvalue<-HiCbinpairs_data$F[HiCbinpairs_data$Group==chosen_group[b]]
    #     hist(Fvalue,
    #          freq=F,
    #          breaks=30)
    #     lines(density(theory_distribution),col="red")
    #   }
    # }
  }
  
  
  if(weight_method=="pi0"){
    rk=aggregate(x = HiCbinpairs_data$P, by= list(HiCbinpairs_data$Group), FUN = length)
    nullratio<-aggregate(x = HiCbinpairs_data$P, by= list(HiCbinpairs_data$Group), FUN = pi0est1)
    pi0_pre<-nullratio$x
    weight<-(1-pi0_pre)/pi0_pre
    weight<-weight/sum(weight*rk$x)*nrow(HiCbinpairs_data)
    weight<-(1-gamma)*weight+gamma
    weight<-data.frame(cbind(rk$Group.1,rk$x,weight))
    
    for (i in 1:nrow(rk)){
      HiCbinpairs_data$weight[HiCbinpairs_data$Group==weight$V1[i]]<-rep(weight$weight[i],weight$V2[i])
    }
    weighted_P<-HiCbinpairs_data$P/HiCbinpairs_data$weight
  }
  
  adj_w_p<-c()
  HiCbinpairs_data<-data.frame(cbind(HiCbinpairs_data,weighted_P))
  HiCbinpairs_data<-data.frame(HiCbinpairs_data[order(HiCbinpairs_data$weighted_P),])
  for (i in 1:nrow(HiCbinpairs_data)){
    adj_w_p[i]<-nrow(HiCbinpairs_data)*HiCbinpairs_data$weighted_P[i]/i
    if (adj_w_p[i]>1){
      adj_w_p[i]=1
    }
  }
  HiCbinpairs_data<-data.frame(cbind(HiCbinpairs_data,adj_w_p))
  return(HiCbinpairs_data)
}
