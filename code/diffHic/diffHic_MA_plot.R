diffHic_MA_plot<-function(data,nrep,bin_size,scalesize,TSetKey,seqdepth){
  useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
  inter.frame <- as.data.frame(interactions(data))[,useful.cols]
  inter_pairs<-paste0(inter.frame[,6]-bin_size/2,"_",inter.frame[,3]-bin_size/2)
  truesig<-which(inter_pairs %in% TSetKey)
  
  ab<-aveLogCPM(asDGEList(data))
  # Order avg logCPM
  o<-order(ab)
  # Calculate counts per million
  adj.counts<-cpm(asDGEList(data),log=TRUE)
  # Calculate M
  mval<- adj.counts[,6]-adj.counts[,2]
  #plot(ab, mval,xlab="A",ylab="M",main="Spike-in GM2 vs. Original GM2",pch=16)
  smoothScatter(ab, mval,xlab="A",ylab="M",main=paste0("no filter Spike-in vs. Original scale",scalesize,"seqx",seqdepth))
  points(ab[truesig],mval[truesig],col="orange",pch = 46)
  fit<-loessFit(x=ab,y=mval)# Add loess fit to MA plot
  lines(ab[o], fit$fitted[o],col="red")
  
  #filtering
  dist.keep <-  pairdist(data)>2*bin_size
  data<-data[dist.keep,]
  inter_pairs<-inter_pairs[dist.keep]
  ave.ab <- aveLogCPM(asDGEList(data))
  count.keep <- ave.ab>= aveLogCPM(2,lib.size=mean(data$totals))
  data<-data[count.keep,]
  inter_pairs<-inter_pairs[count.keep]
  truesig<-which(inter_pairs %in% TSetKey)
  
  # Calculate A
  ab<-aveLogCPM(asDGEList(data))
  # Order avg logCPM
  o<-order(ab)
  # Calculate counts per million
  adj.counts<-cpm(asDGEList(data),log=TRUE)
  # Calculate M
  mval<- adj.counts[,6]-adj.counts[,2]
  #plot(ab, mval,xlab="A",ylab="M",main="Spike-in GM2 vs. Original GM2",pch=16)
  smoothScatter(ab, mval,xlab="A",ylab="M",main=paste0("Spike-in vs. Original scale",scalesize,"seqx",seqdepth))
  points(ab[truesig],mval[truesig],col="orange",pch = 46)
  fit<-loessFit(x=ab,y=mval)# Add loess fit to MA plot
  lines(ab[o], fit$fitted[o],col="red")
  
 #  #Normal
 #  data <- normOffsets(data, method="loess", se.out=TRUE)
 # # adj.counts<-cpm(asDGEList(data),log=TRUE)
 #  adj.counts<-log2(assay(data)+0.5)- assay(data,"offset")/log(2)
 #  mval<- adj.counts[,6]-adj.counts[,2]
 #   smoothScatter(ab, mval,xlab="A",ylab="M",main=paste0("Normalized Spike-in vs. Original scale",scalesize))
 #  fit<-loessFit(x=ab,y=mval)# Add loess fit to MA plot
 #  lines(ab[o], fit$fitted[o],col="red")
}

diffHic_MA_box_plot<-function(data,nrep,bin_size,scalesize,TSetKey,seqdepth){
  useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
  inter.frame <- as.data.frame(interactions(data))[,useful.cols]
  inter_pairs<-paste0(inter.frame[,6]-bin_size/2,"_",inter.frame[,3]-bin_size/2)
  truesig<-rep(0,length(inter_pairs))
  truesig[inter_pairs %in% TSetKey]<-1
  truesig<-as.factor(truesig)
  
  ab<-aveLogCPM(asDGEList(data))
  # Order avg logCPM
  o<-order(ab)
  # Calculate counts per million
  adj.counts<-cpm(asDGEList(data),log=TRUE)
  # Calculate M
  if(abslt){
    mval<- abs(adj.counts[,6]-adj.counts[,2])
  }else{
    mval<- adj.counts[,6]-adj.counts[,2]
  }
  
  ab.stratify<-cut(ab,c(-4,-2,0,2,4,7))
  true_vs_nontrue<-data.frame(ab.stratify,mval,truesig)
  p1<-ggplot(true_vs_nontrue, aes(truesig, mval,fill=ab.stratify)) + #
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ab.stratify) +
    theme_bw()+
    labs(title = paste0("no filter M-A boxplot true vs non-true scale",scalesize,"seqx",seqdepth," abs(M)=",abslt))
  #ggsave("/p/keles/fandingzhou/volumeA/DCI/plots/Rplots.png", width = 14, height = 8, dpi=300, compression = "lzw")
  print(p1)
  
  #filtering
  dist.keep <-  pairdist(data)>2*bin_size
  data<-data[dist.keep,]
  truesig<-truesig[dist.keep]
  ave.ab <- aveLogCPM(asDGEList(data))
  count.keep <- ave.ab>= aveLogCPM(2,lib.size=mean(data$totals))
  data<-data[count.keep,]
  truesig<-truesig[count.keep]
  
  # Calculate A
  ab<-aveLogCPM(asDGEList(data))
  ab.stratify<-cut(ab,c(-4,-2,0,2,4,7))
  # Calculate counts per million
  adj.counts<-cpm(asDGEList(data),log=TRUE)
  # Calculate M
  if(abslt){
    mval<- abs(adj.counts[,6]-adj.counts[,2])
  }else{
    mval<- adj.counts[,6]-adj.counts[,2]
  }
  #plot(ab, mval,xlab="A",ylab="M",main="Spike-in GM2 vs. Original GM2",pch=16)
  true_vs_nontrue<-data.frame(ab.stratify,mval,truesig)
  p2<-ggplot(true_vs_nontrue, aes(truesig, mval,fill=ab.stratify)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ab.stratify) +
    theme_bw()+
    labs(title = paste0("M-A boxplot true vs non-true scale",scalesize,"seqx",seqdepth," abs(M)=",abslt))
  print(p2)
}

diffHic_MA_hist<-function(data,nrep,bin_size,scalesize,TSetKey,seqdepth){
  useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
  inter.frame <- as.data.frame(interactions(data))[,useful.cols]
  inter_pairs<-paste0(inter.frame[,6]-bin_size/2,"_",inter.frame[,3]-bin_size/2)
  truesig<-rep(0,length(inter_pairs))
  truesig[inter_pairs %in% TSetKey]<-1
  truesig<-as.factor(truesig)
  
  ab<-aveLogCPM(asDGEList(data))
  # Order avg logCPM
  o<-order(ab)
  # Calculate counts per million
  adj.counts<-cpm(asDGEList(data),log=TRUE)
  # Calculate M
  mval<- adj.counts[,6]-adj.counts[,2]
  ab.stratify<-cut(ab,c(-4,-2,0,2,4,7))
  ab.catergory<-unique(ab.stratify)
  for(i in 1:length(ab.catergory)){
    hist(mval[truesig==1&ab.stratify==ab.catergory[i]],main = paste0("no filter M-A histogram truesig scale",scalesize,"seqx",seqdepth," A in ",ab.catergory[i]))
    hist(mval[truesig==0&ab.stratify==ab.catergory[i]],main = paste0("no filter M-A histogram truesig scale",scalesize,"seqx",seqdepth," A in ",ab.catergory[i]))
  }
  
  
  dist.keep <-  pairdist(data)>2*bin_size
  data<-data[dist.keep,]
  truesig<-truesig[dist.keep]
  ave.ab <- aveLogCPM(asDGEList(data))
  count.keep <- ave.ab>= aveLogCPM(2,lib.size=mean(data$totals))
  data<-data[count.keep,]
  truesig<-truesig[count.keep]
  
  # Calculate A
  ab<-aveLogCPM(asDGEList(data))
  ab.stratify<-cut(ab,c(-4,-2,0,2,4,7))
  # Calculate counts per million
  adj.counts<-cpm(asDGEList(data),log=TRUE)
  # Calculate M
  mval<- adj.counts[,6]-adj.counts[,2]
  ab.catergory<-unique(ab.stratify)
  for(i in 1:length(ab.catergory)){
    hist(mval[truesig==1&ab.stratify==ab.catergory[i]],main = paste0("M histogram truesig scale",scalesize,"seqx",seqdepth," A in ",ab.catergory[i]))
    hist(mval[truesig==0&ab.stratify==ab.catergory[i]],main = paste0("M histogram truesig scale",scalesize,"seqx",seqdepth," A in ",ab.catergory[i]))
  }
  #plot(ab, mval,xlab="A",ylab="M",main="Spike-in GM2 vs. Original GM2",pch=16)
}
