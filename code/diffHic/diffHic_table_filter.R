diffHic_table_filter<-function(data,binpairs,nrep,bin_size){
  useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
  inter.frame <- as.data.frame(interactions(data))[,useful.cols]
  inter_pairs<-paste0(inter.frame[,6]-bin_size/2,"_",inter.frame[,3]-bin_size/2)
  truesig<-which(inter_pairs %in% binpairs)
  #filtering
  data<-data[truesig,]
  
  #normalization
  data <- normOffsets(data, method="loess", se.out=TRUE)
  
  #DCI
  design <- model.matrix(~factor(c(rep("Group1",nrep),rep("Group2",nrep))))
  colnames(design) <- c("Intercept", "Group2")
  y <- asDGEList(data)
  # Estimate the dispersion
  y <- estimateDisp(y, design)
  # Fit GLM
  fit <- glmQLFit(y, design,robust=TRUE)
  # Perform F test
  result <- glmQLFTest(fit,coef=2)
  # BH method
  adj.p <- p.adjust(result$table$PValue, method="BH")
  inter.frame <- as.data.frame(interactions(data))[,useful.cols]
  results.r <- data.frame(inter.frame, result$table, FDR=adj.p)
  o.r <- order(results.r$PValue)
  HiCbinpairs_data<-results.r[o.r,]
  #HiCbinpairs_data<-cbind(HiCbinpairs$end2-bin_size/2,HiCbinpairs$end1-bin_size/2,HiCbinpairs$F,HiCbinpairs$PValue)
  #colnames(HiCbinpairs_data)<-c("bin_1","bin_2","F","P")
  return(HiCbinpairs_data)
}