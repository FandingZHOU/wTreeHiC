diffHic_table<-function(data,nrep,bin_size,filter_par,reformat=F,long_range=T){
  #filtering
  if(long_range){
    dist.keep <-  pairdist(data)>2*bin_size
    data<-data[dist.keep,]
    
  }
  
  ave.ab <- aveLogCPM(asDGEList(data))
  count.keep <- ave.ab>= aveLogCPM(filter_par,lib.size=mean(data$totals))
  data<-data[count.keep,]

  # trended <- filterTrended(data)
  # trend.keep2 <- trended$abundances > trended$threshold 
  # data<-data[trend.keep2,]
  # summary(trend.keep2)
  
  
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
  useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
  inter.frame <- as.data.frame(interactions(data))[,useful.cols]
  results.r <- data.frame(inter.frame, result$table, FDR=adj.p)
  o.r <- order(results.r$PValue)
  HiCbinpairs<-results.r[o.r,]
  if(reformat){
    HiCbinpairs_data<-data.frame(bin_1=HiCbinpairs$end2-bin_size/2,bin_2=HiCbinpairs$end1-bin_size/2,F=HiCbinpairs$F,P=HiCbinpairs$PValue)
    return(HiCbinpairs_data)
  }else{
    return(HiCbinpairs)
  }
 
}


diffHic_table_filter<-function(data,binpairs,nrep,bin_size,reformat=F){
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
  HiCbinpairs<-results.r[o.r,]
  if(reformat){
    HiCbinpairs_data<-data.frame(bin_1=HiCbinpairs$end2-bin_size/2,bin_2=HiCbinpairs$end1-bin_size/2,F=HiCbinpairs$F,P=HiCbinpairs$PValue)
    return(HiCbinpairs_data)
  }else{
    return(HiCbinpairs)
  }
  
}

use_diffHic<-function(cm,out_path,bin_size,filter_par=NULL,reformat=T,other_filtering=F,filter_file=NULL,long_range=T,TSetKey=NULL,chr){
  if(!dir.exists(paste0(out_path,"/diffHic/originaltable/"))){
    dir.create(paste0(out_path,"/diffHic/originaltable/"))
  }
  if(!dir.exists(paste0(out_path,"/Trueset/"))){
    dir.create(paste0(out_path,"/Trueset/"))
  }
  n_bin<-nrow(cm[[1]])
  nrep<-length(cm)/2
  all.regions <- GRanges(paste0("chr",chr), IRanges(0:(n_bin-1)*bin_size+1, 1:n_bin*bin_size))
  row.indices <- 1:n_bin
  col.indices <- 1:n_bin
  row.regions <- all.regions[row.indices]
  col.regions <- all.regions[col.indices]
  for(j in 1:length(cm)){
    cm[[j]] <- ContactMatrix(as.matrix(cm[[j]]), row.regions, col.regions)
    
  }
  
  if(length(cm)==4){
    data <- mergeCMs(cm[[1]], cm[[2]],cm[[3]],cm[[4]])
  }
  if(length(cm)==8){
    data <- mergeCMs(cm[[1]], cm[[2]],cm[[3]],cm[[4]],cm[[5]], cm[[6]],cm[[7]],cm[[8]])
  }
  if(long_range){
    dist.keep <-  pairdist(data)>2*bin_size
    data<-data[dist.keep,]
  }
  if(!other_filtering){
    HiCbinpairs<-diffHic_table(data,nrep,bin_size,filter_par,reformat,long_range)
    write_tsv(HiCbinpairs,path = paste0(out_path,"/diffHic/originaltable/binpairs_diffHic.tsv"))
    if(reformat){
      binpairs<-paste0(HiCbinpairs$bin_1,"_",HiCbinpairs$bin_2) 
    }else{
      binpairs<-paste0(result$end2-bin_size/2,"_",result$end1-bin_size/2) 
    }
    write(binpairs, file=paste0(out_path,"/Trueset/pairs_keep"))
    if(!is.null(TSetKey)){
      nosmooth<-intersect(TSetKey,binpairs)
      write(nosmooth, file=paste0(out_path,"/Trueset/pairs_TSet"))
    }
  }else{
    filter_pairs<-scan(filter_file,what = character())
    HiCbinpairs<-diffHic_table_filter(data,filter_pairs,nrep,bin_size,reformat)
    write_tsv(HiCbinpairs,path = paste0(out_path,"/diffHic/originaltable/binpairs_diffHic.tsv"))
  }
  return(HiCbinpairs)

}