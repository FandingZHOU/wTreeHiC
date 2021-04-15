#change a function to get f-value
.glm_reformat <- function(result, hic_table, p.method) {
  # create table of location info and p-value results
  result <- cbind(hic_table, result$table)
  colnames(result)[ncol(result)] <-"p.value"
  # adjust p-values
  result$p.adj <- p.adjust(result$p.value, method = p.method)
  # convert logFC from natural log to log2
  result$logFC <- log2(exp(result$logFC))
  return(result)
}

multiHiCcompare_table<-function(sample_list,nrep,Amin,reformat=F){
  GM12878_GMSPIKEIN <- make_hicexp(data_list = sample_list,A.min=Amin,zero.p=1,groups= c(rep(1,nrep),rep(2,nrep))) #A.min=5,
  
  ##Nomalization
  GM12878_GMSPIKEIN <- fastlo(GM12878_GMSPIKEIN)
  
  #GLM test
  d <- model.matrix(~factor(meta(GM12878_GMSPIKEIN)$group))
  GM12878_GMSPIKEIN_glm <- hic_glm(GM12878_GMSPIKEIN,design=d,coef=2)
  HiCbinpairs<-tryCatch(results(GM12878_GMSPIKEIN_glm),error=function(e){GM12878_GMSPIKEIN_glm@comparison})
  HiCbinpairs<-data.frame(HiCbinpairs)
  if(!reformat){
    return(HiCbinpairs)
  }else{
    HiCbinpairs_data<-data.frame(bin_1=HiCbinpairs$region1+bin_size/2,bin_2=HiCbinpairs$region2+bin_size/2,F=HiCbinpairs$F,P=HiCbinpairs$p.value)
    return(HiCbinpairs_data)
  }
}

multiHiCcompare_filter_table<-function(sample_list,filter_pairs,bin_size,nrep,reformat=F){
  GM12878_GMSPIKEIN <- make_hicexp(data_list = sample_list,groups= c(rep(1,nrep),rep(2,nrep)),filter=F) #A.min=5,
  
  #filtering
  hictable<-hic_table(GM12878_GMSPIKEIN)
  binpairs<-paste0(hictable$region1+bin_size/2,"_",hictable$region2+bin_size/2)
  GM12878_GMSPIKEIN@hic_table<-hictable[binpairs %in% filter_pairs,]
  ##Nomalization
  GM12878_GMSPIKEIN <- fastlo(GM12878_GMSPIKEIN)
  
  #GLM test
  d <- model.matrix(~factor(meta(GM12878_GMSPIKEIN)$group))
  GM12878_GMSPIKEIN_glm <- hic_glm(GM12878_GMSPIKEIN,design=d,coef=2)
  HiCbinpairs<-tryCatch(results(GM12878_GMSPIKEIN_glm),error=function(e){GM12878_GMSPIKEIN_glm@comparison})
  HiCbinpairs<-data.frame(HiCbinpairs)
  if(!reformat){
    return(HiCbinpairs)
  }else{
    HiCbinpairs_data<-data.frame(bin_1=HiCbinpairs$region1+bin_size/2,bin_2=HiCbinpairs$region2+bin_size/2,F=HiCbinpairs$F,P=HiCbinpairs$p.value)
    return(HiCbinpairs_data)
  }
  
}


use_multiHiCcompare<-function(input_path,out_path,res,nrep,Amin=NULL,other_filtering=F,filter_file=NULL,reformat=T,long_range=T,TSetKey=NULL){
  options("scipen"=100)
  if(!dir.exists(paste0(out_path,"/multiHiCcompare/originaltable/"))){
    dir.create(paste0(out_path,"/multiHiCcompare/originaltable/"))
  }
  if(!dir.exists(paste0(out_path,"/Trueset/"))){
    dir.create(paste0(out_path,"/Trueset/"))
  }
  sample_list <- list()
  for(j in 1:(nrep*2)){
    data_j<-read_tsv(input_path[j])
    #filter out short-range interactions
    if(long_range){
      original_list<-cbind(rep(1,nrow(data_j)),data_j[,2]-res/2,data_j[,4]-res/2,data_j[,5])
      sample_list[[j]]<-original_list[abs(original_list[,2]-original_list[,3])>res*2,]
    }else{
      sample_list[[j]]<-original_list
    }
    colnames(sample_list[[j]])<-c('chr','region1','region2','IF')
  }
  
  if(!other_filtering){
    HiCbinpairs<-multiHiCcompare_table(sample_list,nrep,Amin,reformat)
    write_tsv(HiCbinpairs,path = paste0(out_path,"/multiHiCcompare/originaltable/binpairs_Multi_glm"))
    bin_size=res
    if(reformat){
      binpairs<-paste0(HiCbinpairs$bin_1,"_",HiCbinpairs$bin_2) 
    }else{
      binpairs<-paste0(HiCbinpairs$region1+bin_size/2,"_",HiCbinpairs$region2+bin_size/2)
    }
    write(binpairs, file=paste0(out_path,"/Trueset/pairs_keep"))
    if(!is.null(TSetKey)){
      nosmooth<-intersect(TSetKey,binpairs)
      write(nosmooth, file=paste0(out_path,"/Trueset/pairs_TSet"))
    }
  }else{
    filter_pairs<-scan(filter_file,what = character())
    HiCbinpairs<-multiHiCcompare_filter_table(sample_list,filter_pairs,bin_size=res,nrep,reformat)
    write_tsv(HiCbinpairs,path = paste0(out_path,"/multiHiCcompare/originaltable/binpairs_Multi_glm"))
  }
  return(HiCbinpairs)
}
                              
                                
                              
                              


