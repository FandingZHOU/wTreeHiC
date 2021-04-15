glm.reformat<-function(inputfile,method){
  if(file.exists(inputfile)){
    if(method in c("diffHic","multiHiCcompare")){
      HiCbinpairs<-read_tsv(inputfile,col_names=TRUE)
      if(method =="multiHiCcompare"){
        HiCbinpairs_data<-cbind(HiCbinpairs$region1+bin_size/2,HiCbinpairs$region2+bin_size/2,HiCbinpairs$F,HiCbinpairs$p.value,rep(0,nrow(HiCbinpairs)),rep(0,nrow(HiCbinpairs)))
      }else{
        HiCbinpairs_data<-cbind(HiCbinpairs$end2-bin_size/2,HiCbinpairs$end1-bin_size/2,HiCbinpairs$F,HiCbinpairs$PValue)
      }
      colnames(HiCbinpairs_data)<-c("bin_1","bin_2","F","P")
      
    }else{
      HiCbinpairs<-results(readRDS(inputfile))
      tset<-strsplit(x=rownames(HiCbinpairs), split=":")
      HiCbinpairs_data<-cbind(sapply(tset, binA,bin_size,T),sapply(tset, binB,bin_size,T),HiCbinpairs$stat,HiCbinpairs$pvalue)
      colnames(HiCbinpairs_data)<-c("bin_1","bin_2","Z","P")
      
    }
    return(HiCbinpairs_data)
  }
  
  
}