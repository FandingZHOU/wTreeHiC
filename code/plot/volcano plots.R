options("scipen"=100)

dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
in_path=c("/p/keles/fandingzhou/volumeA/DCI/ABcompartment")
in_path1=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data/originaltable","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/originaltable/filter")
out_path=c("/p/keles/fandingzhou/volumeA/DCI/multiHiCcompare/data","/p/keles/fandingzhou/volumeA/DCI/diffHic/data/filter")
func_path="/p/keles/fandingzhou/volumeA/DCI/code/weight/"
meth=c("multiHiCcompare","diffHic")
seqdepth<-c(1,3,5)
group<-c("A","GM12","GM13","GM14")
doc_name<-c("_Multi_glm","_diffHic.tsv")
group_size<-c("regular","large","small")
scalesize<-c(2,3,5,10)
bin_size=40000


library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(dplyr,lib.loc = dict_path)# 

library(ggplot2,lib.loc = dict_path)

TSetKey<-unique(scan(paste0(out_path[2],"/TSetKey"),what = character()))


diffHic_volvano<-function(HiCbinpairs,TSetKey,scalesize,seqdepth){
allbinpair<-paste0(HiCbinpairs$end2-bin_size/2,"_",HiCbinpairs$end1-bin_size/2)
volcano<-subset(HiCbinpairs,select = c(FDR,logFC))
threshold<-as.factor(allbinpair %in% TSetKey)
r03=ggplot(volcano,aes(logFC,-log(FDR,10),colour=threshold))+geom_point()
r04=r03+labs(title=paste0("scale",scalesize,"seqx",seqdepth))+theme(plot.title = element_text(hjust = 0.5))+xlim(-10,10)
r05=r04+geom_vline(xintercept=c(-1.5,1.5),linetype="dotted",size=1)+geom_hline(yintercept=-log(0.05,10),col="blue") +
  xlab("log2 fold change") + ylab("-log10 adj p-value") 
print(r05)
}


MultiHiCcompare_volvano<-function(HiCbinpairs,TSetKey,scalesize,seqdepth){
  allbinpair<-paste0(HiCbinpairs$region1+bin_size/2,"_",HiCbinpairs$region2+bin_size/2)
  volcano<-subset(HiCbinpairs,select = c(p.value,logFC))
  threshold<-as.factor(allbinpair %in% TSetKey)
  r03=ggplot(volcano,aes(logFC,-log(p.adjust(p.value,"BH"),10),colour=threshold))+geom_point()
  r04=r03+labs(title=paste0("scale",scalesize,"seqx",seqdepth))+theme(plot.title = element_text(hjust = 0.5))+xlim(-10,10)
  r05=r04+geom_vline(xintercept=c(-1.5,1.5),linetype="dotted",size=1)+geom_hline(yintercept=-log(0.05,10),col="blue") +
    xlab("log2 fold change") + ylab("-log10 adj p-value") 
  print(r05)
}

m=1
n=1
HiCbinpairs<-read_tsv(paste0(in_path1[m],"/binpairs_GMvsGMspikein_seqx",seqdepth[n],doc_name[m]),col_names=TRUE)
#diffHic_volvano(HiCbinpairs,TSetKey,1,seqdepth[n])
MultiHiCcompare_volvano(HiCbinpairs,TSetKey,1,seqdepth[n])
for(ss in 1:4){
  #for (n in 1:3){
    HiCbinpairs<-read_tsv(paste0(in_path1[m],"/binpairs_GMvsGMspikein_scale",scalesize[ss],"_seqx",seqdepth[n],doc_name[m]),col_names=TRUE)
    #smoothScatter(-log(HiCbinpairs$PValue,10),abs(HiCbinpairs$logFC),xlab="-log10pvalue",ylab="|logFC|",main=paste0(meth[m],"_scale",scalesize[ss],"_seqx",seqdepth[n]))
    #plot(-log(HiCbinpairs$PValue,10),HiCbinpairs$logFC)
    #diffHic_volvano(HiCbinpairs,TSetKey,scalesize[ss],seqdepth[n])
    MultiHiCcompare_volvano(HiCbinpairs,TSetKey,scalesize[ss],seqdepth[n])
  #}
}
