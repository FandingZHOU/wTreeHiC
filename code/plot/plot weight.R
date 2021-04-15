options("scipen"=100)

meth=c("multiHiCcompare","diffHic","HiCDCPlus")

dict_path="/p/keles/fandingzhou/volumeA/software/rpackages"
bin_size<-40000
weight_method<-c("xi","pi0")
#install.packages('utf8',lib = dict_path)#reshape2,fansi'
library('utf8',lib.loc = dict_path)
library(cli,lib.loc = dict_path)
library(fansi,lib.loc = dict_path)
library(crayon,lib.loc = dict_path)
library(readr,lib.loc = dict_path)
library(data.table,lib.loc = dict_path)# 
library(ggplot2)
library(cowplot,lib.loc = dict_path)

seed=4
n=2
k=1
par(mfrow=c(2,2))
for(test in c("2vs2","4vs4")){
  for(m in 1:3){
    for(foldchange in c(2,3,4,6)){
      in_path4=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/",test,"/filter_as_",meth[n])
      HiCbinpairs_data<-read_tsv(paste0(in_path4,"/",meth[m],"/weighttable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],".tsv"))
      hist(HiCbinpairs_data$P,xlab = "pvalue",main=paste0(test,"_",meth[m],"_",foldchange))
    }
  }
}


#for(n in 1:2){
n=2
seed=4
in_path4=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/4vs4/filter_as_",meth[n])
for(test in c("2vs2","4vs4")){
  for(m in 1:3){
    for (foldchange in c(3,4,6)){
      Truesig=scan(paste0(in_path4,"/Trueset/GMvsGMFC",foldchange,"_TSet"),what = character())
      pic<-list()
      for(k in 1:2){
        in_path2=paste0("/p/keles/fandingzhou/volumeA/DCI/filter/seed",seed,"/",test,"/filter_as_",meth[n],"/",meth[m])
        HiCbinpairs_data<-read_tsv(paste0(in_path2,"/weighttable/",meth[m],"_binpairs_GMvsGMFC",foldchange,"_weightedp_",weight_method[k],".tsv"))
        binpair<-paste0(HiCbinpairs_data$bin_1,"_",HiCbinpairs_data$bin_2)
        truesig<-binpair %in% Truesig
        wpoverp<-data.frame(weight=HiCbinpairs_data$weight,truesig=as.factor(truesig))
        #wpoverp<-wpoverp[wpoverp$weight<20,]

        pic[[k]]<-ggplot(wpoverp, aes(weight, fill = truesig)) + 
               geom_histogram(alpha = 0.5, aes(), position = 'identity') +
          xlab(paste0("original p/weighted p   ",test,"_",meth[m],"_FC",foldchange))+xlim(-0.5,5)
          #labs(title = paste0("filter_as_",meth[n],"_",meth[m],"_FC",foldchange))
        
        # true_vs_nontrue<-data.frame(ab.stratify,wpoverp,truesig)
        # p1<-ggplot(true_vs_nontrue, aes(truesig, mval,fill=ab.stratify)) + #
        #   geom_boxplot(outlier.shape = NA) +
        #   facet_wrap(~ab.stratify) +
        #   theme_bw()+
        #   labs(title = paste0("no filter M-A boxplot true vs non-true scale",scalesize,"seqx",seqdepth," abs(M)=",abslt))
        # ggsave("/p/keles/fandingzhou/volumeA/DCI/plots/Rplots.png", width = 14, height = 8, dpi=300, compression = "lzw")
        
      }  
      print(plot_grid(pic[[1]], pic[[2]],
                labels = c("w1", "w2"),
                ncol = 1, nrow = 2))
    }
  }
}
#}