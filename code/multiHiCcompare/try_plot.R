threshold<-1:4
detected_fdr<-c(rep(1,4),rep(2,4))
theoretical_fdr<-c(c(1:4),c(1:4))
supp=c(rep("orignal",length(threshold)),rep("weighted",length(threshold)))
supp1=rep(supp,4)
detected_fdr1=rep(detected_fdr,4)
theoretical_fdr1=rep(theoretical_fdr,4)
cate<-c(rep("multi",16),rep("diff",16))
seq_d<-rep(c(rep("seqx1",8),rep("seqx2",8)),2)
tgg=data.frame(supp1,theoretical_fdr1,detected_fdr1,cate,seq_d)

ggplot(tgg, aes(x=theoretical_fdr1, y=detected_fdr1, color=supp1,shape=supp1)) + geom_line(size=1) +geom_point(size=2)+facet_grid(seq_d~cate)
length(cate)
rbind(data.frame(),tgg)
