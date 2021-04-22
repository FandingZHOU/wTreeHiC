binA<-function(x,bin_size,tomid=F){
  if(tomid==F){
    return(as.numeric(x[1])-bin_size/2)
  }else{
    return(as.numeric(x[1])+bin_size/2)
  }
}
binB<-function(x,bin_size,tomid=F){
  if(tomid==F){
    return(as.numeric(x[2])-bin_size/2)
  }else{
    return(as.numeric(x[2])+bin_size/2)
  }
}

hicmatrix<-function(inputfile,out_path,n_bin=NULL,bin_size,ch1_length=NULL,chr){
  out_path=paste0(out_path,"/hicmatrix")
  if(!dir.exists(out_path)){
    dir.create(out_path)
  }
  
  if(is.null(ch1_length)){
    all.regions<-data.frame(cbind(rep(paste0("chr",chr),n_bin+1),seq(1,n_bin*bin_size+1,bin_size),c(seq(bin_size,n_bin*bin_size,bin_size)),1:(n_bin+1)))
  }else{
    n_bin=floor(ch1_length/bin_size)
    all.regions<-data.frame(cbind(rep(paste0("chr",chr),n_bin+1),seq(1,n_bin*bin_size+1,bin_size),c(seq(bin_size,n_bin*bin_size,bin_size),ch1_length),1:(n_bin+1)))
  }
  data.table::fwrite(all.regions,file=paste0(out_path,"/bedfile"),sep = "\t",col.names = F)
  out_path=out_path3
  nrep<-length(inputfile)/2
  for (j in 1:(nrep*2)){
    binpair_data<-read_tsv(inputfile[j])
    binpair_data<-binpair_data[,c(2,4,5)]
    binpair_data[,1]<-(binpair_data[,1]+bin_size/2)/bin_size
    binpair_data[,2]<-(binpair_data[,2]+bin_size/2)/bin_size
    if(j<=nrep){
      data.table::fwrite(binpair_data,paste0(out_path,"/G1_rep",j),sep = "\t",col.names = F)
      
    }else{
      data.table::fwrite(binpair_data,paste0(out_path,"/G2_rep",j-nrep),sep = "\t",col.names = F)
      
    }
  }
  return(nrep)
}
##inpath contains both bedfile and IF matrix files
hicdcmatrix<-function(in_path,out_path,nrep,chr,n_bin,bin_size,ch1_length=NULL,ncore=2,long_range=T,gen = "Mmusculus",
                      gen_ver = "mm9"){
  out_path=paste0(out_path,"/hicmatrix")
  if(!dir.exists(out_path)){
    dir.create(out_path)
  }
  
  indexfile<-data.frame()
  if(is.null(ch1_length)){
    ch1_length=bin_size*n_bin
  }
  for(j in 1:(nrep*2)){
    if(j<=nrep){
      Gp<-"G1"
      repl=j
    }else{
      Gp<-"G2"
      repl=j-nrep
    }
    gi_list<-generate_binned_gi_list(bin_size,chrs=paste0('chr',chr),Dthreshold = ch1_length,gen,
                                     gen_ver )
    # gi_list<-generate_binned_gi_list(bin_size,chrs=paste0('chr',chr),Dthreshold = ch1_length, gen = "Mmusculus",
    #                                  gen_ver = "mm9")
    gi_list<-add_hicpro_matrix_counts(
      gi_list,
      paste0(in_path,"/bedfile"),
      paste0(in_path,"/",Gp,"_rep",repl),
      chrs=rep("chr",chr)
    )
    gi_list<-expand_1D_features(gi_list)
    #run HiC-DC+ on 2 cores
    set.seed(1010) #HiC-DC downsamples rows for modeling
    if(long_range==T){
      gi_list<-HiCDCPlus(gi_list,Dmin = 3*bin_size,Dmax = ch1_length)
    }else{
      gi_list<-HiCDCPlus(gi_list,Dmin = 0,Dmax = ch1_length)
    }
    gi_list_qval<-gi_list[[1]]$qvalue
    gi_list[[1]]<-gi_list[[1]][(!is.na(gi_list_qval))]
    indexfile<-unique(rbind(indexfile,
                            as.data.frame(gi_list[[1]][gi_list[[1]]$qvalue<=0.05])[c('seqnames1','start1','start2')]))
    saveRDS(gi_list,paste0(out_path,"/",Gp,"_rep",repl,"sig_interact.rds"))
  }
  colnames(indexfile)<-c('chr','startI','startJ')
  data.table::fwrite(indexfile,
                     paste0(out_path,'/filter_indices.txt.gz'),
                     sep='\t',row.names=FALSE,quote=FALSE)
  
}
hicdcdiff2<-function(input_paths,output_path,filter_path,bin_size,chr,ch_length,long_range=T,nrep,reformat=T,...){
     filter1<-scan(filter_path,what = character())
     filterset<-strsplit(x=filter1, split="_")
      indexfile<- data.frame()
      indexfile<-cbind(rep(paste0('chr',chr),length(filter1)),sapply(filterset, binA,bin_size),sapply(filterset, binB,bin_size))
      colnames(indexfile)<-c('chr','startI','startJ')
      data.table::fwrite(indexfile,
                        paste0(output_path,'/filter_indices.txt.gz'),
                        sep='\t',row.names=FALSE,quote=FALSE)
      output_path2=paste0(output_path,"/HiCDCPlus/originaltable/")
      if(!dir.exists(output_path2)){
          dir.create(output_path2)
        }
      filter_file=paste0(output_path,'/filter_indices.txt.gz')
      if(long_range){
          Dmin=3*bin_size
        }else{
            Dmin=0
        }
      
      input_files<-list(Gp1=paste0(input_paths,"/G1_rep",1:nrep,"sig_interact.rds"),
                        Gp2=paste0(input_paths,"/G2_rep",1:nrep,"sig_interact.rds"))
      hicdcdiff(input_paths=input_files,
                              output_path=output_path2,
                              filter_file=filter_file,
                              Dmin =Dmin,
                              Dmax = ch_length,
                              binsize=bin_size,
                              granularity=bin_size,
                              DESeq.save=TRUE,
                              diagnostics=F,
                              reformat = reformat
                    )
   
      }
.readDiffInputFiles <- function(conds, input_paths, sigs, chrom, Dmin, Dmax, dband, bin_type, binsize) {
  retlist <- list()
  normfac <- NULL
  for (cond in conds) {
    numrep <- length(input_paths[[cond]])
    for (i in seq(1, numrep)) {
      prefix <- input_paths[[cond]][i]
      if (!grepl("\\.hic$", prefix, ignore.case = TRUE)){
        if (is.list(prefix)&methods::is(prefix[[1]], "GInteractions")){
          gi_list_validate(prefix)
          normfac_add <- prefix[[chrom]]
          rm(prefix)
        } else {
          gi_list <- gi_list_read(path.expand(prefix))
          gi_list_validate(gi_list)
          normfac_add <- gi_list[[chrom]]
          rm(gi_list)
        }
        normfac_add <- data.frame(chr = chrom, startI = BiocGenerics::start(InteractionSet::anchors(normfac_add)$first), 
                                  startJ = BiocGenerics::start(InteractionSet::anchors(normfac_add)$second), counts = mcols(normfac_add)$counts, 
                                  D = mcols(normfac_add)$D, stringsAsFactors = FALSE)
        normfac_add <- normfac_add %>% dplyr::filter(.data$counts > 0)
      } else if (grepl("\\.hic$", prefix, ignore.case = TRUE)) {
        if (bin_type == "Bins-uniform") {
          if (!.Platform$OS.type=="windows"){
            normfac_add <- tryCatch(
              straw(norm = "NONE", fn = path.expand(prefix), bs = binsize, ch1 = gsub("chr", "", chrom), ch2 = gsub("chr", "", chrom), 
                    u = "BP"),
              error=function(e){
                tryCatch(straw(norm = "NONE", fn = path.expand(prefix), bs = binsize, ch1 = chrom, ch2 = chrom, 
                               u = "BP"),
                         error=function(e){
                           straw_dump(norm = "NONE",fn=path.expand(prefix),bs=binsize,ch1=gsub("chr", "", chrom),ch2=gsub("chr", "", chrom),u="BP")   
                         })
              })
          }else{
            normfac_add<-straw_dump(norm = "NONE",fn=path.expand(prefix),bs=binsize,ch1=gsub("chr", "", chrom),ch2=gsub("chr", "", chrom),u="BP")   
            gc(reset=TRUE,full=TRUE)
          }
        } else {
          if (!.Platform$OS.type=="windows"){
            normfac_add <- tryCatch(
              straw(norm = "NONE", fn = path.expand(prefix), bs = binsize, ch1 = gsub("chr", "", chrom), ch2 = gsub("chr", "", chrom), 
                    u = "FRAG"),
              error=function(e){
                tryCatch(straw(norm = "NONE", fn = path.expand(prefix), bs = binsize, ch1 = chrom, ch2 = chrom, 
                               u = "FRAG"),
                         error=function(e){
                           straw_dump(norm = "NONE",fn=path.expand(prefix),bs=binsize,ch1=gsub("chr", "", chrom),ch2=gsub("chr", "", chrom),u="FRAG")   
                         })
              })
          }else{
            normfac_add<-straw_dump(norm = "NONE",fn=path.expand(prefix),bs=binsize,ch1=gsub("chr", "", chrom),ch2=gsub("chr", "", chrom),u="FRAG")   
            gc(reset=TRUE,full=TRUE)
          }
        }
        colnames(normfac_add) <- c("startI", "startJ", "counts")
        normfac_add <- normfac_add %>% dplyr::mutate(chr = chrom, D = abs(.data$startI - .data$startJ))
      } else {
        stop(paste0("File not found relating to ", prefix))
      }
      normfac_add <- normfac_add %>% dplyr::filter(.data$D >= Dmin & .data$D <= Dmax) %>% dplyr::mutate(Dband = as.numeric(cut(.data$D, 
                                                                                                                               breaks = dband, labels = seq(1, (length(dband) - 1), 1), include.lowest = TRUE))) %>% dplyr::select(.data$Dband, 
                                                                                                                                                                                                                                   .data$chr, .data$startI, .data$startJ, .data$counts)
      # rename counts column
      colnames(normfac_add)[colnames(normfac_add) %in% "counts"] <- paste0(cond, ".", i)
      if (is.null(normfac)) {
        normfac <- normfac_add
      } else {
        # join counts across conditions and replicates
        normfac <- dplyr::full_join(normfac, normfac_add)
        normfac<-dplyr::mutate_all(normfac, ~replace(.,is.na(.),0))
      }
    }
  }
  # calculate geometric mean of all count columns
  countcols <- colnames(normfac)[!colnames(normfac) %in% c("chr", "startI", "startJ", "Dband")]
  normfac$geomean <- exp(rowMeans(log(normfac[countcols] + 0.5)))
  # calculate factors: count/geometric mean
  for (countcol in countcols) {
    normfac[, paste0(countcol, ".fac")] <- (normfac[, countcol] + 0.5)/normfac$geomean
  }
  # calculate median factor by distance band and store it as normfac.final
  countcols <- paste0(countcols, ".fac")
  normfac.final <- normfac %>% dplyr::group_by(.data$Dband) %>% dplyr::summarise_at(.vars = countcols, .funs = stats::median)
  normfac.final$geomean <- exp(rowMeans(log(normfac.final[countcols])))
  # join back with normfac with normfac.final columns .fac replaced with .norm and geomean removed
  colnames(normfac.final) <- gsub(".fac", ".norm", colnames(normfac.final))
  normfac <- dplyr::left_join(normfac, normfac.final %>% dplyr::select(-.data$geomean))
  # significant bins to be filtered
  sigbins <- sigs %>% dplyr::filter(.data$chr == chrom)
  sigbins <- paste0(sigbins$startI, ":", sigbins$startJ)
  retlist[["normfac"]] <- normfac
  retlist[["normfac.final"]] <- normfac.final
  retlist[["sigbins"]] <- sigbins
  return(retlist)
}
# 
# bin_type = "Bins-uniform"
# binsize = bin_size
# granularity = bin_size
# chrs = NULL
# Dmax = ch_length
# diagnostics = FALSE
# DESeq.save = T
# fitType = "local"
hicdcdiff <- function(input_paths, filter_file, output_path, bin_type = "Bins-uniform", binsize = 5000, granularity = 5000, 
                      chrs = NULL, Dmin = 0, Dmax = 2e+06, diagnostics = FALSE, DESeq.save = FALSE, fitType = "local",reformat=F) {
  options(scipen = 100, digits = 4, warn = -1)
  conds <- names(input_paths)
  granularity <- ifelse(bin_type == "Bins-uniform" & granularity == 5000, binsize, granularity)
  dband <- seq(0, Dmax, granularity)
  dband[1] <- Dmin
  dband <- unique(sort(dband[dband >= Dmin]))
  deseq2paths <- NULL
  outputpaths <- NULL
  plotpaths <- NULL
  sigs <- data.table::fread(filter_file)
  #set memory limit to max if i386
  if (.Platform$OS.type=='windows'&Sys.getenv("R_ARCH")=="/i386") {
    gc(reset=TRUE,full=TRUE)
    utils::memory.limit(size=4095)
  }
  # set default chrs if need be
  if (is.null(chrs)) {
    # get list of chromosomes from sigs
    if (!"chr" %in% colnames(sigs)) 
      stop("No column named 'chr' in filter_file.")
    chrs <- sort(unique(sigs$chr))
  }
  
  for (chr in chrs) {
    plotpaths_add <- NULL
    outputpaths_add <- NULL
    deseq2paths_add <- NULL
    retlist <- .readDiffInputFiles(conds, input_paths, sigs, chr, Dmin, Dmax, dband, bin_type, binsize)
    if (diagnostics) {
      plotpaths_add <- .plotNormalizationFactors(retlist[["normfac.final"]], binsize, chr, dband, conds, output_path)
    }
    countcols <- c()
    for (cond in conds) {
      countcols <- c(countcols, paste0(cond, ".", seq(length(input_paths[[cond]]))))
    }
    count.mat <- retlist[["normfac"]][countcols]
    rownames(count.mat) <- paste0(retlist[["normfac"]]$startI, ":", retlist[["normfac"]]$startJ)
    normcols <- paste0(countcols, ".norm")
    normfac.mat <- retlist[["normfac"]][normcols]
    rownames(normfac.mat) <- rownames(count.mat)
    sampleinfo <- data.frame(condition = rep(conds, times = unlist(lapply(input_paths, length))))
    rownames(sampleinfo) <- colnames(count.mat)
    sampleinfo$samples <- rownames(sampleinfo)
    count.mat <- as.matrix(count.mat)
    normfac.mat <- as.matrix(normfac.mat)
    ## run deseq2
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.mat, colData = sampleinfo, design = ~condition)
    dds$condition <- stats::relevel(dds$condition, ref = conds[1])
    DESeq2::normalizationFactors(dds) <- normfac.mat
    keep <- rownames(DESeq2::counts(dds)) %in% retlist[["sigbins"]]
    dds <- dds[keep, ]
    dds <- DESeq2::estimateDispersionsGeneEst(dds)
    DESeq2::dispersions(dds) <- tryCatch({
      dds <- DESeq2::estimateDispersionsFit(dds, fitType = fitType)
      GenomicRanges::mcols(dds)$dispFit},
      error=function(e) {
        msg<-paste0("fitType failed to fit for ",chr,". Overriding with gene estimates")
        message(msg)
        return(mcols(dds)$dispGeneEst)
      })
    dds <- DESeq2::nbinomWaldTest(dds)
    
    # save DESeq file
    if (DESeq.save) {
      deseq2path <- paste0(output_path,"/", chr, "_DESeq2_obj.rds")
      deseq2output <- path.expand(deseq2path)
      deseq2outputdir<-gsub("/[^/]+$", "",deseq2output)
      if (deseq2outputdir==deseq2output){
        deseq2outputdir<-gsub("\\[^\\]+$", "",deseq2output)
      }
      if (deseq2outputdir==deseq2output){
        deseq2outputdir<-gsub("\\\\[^\\\\]+$", "",deseq2output)
      }
      if (!deseq2outputdir==deseq2output&!dir.exists(deseq2outputdir)){
        dir.create(deseq2outputdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
      }
      if(!reformat){
        saveRDS(dds, deseq2path)
      }else{
        HiCbinpairs<-results(dds)
        tset<-strsplit(x=rownames(HiCbinpairs), split=":")
        HiCbinpairs_data<-data.frame(cbind(sapply(tset, binA,bin_size,T),sapply(tset, binB,bin_size,T),HiCbinpairs$stat,HiCbinpairs$pvalue))
        colnames(HiCbinpairs_data)<-c("bin_1","bin_2","Z","P")
        deseq2path2<-paste0(output_path, "/",chr, "_DESeq2_table")
        write_tsv(HiCbinpairs_data, deseq2path2)
      }
      
      deseq2paths_add <- c(deseq2paths_add, deseq2path)
    }
    deseq2paths <- c(deseq2paths, deseq2paths_add)
    outputpaths <- c(outputpaths, outputpaths_add)
    plotpaths <- c(plotpaths, plotpaths_add)
  }
  return(list(deseq2paths = deseq2paths, outputpaths = outputpaths, plotpaths = plotpaths))
}