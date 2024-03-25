Calculate_Effects<-function(All_data,Gene_Data,SNP_Data){
  library(lsr)
  library(effects)
  Loci_pick<-All_data[,c("snps","gene")] %>% unique #, "chr","genestart","geneend","snploc"
  SNP_all<-read.delim(sep="\t", SNP_Data)
  Exp_all<-read.delim(sep="\t", Gene_Data)
  
  Gene_data<-Exp_all[ match(Loci_pick$gene, Exp_all$gene_id),] 
  SNP_data<-SNP_all[ match(Loci_pick$snps, SNP_all$snpid),] 
  stopifnot(all(colnames(Gene_data)[-1]==colnames(SNP_data)[-1] ) )
  
  N<-rowSums(!is.na(SNP_data[,-1]))
  MAF<-rowSums(SNP_data[,-1])/2/N
  
  Gene_data<-as.list(as.data.frame(t(Gene_data[,-1])))%>% setNames(Gene_data$gene_id)
  SNP_data<-as.list(as.data.frame(t(SNP_data[,-1])))%>% setNames(SNP_data$snpid)
  
  AOV<-lapply(1:length(Gene_data), function(i){    #nrow(Gene_data)
    
    # for(i in 1:length(Gene_data)){
    #  cat(i,"\t")
    data_sub<-data.frame(gene=Gene_data[[i]], SNP=factor(SNP_data[[i]]))
    model.1 <- aov( gene~SNP, data=data_sub )
    #AOV effect size
    AOV_Effect<-etaSquared( model.1 )[1,2]
    #Mean for each level
    eff=effect( term = "SNP", mod = model.1 )
    eff=as.data.frame(summary(eff)[c(3,5,7)])
    eff$SE<-eff$effect-eff$lower
    #max_Effect<-max(abs(eff$effect))*sign(eff$effect)
    #min_Effect<-min(abs(eff$effect))*sign(eff$effect)
    
    Effect<-(eff[nrow(eff),"effect"] - eff[1,"effect"])/(nrow(eff)-1)
    SE=sqrt(eff[nrow(eff),"SE"]^2 + eff[1,"SE"]^2)
    
    data.frame(gene=names(Gene_data)[i],snps=names(SNP_data)[i],AOV_Effect=AOV_Effect, Effect=Effect, SE=SE)
    
  }
  ) 
  
  AOV<-do.call(rbind.data.frame,AOV)
  AOV$N<-N
  AOV$MAF<-MAF
  All_data<-merge(All_data, AOV, by=c("snps","gene"))
  return(All_data)
}

MatrixEqtl_wrapper<-function(Gene_Data, Gene_loc,SNP_Data,SNP_loc,Out_PATH,Out_suffix="",cisDist = 2e5,useModel = modelANOVA, threshold=0.005 ){
  require(MatrixEQTL)
  
  suppressWarnings(dir.create(paste0(Out_PATH)))
  
  
  SNP_head <- read.table(sep="\t",row.names = 1, paste0(SNP_Data), nrow=1)
  Gene_head <- read.table(sep="\t",row.names = 1, paste0(Gene_Data), nrow=1)
  if(!all(colnames(SNP_head)==colnames(Gene_head))){
    stop("Genotype and Expression headers don't match!")
  }
  
  
  
  ## Load genotype data
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;
  snps$fileSkipColumns = 1;
  snps$fileSliceSize = 2000; #Suggested 2000
  snps$LoadFile(paste0(SNP_Data));
  ## Load gene expression data
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;
  gene$fileSkipColumns = 1;
  gene$fileSliceSize = 2000; #Suggested 2000
  gene$LoadFile(paste0(Gene_Data));
  
  ## Load covariates (NULL for this)
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  
  
  ## locatoin file for cis identificaiton
  snpspos = read.table(paste0(SNP_loc), header = TRUE, stringsAsFactors = FALSE)[,c("snpid","chr","pos")];
  genepos = read.table(paste0(Gene_loc), header = TRUE, stringsAsFactors = FALSE)[,c("Gene_Name","Chr","Start","End")];
  genepos<-genepos[match(rownames(gene) , genepos$Gene_Name),]
  
  
  ## Run the analysis
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = '',
    pvOutputThreshold = 0,
    useModel = useModel,
    errorCovariance = numeric(),
    verbose = F,
    output_file_name.cis  = tempfile(),
    pvOutputThreshold.cis = 1,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = F);
  
  
  
  
  me$cis$eqtls$chr<-snpspos$chr[match(me$cis$eqtls$snps, snpspos$snpid)]
  me$cis$eqtls$snploc<-snpspos$pos[match(me$cis$eqtls$snps, snpspos$snpid)]
  me$cis$eqtls$genestart<-genepos$Start[match(me$cis$eqtls$gene, genepos$Gene_Name)]
  me$cis$eqtls$geneend<-genepos$End[match(me$cis$eqtls$gene, genepos$Gene_Name)]
  
  snpspos = read.table(paste0(SNP_loc), header = TRUE, stringsAsFactors = FALSE)[,c("snpid","chr","pos","allele")]
  me$cis$eqtls$allele<-snpspos$allele[match(me$cis$eqtls$snps,snpspos$snpid )]
  
  eQTL_all<-me$cis$eqtls
  saveRDS(eQTL_all, paste0( Out_PATH,"eQTL_All_",Out_suffix,".RDS"))
  
  eQTL_all<-subset(eQTL_all, pvalue<=threshold)
  cat(nrow(eQTL_all),"of eQTL are significant!\n")
  
  Data_w_Effects<-Calculate_Effects(All_data=eQTL_all, 
                                    Gene_Data=Gene_Data,
                                    SNP_Data=SNP_Data
  )
  
  
  write.table(Data_w_Effects,sep="\t", row.names = F,  paste0( Out_PATH,"eQTL_sig_",Out_suffix,".tsv"))
  cat("Data saved as ",paste0( Out_PATH,"eQTL_sig_",Out_suffix,".tsv") )
  return(NULL)
}




QQplot<-function(eQTL_data,Out_PATH, Out_suffix=""){
  eQTL_data<- readRDS(eQTL_data)
  chisq <- qchisq(1-eQTL_data$pvalue,1)
  lambda = round(median(chisq)/qchisq(0.5,1),2)
  Expected=-log10(qunif(ppoints(nrow(eQTL_data))))  #, max=max(-log10(Corr_result$p))
  
  suppressWarnings( dir.create(paste0(Out_PATH),recursive = T))
  #pdf(paste0(PATH,sub_PATH,"Plots/correlation/",paste0("Correlation_qqplot_Snps",Covariate_name_addup,sep="_") ,".pdf"))
  jpeg(paste0(Out_PATH,"/QQplot_",Out_suffix,".jpg"),width = 1260, height = 1260, res=300)
  qqplot(Expected,-log10(eQTL_data$pvalue),pch=16,cex=1, xlab="Expected p-quantile",ylab="Observed p-quantile",main="QQ-plot", sub=paste("Median-lambda = ",lambda))
  abline(a=0,b=1,col="steelblue",lwd=2)
  dev.off()
}


Manhattan_plot<-function(eQTL_data,Out_PATH, Out_suffix="", max_pvalue=0.1, sig_pvalue=0.005){
  require(qqman)
  eQTL_data<- readRDS(eQTL_data)
  eQTL_data$eQTLid<-paste0(eQTL_data$snps,"/",eQTL_data$gene)
  suppressWarnings( dir.create(paste0(Out_PATH),recursive = T))
  par(mar=rep(0.3,4))
  jpeg(paste0(Out_PATH,"/Mannhattan_",Out_suffix,".jpg"),res=200,width=2000,height=800)
  
  manhattan(eQTL_data[eQTL_data$pvalue<max_pvalue,],ylim=c(0,ceiling(max(-log10(eQTL_data$pvalue) ))),
            chr = "chr",bp="snploc",p="pvalue",snp="eQTLid",col=palette(), annotatePval = NULL, 
            genomewideline = FALSE, suggestiveline = -log10(sig_pvalue))
  dev.off()
}


LDmatrix_catch<-function(snps, pop = "CEU", r2d = "r2", token = NULL, file = FALSE, genome_build = genome_build_sel){
  require(LDlinkR)
  tryCatch({
    LDmatrix(snps=snps,  pop = pop,   r2d = r2d,  token = token, file = file, genome_build = genome_build  )
  },
  error = function(err){
    message(err)
    LDmatrix_catch(snps=snps,  pop = pop,  r2d = r2d,   token = token,  file = file,  genome_build = genome_build)
  }
  )
}


GWAS_Olap_wLD<-function(GWAS_file,eQTL_file, LD_PATH,Extend_range=2e5,token_current,All_Pop="EUR",GWAS_pcut=0.005,
                        genome_build_sel="grch37", API_segment_length=500, rm_temp_file=T, r2_threshold=0.8){
  warning("This function is written for the purpose of code storage. 
          Running LDlink API is very very slow for more than a few loci. 
          Please calculate your own LD or contact LDlink admin for massive access of data if you need!\n
          However, this function automatically track query process so if the process stopped, 
          just restart and it will continue without re-downloading everything again")
  #install.packages("LDlinkR")
  require(LDlinkR)
  require(GenomicRanges)
  suppressWarnings( dir.create(LD_PATH, recursive = T) )
  closeAllConnections()
  # Group snps into sub-ranges for LD extraction
  hg19_seqinfo<-Seqinfo(genome="hg19")
  
  GWAS <- read.delim(sep = "\t", GWAS_file)
  GWAS<-GWAS[GWAS$pvalue<=GWAS_pcut,]
  
  GWAS_GR<-makeGRangesFromDataFrame(GWAS, keep.extra.columns = T,ignore.strand = T, seqinfo = hg19_seqinfo, seqnames.field = "CHR", start.field = "POS", end.field = "POS")
  
  eQTL_data<-read.delim(eQTL_file, sep="\t")
  eQTL_data$chr<-gsub("^chr","", eQTL_data$chr)
  eQTL_data$chr<-paste0("chr",eQTL_data$chr)
  eQTL_data$start<-eQTL_data$snploc - Extend_range
  eQTL_data$end<-eQTL_data$snploc + Extend_range
  
  eQTL_GR<-makeGRangesFromDataFrame(eQTL_data, keep.extra.columns = T,ignore.strand = T, seqinfo = hg19_seqinfo, seqnames.field = "chr", start.field = "start", end.field = "end")%>%trim()
  
  eQTL_GR_reduced<-GenomicRanges::reduce(eQTL_GR)
  eQTL_GR_reduced_df<-data.frame(eQTL_GR_reduced, range=1:length(eQTL_GR_reduced))
  
  Olap<-findOverlaps(GWAS_GR, eQTL_GR_reduced)
  stopifnot(all(!duplicated(Olap@from)))
  GWAS<-GWAS[Olap@from,]
  GWAS$Range<-eQTL_GR_reduced_df$range[Olap@to]
  #saveRDS(data_list, paste0(PATH, SUB_PATH,"GWAS_w_eQTL_matched_Range.RDS"))
  
  Olap<-findOverlaps(eQTL_GR, eQTL_GR_reduced)
  stopifnot(all(!duplicated(Olap@from)))
  stopifnot(all(Olap@from == 1:nrow(eQTL_data)))
  eQTL_data$Range<-eQTL_GR_reduced_df$range[Olap@to]
  #saveRDS(eQTL_data, paste0(PATH, SUB_PATH,"eQTL_w_GWAS_matched_Range.RDS"))
  
  SNP_all<-rbind.data.frame(GWAS[,c("snpid","CHR", "Range")] %>% setNames(c("SNPS","CHR_ID", "Range")) , eQTL_data[,c("snps","chr", "Range")] %>% setNames(c("SNPS","CHR_ID", "Range")))
  SNP_all<-unique(SNP_all)
  #saveRDS(SNP_all, paste0(PATH, SUB_PATH,"SNP_all.RDS"))
  
  
  
  All_Ranges <- unique(GWAS$Range)%>% sort
  #All_Pop<-list_pop()$super_pop_code %>% unique
  
  Done_ranges<-list.files(LD_PATH) 
  Done_ranges<-Done_ranges[grepl("SNP_LD_sub_.*\\.RDS",Done_ranges)]
  Done_ranges<-gsub("SNP_LD_sub_|\\.RDS","",Done_ranges) %>% as.numeric %>% na.omit
  Missing_Range<-All_Ranges[!All_Ranges %in% Done_ranges]
  
  
  for(Range_i in Missing_Range){ #Range_i=1
    
    GWAS_sub<-GWAS[GWAS$Range%in%Range_i,"snpid"] %>% unique
    eQTL_sub<-eQTL_data[eQTL_data$Range%in%Range_i,"snps"] %>% unique
    GWAS_sub<-GWAS_sub[!GWAS_sub%in%eQTL_sub] 
    
    
    if((length(GWAS_sub)+length(eQTL_sub))<=API_segment_length ){
      # zz <- file(paste0(PATH, SUB_PATH, "Errorlog/errortag_",Range_i,".Rout"), open = "wt")
      # sink(zz,type = "message") #  c("output"))
      # sink(zz,type = "output") #  c("output"))    
      # cat("Thread", This_thread,"\n")
      cat("List", Range_i,"\n") #,"/",segments,
      
      LD_sub<-lapply(All_Pop, function(x){
        
        LD_sub_i<-LDmatrix_catch(
          c(GWAS_sub,eQTL_sub),
          pop = x,#"EUR",# "CEU",
          r2d = "r2",
          token = token_current,
          file = FALSE,
          genome_build = genome_build_sel
        )
        if( all(!(eQTL_sub %in% colnames(LD_sub_i))) ){
          data.frame(matrix(NA,0,4))%>%setNames(c("POP", "RS_number","eQTL_SNP","r2"))
        }else{
          LD_sub_i<-LD_sub_i[,colnames(LD_sub_i) %in% c("RS_number",eQTL_sub)]
          LD_sub_i<-melt(LD_sub_i, id.vars = "RS_number", variable.name = "eQTL_SNP", value.name = "r2")
          data.frame(POP=x, LD_sub_i)
        }
      })
      LD_sub<-do.call(rbind.data.frame,LD_sub)
      if (nrow(LD_sub)>0){
        LD_sub<-dcast(LD_sub, eQTL_SNP+RS_number~POP, value.var = "r2")
        LD_sub<-data.frame(Range=Range_i,LD_sub)
      }
      saveRDS(LD_sub, paste0(LD_PATH, "/SNP_LD_sub_",Range_i,".RDS"))
      closeAllConnections()
    }else{
      SNP_add_ct<-API_segment_length-length(eQTL_sub)
      segments_sub<-ceiling(length(GWAS_sub)/SNP_add_ct)
      separator<-rep(1:segments_sub, each=SNP_add_ct)[1:length(GWAS_sub)]
      GWAS_sub_list<-split(GWAS_sub, f = separator)
      
      Done_part_current<-list.files(LD_PATH)
      Done_part_current<-Done_part_current[grepl(paste0("SNP_LD_sub_",Range_i,"_part.*\\.RDS"),Done_part_current)]
      Done_part_current<-gsub("\\.RDS","", gsub(paste0("SNP_LD_sub_",Range_i,"_part"),"",Done_part_current)) %>% as.numeric
      
      if(all(1:segments_sub %in% Done_part_current)){
        LD_sub <- lapply(1:segments_sub, function(x)
          readRDS(paste0(LD_PATH, "/SNP_LD_sub_",Range_i,"_part",x, ".RDS"))
        )
        #LD_sub <- ListToDataframe_mixed(LD_sub, "r")
        LD_sub <- do.call(rbind.data.frame,LD_sub)
        saveRDS(LD_sub, paste0(LD_PATH, "/SNP_LD_sub_",Range_i, ".RDS"))  
        closeAllConnections()
      }else{
        
        # zz <- file(paste0(PATH, SUB_PATH, "Errorlog/errortag_",Range_i,".Rout"), open = "at")
        # sink(zz,type = "message",append = T) #  c("output"))
        # sink(zz,type = "output", append = T) #  c("output"))    
        #cat("Thread", This_thread,"\n")
        cat("List", Range_i,"\n") #"/",segments,
        
        if(! (1 %in% Done_part_current) ){
          i=1 #Do i=1 to see if GWAS SNP is in the database at all.
          cat("sub_list", i,"/",length(GWAS_sub_list))
          GWAS_sub_sub<-GWAS_sub_list[[i]]
          LD_sub_sub<-lapply(All_Pop, function(x){
            LD_sub_i<-LDmatrix_catch(
              c(eQTL_sub,GWAS_sub_sub),
              pop = x,#"EUR",# "CEU",
              r2d = "r2",
              token = token_current,
              file = FALSE,
              genome_build = genome_build_sel
            )
            if( all(!(eQTL_sub %in% colnames(LD_sub_i))) ){
              data.frame(matrix(NA,0,4))%>%setNames(c("POP", "RS_number","eQTL_SNP","r2"))
            }else{
              LD_sub_i<-LD_sub_i[,colnames(LD_sub_i) %in% c("RS_number",eQTL_sub)]
              LD_sub_i<-melt(LD_sub_i, id.vars = "RS_number", variable.name = "eQTL_SNP", value.name = "r2")
              data.frame(POP=x, LD_sub_i)
            }
          })
          LD_sub_sub<-do.call(rbind.data.frame,LD_sub_sub)
          if(nrow(LD_sub_sub)>0){ # Only proceed if GWAS SNP is in the database
            LD_sub_sub<-dcast(LD_sub_sub, eQTL_SNP+RS_number~POP, value.var = "r2")
            LD_sub_sub<-data.frame(Range=Range_i,LD_sub_sub)
            saveRDS(LD_sub_sub, paste0(LD_PATH, "/SNP_LD_sub_",Range_i,"_part",i, ".RDS"))  
            
            
            for(i in 2:segments_sub) {
              cat("sub_list", i,"/",segments_sub)
              GWAS_sub_sub<-GWAS_sub_list[[i]]
              LD_sub_sub<-lapply(All_Pop, function(x){
                
                LD_sub_i<-LDmatrix_catch(
                  c(eQTL_sub,GWAS_sub_sub),
                  pop = x,#"EUR",# "CEU",
                  r2d = "r2",
                  token = token_current,
                  file = FALSE,
                  genome_build = genome_build_sel
                )
                if( all(!(eQTL_sub %in% colnames(LD_sub_i))) ){
                  data.frame(matrix(NA,0,4))%>%setNames(c("POP", "RS_number","eQTL_SNP","r2"))
                }else{
                  LD_sub_i<-LD_sub_i[,colnames(LD_sub_i) %in% c("RS_number",eQTL_sub)]
                  LD_sub_i<-melt(LD_sub_i, id.vars = "RS_number", variable.name = "eQTL_SNP", value.name = "r2")
                  data.frame(POP=x, LD_sub_i)
                }
              })
              LD_sub_sub<-do.call(rbind.data.frame,LD_sub_sub)
              LD_sub_sub<-dcast(LD_sub_sub, eQTL_SNP+RS_number~POP, value.var = "r2")    
              LD_sub_sub<-data.frame(Range=Range_i,LD_sub_sub)
              saveRDS(LD_sub_sub, paste0(LD_PATH, "/SNP_LD_sub_",Range_i,"_part",i, ".RDS"))  
            }
            
            Done_part_current<-list.files(LD_PATH) 
            Done_part_current<-Done_part_current[grepl(paste0("SNP_LD_sub_",Range_i,"_part.*\\.RDS"),Done_part_current )]
            
            Done_part_current<-gsub("\\.RDS","", gsub(paste0("SNP_LD_sub_",Range_i,"_part"),"",Done_part_current)) %>% as.numeric
            if(all(1:segments_sub) %in% Done_part_current){
              LD_sub <- lapply(1:segments_sub, function(x)
                readRDS(paste0(LD_PATH, "/SNP_LD_sub_",Range_i,"_part",x, ".RDS"))
              )
              #LD_sub <- ListToDataframe_mixed(LD_sub, "r")
              LD_sub <- do.call(rbind.data.frame,LD_sub)
              saveRDS(LD_sub, paste0(LD_PATH, "/SNP_LD_sub_",Range_i, ".RDS"))  
              closeAllConnections()
            }
            
          }else{
            cat("No r2 return for GWAS SNP", Range_i,"Jump out of loop!\n")
            saveRDS(LD_sub, paste0(LD_PATH, "/SNP_LD_sub_",Range_i,".RDS"))  
            closeAllConnections()
          }
        }else{
          for(i in (2:segments_sub)[!(2:segments_sub) %in% Done_part_current] ) {
            cat("sub_list", i,"/",segments_sub)
            GWAS_sub_sub<-GWAS_sub_list[[i]]
            LD_sub_sub<-lapply(All_Pop, function(x){
              
              LD_sub_i<-LDmatrix_catch(
                c(eQTL_sub,GWAS_sub_sub),
                pop = x,#"EUR",# "CEU",
                r2d = "r2",
                token = token_current,
                file = FALSE,
                genome_build = genome_build_sel
              )
              if( all(!(eQTL_sub %in% colnames(LD_sub_i))) ){
                data.frame(matrix(NA,0,4))%>%setNames(c("POP", "RS_number","eQTL_SNP","r2"))
              }else{
                LD_sub_i<-LD_sub_i[,colnames(LD_sub_i) %in% c("RS_number",eQTL_sub)]
                LD_sub_i<-melt(LD_sub_i, id.vars = "RS_number", variable.name = "eQTL_SNP", value.name = "r2")
                data.frame(POP=x, LD_sub_i)
              }
            })
            LD_sub_sub<-do.call(rbind.data.frame,LD_sub_sub)
            LD_sub_sub<-dcast(LD_sub_sub, eQTL_SNP+RS_number~POP, value.var = "r2")    
            LD_sub_sub<-data.frame(Range=Range_i,LD_sub_sub)
            saveRDS(LD_sub_sub, paste0(LD_PATH, "/SNP_LD_sub_",Range_i,"_part",i, ".RDS"))  
          }
          
          Done_part_current<-list.files(paste0(LD_PATH)) 
          Done_part_current<-Done_part_current[grepl(paste0("SNP_LD_sub_",Range_i,"_part.*\\.RDS"),Done_part_current)]
          Done_part_current<-gsub("\\.RDS","", gsub(paste0("SNP_LD_sub_",Range_i,"_part"),"",Done_part_current)) %>% as.numeric
          if(all(1:segments_sub) %in% Done_part_current){
            LD_sub <- lapply(1:segments_sub, function(x)
              readRDS(paste0(LD_PATH, "/SNP_LD_sub_",Range_i,"_part",x, ".RDS"))
            )
            #LD_sub <- ListToDataframe_mixed(LD_sub, "r")
            LD_sub <- do.call(rbind.data.frame,LD_sub)
            saveRDS(LD_sub, paste0(LD_PATH, "/SNP_LD_sub_",Range_i, ".RDS"))  
            closeAllConnections()
          }
        }
      }      
    }
  }
  #}
  closeAllConnections()
  #}
  
  #stopCluster(cl)
  if(rm_temp_file){
    Done_part_current<-list.files(paste0(LD_PATH)) 
    Done_part_current<-Done_part_current[grepl(paste0("_part.*\\.RDS"),Done_part_current)]  
    unlink(paste0(LD_PATH,Done_part_current))
  }
  
  
  Done_part_current<-list.files(paste0(LD_PATH))   
  Done_part_current<-Done_part_current[grepl(paste0("SNP_LD_sub_.*\\.RDS"),Done_part_current)& !grepl(paste0("part*\\.RDS"),Done_part_current)]  
  
  Results<-lapply(Done_part_current, function(x) readRDS(paste0(LD_PATH,"",x)) )
  Results<-Results[unlist(lapply(Results,function(x) !is.null(x)))]
  #Results<-ListToDataframe_mixed(Results, "r") %>% unique
  Results<-do.call(rbind.data.frame, Results) %>% unique
  Results$eQTL_SNP<-paste(Results$eQTL_SNP)
  #saveRDS(Results, paste0(PATH, SUB_PATH,"eQTL_LDLink_All.RDS"))
  Results<-Results[!is.na(Results$EUR),]
  Results08<-Results[Results[,All_Pop]>=r2_threshold,]
  
  
  
  
  #Add missing GWAS SNPs back
  LDlink_missed_SNP<-GWAS$snpid[GWAS$snpid%in%eQTL_data$snps & !GWAS$snpid%in%Results08$RS_number ] %>% unique
  if(length(LDlink_missed_SNP)>0){
    LDlink_missed_SNP_dataframe<-data.frame(GWAS$Range[match(LDlink_missed_SNP, GWAS$snpid)], LDlink_missed_SNP, LDlink_missed_SNP, 1) %>% setNames(colnames(Results08))
    Results08<-rbind.data.frame(Results08, LDlink_missed_SNP_dataframe)
  }
  
  saveRDS(Results08, paste0(LD_PATH, "/LDLink_results.RDS"))
  #Results08 <- readRDS( paste0(LD_PATH, "/LDLink_results.RDS"))
  
  eQTL_sub<-eQTL_data[eQTL_data$snps%in%Results08$eQTL_SNP,]
  eQTL_sub$Range<-NULL
  GWAS_sub<-GWAS[GWAS$snpid%in%Results08$RS_number,]
  GWAS_sub$Range<-NULL
  
  colnames(GWAS_sub)[-1]<-paste0("GWAS_",colnames(GWAS_sub)[-1])
  colnames(GWAS_sub)<-gsub("GWAS_GWAS","GWAS",colnames(GWAS_sub))
  
  
  Results08<-Results08[Results08$RS_number %in% GWAS_sub$snpid & Results08$eQTL_SNP %in% eQTL_sub$snps,-1]
  colnames(Results08)<-c("snps","snpid", "eQTL_GWAS_SNP_r2_EUR")
  
  data_all<-merge(eQTL_sub, Results08, by="snps")
  data_all<-merge(data_all, GWAS_sub, by="snpid")
  data_all<-dplyr::rename(data_all, "GWAS_snp"="snpid")
  
  write.table(data_all,sep="\t",row.names = F, paste0(PATH, "Output/eQTL_sig_Olap_GWAS.tsv"))
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

Plot_GeneTrack<-function(Genes_of_interest, arrow_freq=NULL, yloc=-10, nameloc=-15,fontsize=5, window=5000,
                         Gene_annotation_file,
                         Exon_annotation_file,
                         NameConver=NULL
){
  require(ggplot2)
  require(dplyr)
  Gene_annotation<-read.delim(sep="\t",stringsAsFactors = F, Gene_annotation_file) 
  Exon_annotation<-read.delim(sep="\t",stringsAsFactors = F, Exon_annotation_file)
  
  if(!is.null(NameConver)){
    Gene_annotation<-Gene_annotation_NameConver(Gene_annotation, NameConver)
    Exon_annotation<-Gene_annotation_NameConver(Exon_annotation, NameConver)
  }
  
  
  Gene_loc<-Gene_annotation[Gene_annotation$Gene_Name%in%Genes_of_interest,]
  Gene_loc$Start_promoter<-Gene_loc$Start-window
  Gene_loc$Chr<-paste0("chr",Gene_loc$Chr)
  #Gene_loc_GR<-makeGRangesFromDataFrame(Gene_loc, keep.extra.columns = T, ignore.strand = F,seqinfo = Seqinfo(genome="hg38"),seqnames.field = "Chr",start.field = "Start_promoter", end.field = "End",strand.field = "strand")
  
  Exon_loc<-Exon_annotation[Exon_annotation$Gene_Name%in%Genes_of_interest,]
  Exon_loc<-Exon_loc[order(Exon_loc$Ensembl_id,Exon_loc$Chr,Exon_loc$Start,Exon_loc$End),]
  Exon_loc$Chr<-paste0("chr",Exon_loc$Chr)
  #Exon_loc_GR<-makeGRangesFromDataFrame(Exon_loc, keep.extra.columns = T, ignore.strand = F,seqinfo = Seqinfo(genome="hg38"),seqnames.field = "Chr",start.field = "Start", end.field = "End",strand.field = "strand")
  Exon_loc_list<-lapply(unique(Exon_loc$Gene_Name), function(x) Exon_loc[Exon_loc$Gene_Name%in%x,] ) %>% setNames(unique(Exon_loc$Gene_Name))
  
  
  #### Generate gene-segments so that arrows can be plotted  
  Gene_loc_segs<-list()
  for(Gene in unique(Gene_loc$Gene_Name)){
    Gene_loc_seg<-Gene_loc[Gene_loc$Gene_Name%in%Gene,]
    
    if(is.null(arrow_freq)){
      arrow_freq_current<- abs(Gene_loc_seg$Start-Gene_loc_seg$End)/6
    }else{
      arrow_freq_current<-arrow_freq
    }
    
    if(Gene_loc_seg$strand=="+"){
      start_seg<-c(seq(Gene_loc_seg$Start,Gene_loc_seg$End-50, by = arrow_freq_current),Gene_loc_seg$End-50) %>% unique
    }else{
      start_seg<-c(seq(Gene_loc_seg$End,Gene_loc_seg$Start+50, by = -arrow_freq_current),Gene_loc_seg$Start+50) %>% unique
      #Gene_loc_seg<-data.frame(Gene_Name=Gene, Start=start_seg[2:(length(start_seg))],End=start_seg[1:(length(start_seg)-1)], strand=Gene_loc_seg$strand)
    }
    Gene_loc_seg<-data.frame(Gene_Name=Gene, Start=start_seg[1:(length(start_seg)-1)],End=start_seg[2:(length(start_seg))], 
                             strand=Gene_loc_seg$strand)
    Gene_loc_segs[[Gene]]<-Gene_loc_seg
  }
  
  Plot<-list()
  for(Gene in names(Gene_loc_segs)){
    #wh <- range(wh, ignore.strand = T)
    Plot[[Gene]]<-ggplot()+
      geom_segment(data=Exon_loc_list[[Gene]], aes_(x=Exon_loc_list[[Gene]]$Start,y=yloc, xend=Exon_loc_list[[Gene]]$End,yend=yloc), color="black", size=8)+
      geom_segment(data=Gene_loc_segs[[Gene]], aes_(x=Gene_loc_segs[[Gene]]$Start,y=yloc, xend=Gene_loc_segs[[Gene]]$End,yend=yloc), color="black", size=1.5,
                   arrow = arrow(length = unit(0.02,"npc")))+
      geom_text(aes_(x=mean(Gene_loc_segs[[Gene]]$Start),y=nameloc, label=Gene), size=fontsize)
  }
  return(list(Plot=Plot,Gene_loc=Gene_loc,Exon_loc_list=Exon_loc_list,Gene_loc_segs=Gene_loc_segs ))
}



theme_neat<-function(text_size,...){
  require(ggplot2)
  themes<-theme_bw()+theme(
    text = element_text(size=text_size),#, family="Arial"),
    axis.text = element_text(size=text_size),#, family="Arial"), 
    legend.text=element_text(size=text_size),#, family="Arial"), 
    panel.grid= element_blank(),
    panel.border = element_rect(color="black",fill=NA),...)
  return(themes)
}


Prepare_data_for_circosplot<-function(eQTL_file,Loci_file,SNP_data, OUT_PATH, cisDist=2e5){
  eQTL_data<-read.delim(sep="\t", stringsAsFactors = F, eQTL_file)
  Loci_pick<-read.delim(sep="\t", stringsAsFactors = F, Loci_file)
  
  snp_data_all<-read.table(sep="\t",header=T,SNP_data)
  rownames(snp_data_all)<-snp_data_all$snpid
  snp_data_all$snpid<-NULL
  
  eQTL_sub<-eQTL_data
  eQTL_sub$chr<-gsub("chr","",eQTL_sub$chr)
  eQTL_sub$chr<-paste0("chr",eQTL_sub$chr)
  eQTL_sub$snploc<-as.numeric(eQTL_sub$snploc)
  eQTL_sub$genestart<-as.numeric(eQTL_sub$genestart)
  eQTL_sub$geneend<-as.numeric(eQTL_sub$geneend)
  
  Loci<-Loci_pick 
  Loci$eQTL<-paste(Loci$gene, Loci$snps,sep="-")
  Loci$chr<-gsub("chr","",Loci$chr)
  Loci$chr<-paste0("chr",Loci$chr)
  Loci$GWAS_col <- c("darkred", "chocolate1","mistyrose","gold","darkolivegreen1")[factor(Loci$GWAS_PT)]
  
  eQTL_sub<-eQTL_sub[eQTL_sub$gene%in%Loci$gene,]
  eQTL_sub$eQTL <- paste(eQTL_sub$gene, eQTL_sub$snps, sep="-")
  GWAS <- Loci[,c("snps","GWAS_PT","GWAS_col")]%>%unique
  eQTL_sub<-merge(eQTL_sub,GWAS, by="snps", all.x=T)
  
  eQTL_sub$Selected <- eQTL_sub$eQTL %in% Loci$eQTL
  eQTL_sub$genestart_low<-eQTL_sub$genestart - cisDist
  eQTL_sub$geneend_high<-eQTL_sub$geneend + cisDist
  eQTL_sub<-eQTL_sub[order(eQTL_sub$chr,eQTL_sub$snploc),]
  
  Gene_Range<- eQTL_sub[,c("gene","chr", "genestart_low","geneend_high")] %>%unique
  Gene_Range$Range<-1:nrow(Gene_Range)
  eQTL_sub$Range <- Gene_Range$Range[match(eQTL_sub$gene, Gene_Range$gene)]
  
  
  eQTL_sub$cor<-NA_real_
  #Corr_data<-list()
  for(Range in unique(eQTL_sub$Range)){ #Range=1
    eQTL_sub_sub<-eQTL_sub[eQTL_sub$Range%in%Range,]
    snp_data<-snp_data_all[rownames(snp_data_all)%in%eQTL_sub_sub$snps,] %>% t
    snp_select<-eQTL_sub_sub[eQTL_sub_sub$Selected,]
    snp_select<-paste(snp_select[order(snp_select$pvalue),][which(snp_select$pvalue==min(snp_select$pvalue))[1], "snps"])
    cor_result<-cor(snp_data[,snp_select], snp_data, use = "pairwise.complete.obs")
    eQTL_sub$cor<-case_when(
      eQTL_sub$snps%in%colnames(cor_result) ~ cor_result[,match(eQTL_sub$snps, colnames(cor_result) )], TRUE ~ eQTL_sub$cor
    )
  }
  
  eQTL_sub$col<-case_when(
    eQTL_sub$Selected ~ "purple", 
    eQTL_sub$cor > 0.8 ~ "red",
    eQTL_sub$cor > 0.6 ~ "orange",
    eQTL_sub$cor > 0.4 ~ "green",
    eQTL_sub$cor > 0.2 ~ "cyan",
    TRUE ~ "grey75"
  )
  
  eQTL_sub$size<-case_when(
    eQTL_sub$Selected ~ 1.5, 
    eQTL_sub$cor>0.2 ~ 1.1,
    TRUE ~ 0.5
  )
  
  write.table(eQTL_sub, quote = F, row.names = F, sep="\t", paste0(OUT_PATH, "/eQTL_all_wCorr.tsv" ))
}


Regional_plot<-function(eQTL_file, Gene_annotation_file,Exon_annotation_file,
                        cisDist=2e+05, yloc = -80, nameloc = -100, eQTL_size_scale=3.33){
  require(GenomicRanges)
  seq_info<-Seqinfo(genome="hg19")
  eQTL<-read.delim(eQTL_file, sep="\t", stringsAsFactors = F)
  eQTL$size<-eQTL$size*eQTL_size_scale
  
  col_CHM = structure( names=c("TSS", "Enhancer","Enhancer_Mixed","Transcribed","Repressed","Ambiguous"), c(gg_color_hue(5),"grey30") )
  
  eQTL<-split(eQTL, f=eQTL$gene)
  
  Gene_annotation<-read.delim(sep="\t", stringsAsFactors = F, Gene_annotation_file)
  Exon_annotation<-read.delim(sep="\t", stringsAsFactors = F, Exon_annotation_file)
  
  
  All_Gene<-names(eQTL)
  Gene_loc<-Gene_annotation[Gene_annotation$Gene_Name%in%All_Gene,]
  Gene_loc$Start_ext<-Gene_loc$Start-cisDist
  Gene_loc$End_ext<-Gene_loc$End+cisDist
  Gene_loc$Chr<-gsub("chr","",Gene_loc$Chr)
  Gene_loc$Chr<-paste0("chr",Gene_loc$Chr)
  Gene_loc_GR<-makeGRangesFromDataFrame(Gene_loc, keep.extra.columns = T, ignore.strand = F,seqinfo = seq_info,seqnames.field = "Chr",start.field = "Start_ext", end.field = "End_ext",strand.field = "strand")
  
  
  
  GeneTrack<-Plot_GeneTrack(All_Gene,yloc = yloc, nameloc = nameloc,
                            Gene_annotation_file = Gene_annotation_file, 
                            Exon_annotation_file = Exon_annotation_file, arrow_freq = NULL)
  
  
  
  Plots<-list()
  
  for(Gene in All_Gene){ 
    Gene_loc_GR_sub<-Gene_loc_GR[Gene_loc_GR$Gene_Name%in%Gene]
    eQTL_sub<-eQTL[[Gene]]
    eQTL_sub$logP<- -log10(eQTL_sub$pvalue)
    eQTL_sub_GR<-makeGRangesFromDataFrame(eQTL_sub, seqinfo = seq_info, seqnames.field = "chr", start.field = "snploc", end.field = "snploc",ignore.strand = T, keep.extra.columns = T )
    
    
    eQTL_sub_top<-eQTL_sub[eQTL_sub$col=="purple",,drop=F]
    
    ChIP_df_sub<-eQTL_sub[eQTL_sub$ChIP,]
    ChIP_df_sub$y<-15 
    ChIP_df_sub$col<-"blue"
    CHM_GR_sub<-eQTL_sub[!is.na(eQTL_sub$CHM),]
    
    
    Plots[[Gene]]<-
      GeneTrack$Plot[[Gene]]+
      
      geom_point(data=eQTL_sub, aes_(x=eQTL_sub$snploc, y=eQTL_sub$logP*100), color=eQTL_sub$col, size=eQTL_sub$size)+ #eQTL track
      geom_point(data=eQTL_sub_top, aes_(x=eQTL_sub_top$snploc, y=eQTL_sub_top$logP*100), color=eQTL_sub_top$col, size=eQTL_sub_top$size)+ #TopSNP always on top
      geom_segment(data=ChIP_df_sub, aes_(x=ChIP_df_sub$ChIP_start-5, xend=ChIP_df_sub$ChIP_end+5, y=-25, yend=-25),color=ChIP_df_sub$col, size=6)+ #ChIPseq track (All ChIP)
      geom_hline(yintercept = -5, size=1, color="black")+
      theme_neat(15)+
      scale_y_continuous(breaks = seq(0,2000, 100), labels = seq(0,20, 1) ) + 
      labs(title=Gene, x=eQTL_sub$chr[1], y="-log10(p)" )
    
    if(nrow(CHM_GR_sub)>0){
      Plots[[Gene]]<- Plots[[Gene]] + geom_segment(data=CHM_GR_sub, aes_(x=CHM_GR_sub$snploc-5, xend=CHM_GR_sub$snploc+5, y=-45, yend=-45),color=col_CHM[match(CHM_GR_sub$CHM,names(col_CHM))], size=6) #CHM Track All-short
    }  
    
    
    
  }
  
  pdf(paste0(PATH, "Plots/Regional_plots.pdf"), width=8, height=6)
  print(Plots)
  dev.off()
}


Plot_circos<-function(eQTL_file,
                      SNP_file,
                      Exp_file,
                      Plot_PATH,
                      gene_name=T,
                      cisDist=2e+05,
                      CHM_extend=50,
                      CHIP_extend=100,
                      eps=F){
  library(GenomicRanges)
  library(circlize)
  library(dendextend)
  library(ComplexHeatmap)
  
  circos.clear()
  hg19_chr<-Seqinfo(genome="hg19") #readRDS(paste0(PATH, "../Refdata/hg19_seqinfo.RDS"))
  col_CHM = structure( names=c("TSS", "Enhancer","Enhancer_Mixed","Transcribed","Repressed","Ambiguous"), c(gg_color_hue(5),"grey30") )
  col_fun1 = colorRamp2(c(-0.5, 0, 0.5), c("blue3", "white", "red3"))
  
  
  eQTL_data<- read.delim(sep="\t", stringsAsFactors = F, eQTL_file)
  eQTL_data$CHM_col <- col_CHM[match(eQTL_data$CHM,names(col_CHM))]
  snp_data_all<-read.table(sep="\t",header=T,row.names = 1,SNP_file)
  gene_data_all<-read.table(sep="\t",header=T,row.names = 1,Exp_file)
  if(!all(colnames(gene_data_all)==colnames(snp_data_all))){stop("Too bad, gene/snp colnames don't match")}
  
  
  eQTL_sub<-eQTL_data
  eQTL_sub$Range_start<-eQTL_sub$genestart - cisDist
  eQTL_sub$Range_end<-eQTL_sub$geneend + cisDist
  eQTL_sub$logP<- -log10(eQTL_sub$pvalue)
  
  # #Expression fold change heatmap Ranges
  eQTL_sub_sub<-eQTL_sub[eQTL_sub$Selected,]
  eQTL_sub_sub<-lapply(unique(eQTL_sub_sub$Range), function(x) {
    temp<-eQTL_sub_sub[eQTL_sub_sub$Range%in%x,];temp[order(temp$pvalue),][1,]}  )
  eQTL_sub_sub<-do.call(rbind.data.frame,eQTL_sub_sub)
  GWAS<-eQTL_sub_sub[!is.na(eQTL_sub_sub$GWAS_PT),]
  
  Heatmap_data<-list()
  for(i in 1:nrow(eQTL_sub_sub)){
    Heatmap_data[[i]]<-data.frame(eQTL_sub_sub$gene[i],eQTL_sub_sub$snps[i], eQTL_sub_sub$eQTL[i],
                                  scale(t(gene_data_all[paste(eQTL_sub_sub$gene[i]), ])),t(snp_data_all[paste(eQTL_sub_sub$snps[i]), ]), stringsAsFactors = F) %>% setNames(c("gene","snps","eQTL", "Exp","Genotype"))
  }
  Heatmap_data<-do.call(rbind.data.frame, Heatmap_data)
  Heatmap_data<-Heatmap_data[!is.na(Heatmap_data$Genotype),]
  Heatmap_data$Range <- eQTL_sub_sub$Range[match(Heatmap_data$eQTL, eQTL_sub_sub$eQTL)]
  Heatmap_data_summary<-dcast(eQTL ~ Genotype, data=Heatmap_data,value.var = "Exp", fun.aggregate = mean)
  Heatmap_Range<-Heatmap_data_summary$Range<-eQTL_sub_sub$Range[match(Heatmap_data_summary$eQTL, eQTL_sub_sub$eQTL)]
  
  Heatmap_data_summary_0<-recast(Range ~ variable, data=Heatmap_data_summary[,-1],id.var = "Range", fun.aggregate=mean)
  rownames(Heatmap_data_summary_0)<-Heatmap_data_summary_0$Range
  Heatmap_data_summary_0$Range<-NULL
  Heatmap_data_summary_0<-Heatmap_data_summary_0[order.hclust(hclust(d=dist(Heatmap_data_summary_0))),]
  
  Heatmap_data_list<-list()
  for(Range in rownames(Heatmap_data_summary_0)){
    eQTL<-Heatmap_data$eQTL[Heatmap_data$Range%in%Range]
    Heatmap_data_list[[Range]]<-Heatmap_data_summary[Heatmap_data_summary$eQTL%in%eQTL,]
  }
  Heatmap_data_list<-do.call(rbind.data.frame, Heatmap_data_list)
  
  Heatmap_data_summary<-Heatmap_data_list
  colnames(Heatmap_data_summary)[colnames(Heatmap_data_summary) %in% 0:2] <-
    paste0("G",colnames(Heatmap_data_summary[,colnames(Heatmap_data_summary) %in% 0:2]))
  
  
  eQTL_sub$Range<-factor(eQTL_sub$Range, levels = unique(sort(Heatmap_data_summary$Range)))
  
  if(eps){
    setEPS()
    postscript(file = paste0(Plot_PATH,"Circos_plot.eps" ), height=10, width=10)
  }else{
    pdf(paste0(Plot_PATH,"/Circos_plot.pdf" ), height=10, width=10)
  }
  lgd_heatmap = Legend(at = c(-1, -0.5, 0, 0.5, 1), col_fun = col_fun1, title_position = "topleft", title = "Heatmap-\nTrack")
  
  lgd_eQTL = Legend(at = c("Top-SNP", "r2>0.8", "r2>0.6", "r2>0.4", "r2>0.2", "r2<=0.2"), type = "points", 
                    legend_gp = gpar(col = c("purple","red","orange","green","cyan","grey75")), title_position = "topleft", title = "eQTL-Track")
  
  #lgd_peak = Legend(at = c("Gene"), type = "lines", legend_gp = gpar(col = c("blue3"),lwd=2), title_position = "topleft",  title = "Gene-Track")
  lgd_peak = Legend(at = c("ChIP-Peak", "Gene"), type = "lines", legend_gp = gpar(col = c("red3","skyblue"),lwd=2), title_position = "topleft",  title = "ChIP/Gene-Track")
  
  lgd_CHM = Legend(at = c("TSS","Enhancer","Enhancer_Mixed","Transcribed","Repressed", "Ambiguous"), type = "lines", legend_gp = gpar(col = c(gg_color_hue(5),"grey30") ,lwd=2), title_position = "topleft",  title = "cis-element")
  
  lgd_GWAS = Legend(at = unique(GWAS$GWAS_PT), type = "lines", legend_gp = gpar(col = unique(GWAS$GWAS_col),lwd=2), title_position = "topleft",  title = "GWAS")
  
  lgd_list_vertical = packLegend(lgd_heatmap, lgd_eQTL, lgd_peak)
  lgd_list_vertical2 = packLegend(lgd_CHM,lgd_GWAS)
  
  circos.clear()
  circos.par("start.degree" = 0, "gap.degree"=0, "cell.padding"=c(0,0,0,0), track.margin=c(0.005, 0.005))
  circos.initialize(factors = eQTL_sub$Range, x=eQTL_sub$snploc)
  
  #GWAS track
  circos.par("track.height" = 0.025)
  circos.track(factors = eQTL_sub$Range, ylim = c(0, 1), bg.col="grey95",bg.border = "grey70")
  for(Range in unique(eQTL_sub$Range) ){ #Range = unique(eQTL_data$Range)[1]
    Loci_sub<-eQTL_sub[eQTL_sub$Range%in%Range & eQTL_sub$Selected,,drop=F]
    Loci_sub<-Loci_sub[order(Loci_sub$snploc),]
    n_GWAS<-length(unique(Loci_sub$GWAS_PT))
    
    xpos <- seq(Loci_sub$Range_start[1],Loci_sub$Range_end[1], (Loci_sub$Range_end[1]-Loci_sub$Range_start[1])/n_GWAS)
    
    
    circos.rect( xleft = xpos[1:n_GWAS], xright = xpos[2:(n_GWAS+1)],
                 ybottom =rep(0,n_GWAS) , ytop = rep(1,n_GWAS),
                 track.index = 1, sector.index = Range,col=unique(Loci_sub$GWAS_col),border = unique(Loci_sub$GWAS_col) )
  }
  
  #"heatmap" track
  circos.par("track.height" = 0.1)
  circos.track(factors = Heatmap_data_summary$Range, ylim = c(0, 1), bg.col="grey20")
  for(Range in paste(unique(Heatmap_data_summary$Range)) ){ 
    Heatmap_data_sub<-Heatmap_data_summary[Heatmap_data_summary$Range%in%Range,]
    
    ypos<-seq(1,0,-1/3)
    ypos<- data.frame(y0=ypos[1:3],y1=ypos[2:4],Genotype=paste0("G",0:2))
    
    xrange<-data.frame(start=min(eQTL_sub[eQTL_sub$Range%in%Range,"snploc"]), end=max(eQTL_sub[eQTL_sub$Range%in%Range,"snploc"]))
    xrange$width<-xrange$end-xrange$start+1
    
    xpos<-seq(xrange$start,xrange$end,(xrange$width-1)/nrow(Heatmap_data_sub))
    
    Heatmap_data_sub_melt<- data.frame(x0=xpos[1:(length(xpos)-1)] , x1=xpos[2:(length(xpos))], Heatmap_data_sub[,paste0("G",0:2)])
    
    Heatmap_data_sub_melt<-melt(Heatmap_data_sub_melt, id=c("x0","x1"))
    
    Heatmap_data_sub_melt<-cbind.data.frame(Heatmap_data_sub_melt, ypos[match(Heatmap_data_sub_melt$variable, ypos$Genotype),])
    Heatmap_data_sub_melt$col<-col_fun1(Heatmap_data_sub_melt$value)
    
    circos.rect( xleft = Heatmap_data_sub_melt$x0
                 , xright = Heatmap_data_sub_melt$x1, 
                 ybottom =Heatmap_data_sub_melt$y0, ytop =Heatmap_data_sub_melt$y1,
                 col=Heatmap_data_sub_melt$col,
                 track.index = 2, sector.index = Range )
    
  }
  
  #eQTL track
  circos.par("track.height" = 0.25)
  circos.track(factors = eQTL_sub$Range, y = eQTL_sub$logP, bg.col="white",bg.border = "grey50", 
               panel.fun = function(x, y) {
                 circos.segments(x0 = CELL_META$xlim[1], x1 = CELL_META$xlim[2], y0=c(1:4),y1=c(1:4), col="grey50", lty=3)
               })
  circos.trackPoints(eQTL_sub$Range, x=eQTL_sub$snploc, y=eQTL_sub$logP, col = "black", bg = eQTL_sub$col, pch = 21, cex = eQTL_sub$size,track.index = 3)
  
  
  # #CHM track
  circos.par("track.height" = 0.025)
  CHM_sub<-eQTL_sub[!is.na(eQTL_sub$CHM),]
  circos.track(factors = CHM_sub$Range, ylim = c(0, 1), bg.col="grey95",bg.border = "grey70")
  
  for(Range in unique(CHM_sub$Range) ){ 
    CHM_sub_sub<-CHM_sub[CHM_sub$Range%in%Range&!is.na(CHM_sub$Range),,drop=F]
    circos.rect( xleft =CHM_sub_sub$snploc-CHM_extend, xright =CHM_sub_sub$snploc+CHM_extend,
                 ybottom =rep(0,nrow(CHM_sub_sub)) , ytop = rep(1,nrow(CHM_sub_sub)),
                 track.index = 4, sector.index = Range,col=CHM_sub_sub$CHM_col,border = CHM_sub_sub$CHM_col )
  }
  
  #ChIP track
  circos.par("track.height" = 0.025)
  ChIP_sub<-eQTL_sub[(eQTL_sub$ChIP),]
  circos.track(factors = ChIP_sub$Range, ylim = c(0, 1), bg.col="grey95",bg.border = "grey70")
  for(Range in unique(ChIP_sub$Range) ){
    ChIP_sub_sub<-ChIP_sub[ChIP_sub$Range%in%Range,,drop=F]
    circos.rect( xleft =ChIP_sub_sub$ChIP_start-CHIP_extend, xright =ChIP_sub_sub$ChIP_end+CHIP_extend,
                 ybottom =rep(0.25,nrow(ChIP_sub_sub)) , ytop =rep(0.75,nrow(ChIP_sub_sub)),
                 track.index = 5, sector.index = Range, col="red3",border = "red3" )
  }
  
  
  
  #Gene track
  circos.par("track.height" = 0.05)
  Gene_Ranges<-eQTL_sub_sub
  circos.track(factors = Gene_Ranges$Range, ylim = c(0, 1), bg.col="grey95",bg.border = "grey70")
  for(Range in unique(Gene_Ranges$Range) ){ 
    Gene_Ranges_sub<-Gene_Ranges[Gene_Ranges$Range%in%Range,]
    
    ypos <- seq(0,1, 1/nrow(Gene_Ranges_sub))
    
    circos.rect( xleft =Gene_Ranges_sub$genestart, xright =Gene_Ranges_sub$geneend, 
                 ybottom =ypos[1:nrow(Gene_Ranges_sub)], ytop =ypos[2: (nrow(Gene_Ranges_sub)+1)], 
                 track.index = 6, sector.index = Range,col="skyblue",border = "skyblue" )
    if(gene_name){circos.text( x =rowMeans(Gene_Ranges_sub[,c("genestart", "geneend")]),
                               y =ypos[1:nrow(Gene_Ranges_sub)],
                               label = Gene_Ranges_sub$gene,
                               track.index = 6, sector.index = Range,
                               col="black", cex=1.5)}
  }
  draw(lgd_list_vertical, x = unit(1, "mm"), y = unit(1, "mm"), just = c("left", "bottom"))
  draw(lgd_list_vertical2, x = unit(255, "mm"), y = unit(10, "mm"), just = c("right", "bottom"))
  
  circos.clear()
  
  
  dev.off()
}