#############################################################
## Short tandem repeat data analysis
## Model 3 Outlier test
##
## This analysis includes the following steps:
## 01. STR count comparision
## 02. STR count comparision for each cutoff for amyloid beta positivity
## 03. Amyloid beta level comparision
## 04. Go term 
#############################################################



######## set up ########
# Input file prepare
{
  library(readxl)
  library(tidyverse)
  library(GenomicRanges)
  library(RNOmni)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  rm(list=ls())
  sessionInfo()
  
  setwd("~/Dropbox/ADWGS/STR/EH/")
  
  rename = dplyr::rename
  select = dplyr::select
  mutate = dplyr::mutate
  
  # Sample QC list by PCA
  pca_rm = read_tsv("eh.v5_coverage_240419_1558_PCA.txt")
  pca_rm_samples = pca_rm %>% filter(pc_remove == 1) %>% pull(sample)
  
  # Phenotype data
  pheno = read_tsv("~/Dropbox/SMC_AD_WGS_paper/Data/WGS_1824_phenotype_update.240416.txt") %>% dplyr::rename(sample = IID)
  pheno = pheno %>% select(sample,RCL40_visual, batch ,age,sex, DX,PC1,PC2,PC3,APOE, RCL_global, KlunkCL )
  pheno = pheno %>% filter(!sample %in% pca_rm_samples)
  pheno = pheno %>% filter(sample != 'WGS_1622')
  
  # Target chormosomes
  chr_list = paste0("chr",1:22)
  
  # bed file
  dat.annot<-read.delim("eh.v5_w_gangstr.v13.polymorphic_240419.bed", sep="\t", stringsAsFactors = F) #read in bed file of STRs
  dat.annot = dat.annot %>% filter(chr %in% chr_list)
  dat.annot %>% dim()
  
  # BED file from Guo et al. supplementary table 
  all_bed = readxl::read_xlsx('~/Dropbox/ADWGS/STR/ADGP_preprint/ADSP_supp_1.xlsx', sheet = 1) %>% rename('start' = "STR start position (GRCh38)", 'end' =  "STR end position (GRCh38)", 'ref_length' = "Reference tract length" ) %>% mutate(id = paste(chr ,   start  , motif,sep = '-'))
  
  
  # Genotype info
  dat<-data.table::fread("eh.v5_allele_long_240419.txt.gz", header=T, sep=",", data.table=F) 
  
  # STR coverage
  dat.site_cov<-data.table::fread("eh.v5_coverage_240419.txt.gz", header=T, sep=",", data.table=F) 
  
  # Sample mean coverage
  dat.sample_cov = data.table::fread("~/Dropbox/ADWGS/STR/EH/eh.v5.sample_coverage_240419.txt.gz")
  
}

# Segmental duplication region filtering
{
  dat.background = all_bed %>% filter(chr %in% chr)
  dat.background.se<-makeGRangesFromDataFrame(all_bed, seqnames = "chr", start.field = "start", end.field = "end",keep.extra.columns=TRUE) #Make into GRanges format
  #Subset out segmental duplication regions from background
  dat.segdup<-data.table::fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz", data.table=F, header=F, stringsAsFactors = F)
  dat.segdup<-dat.segdup[,c(2:4)]
  names(dat.segdup)<-c("chr", "start", "stop")
  dat.segdup.gr<-makeGRangesFromDataFrame(dat.segdup, seqnames.field = "chr", start.field = "start", end.field = "stop")
  overlap_index<-data.frame(findOverlaps(dat.background.se, dat.segdup.gr, type="any"))$queryHits
  dat.background$index<-seq(1,nrow(dat.background))
  dat.background$seg_dup<-ifelse(dat.background$index%in%overlap_index, 1, 0)
  dat.background<-subset(dat.background, seg_dup==0, select=-c(seg_dup,index))
  dat.background %>% dim()
  target = dat.background
  
}

# target str
{
  dat.annot.t = dat.annot %>% filter(id %in% target$id)
  dat.t = dat[dat.annot.t$index,]
  dat.site_cov.t = dat.site_cov[dat.annot.t$index,]
}

# phenotype & filtering
{
  sample_list = pheno %>% filter(DX %in% c("DAT","MCI","CU")) %>% pull(sample)
  pheno = pheno %>% mutate( is_case = case_when(
    RCL40_visual == 1 ~ 1,
    RCL40_visual == 0 ~ 0,
    T ~ NA
  ),
  is_female = case_when(
    sex == 2 ~ 1,
    sex == 1 ~ 0,
    T ~ NA
  ))
  
  # Remove bad quality sample
  sample_list <- sample_list[sample_list != "WGS_1622"]
  
  dat.t.test = dat.t %>% select(all_of(sample_list))
  dat.site_cov.t.test = dat.site_cov.t %>% select(all_of(sample_list))
  
}

# make test file
{
  pheno$is_female <- factor(pheno$is_female)
  pheno$batch <- factor(pheno$batch)
  
  dat.t.test.1 = t(dat.t.test)
  dat.t.test.s = dat.t.test.1 %>% as_tibble() %>% mutate('sample' = rownames(dat.t.test.1))
  
  dat.site_cov.t.test.1 = t(dat.site_cov.t.test)
  
  # Sample coverage
  dat.site_cov.t.test.s = dat.site_cov.t.test.1 %>% as_tibble() %>% mutate('sample' = rownames(dat.site_cov.t.test.1))
  dat.site_cov.t.test.s = dat.site_cov.t.test.s %>% left_join(dat.sample_cov)
  
  expan.test = read_tsv("outputs/eh_DBSCAN_outlier.output.tsv")
  
  expan.test$is_female <- as.factor(expan.test$is_female)
  expan.test$batch <- as.factor(expan.test$batch)
  
  expan.test <- expan.test %>% mutate(outlier_length_ref = outlier_length - ref)
  
  
  expan.test <- expan.test %>% left_join(dat.sample_cov)
  expan.test$s_avg_depth <- expan.test$eh.coverage.mean
  
  expan.test <- expan.test %>% mutate(id = paste(chr,pos_start ,motif, sep= '-'))
  
  expan.test.sample<-data.frame(table(expan.test$sample), stringsAsFactors = F)
  names(expan.test.sample)<-c("sample", "outlier_count")
  
  expan.test.sample = expan.test.sample %>% left_join(pheno)
  expan.test.sample <- expan.test.sample %>% left_join(dat.sample_cov)
  expan.test.sample$s_avg_depth <- expan.test.sample$eh.coverage.mean
}

########## 01. STR count comparision ########## 
{
  expan.test.sample$batch <- as.factor(expan.test.sample$batch)
  expan.test.sample$sex <- as.factor(expan.test.sample$sex)
  
  # For amyloid beta
  d1 = lm(is_case ~ outlier_count + is_female + batch + age + s_avg_depth +
            PC1 + PC2 + PC3, data = expan.test.sample )
  coef(summary(d1))
  
  expan.test.sample.DX <- expan.test.sample %>% filter(DX %in% c('DAT','CU')) %>%
    mutate(is_case = case_when(
      DX == 'DAT' ~ 1,
      DX == 'CU' ~ 0, 
      T ~ NA
    ))
  
  d2 = lm(is_case ~ outlier_count + is_female + batch + age + s_avg_depth +
            PC1 + PC2 + PC3, data = expan.test.sample.DX)
  coef(summary(d2))
  

}


########## 02. STR count comparision for each cutoff for amyloid beta positivity ########## 
{
  n_Case = 895 # Amyloid beta positive
  n_Control = 620 # Amyloid beta negative
  
  
  expan.test.sample <- expan.test.sample %>% left_join(pheno %>% select(sample, RCL_global, KlunkCL))
  
  
  threshold_list <- c(10,20,30,40)
  
  result_df <- data.frame(threshold = numeric(),
                          ratio = numeric(), stringsAsFactors = FALSE)
  
  
  for(i in threshold_list){
    if(i == 10){
      expan.test.sample.case = expan.test.sample %>% filter(outlier_count < i & is_case == 1) %>% pull(sample) %>% unique %>% length()
      expan.test.sample.ctrl = expan.test.sample %>% filter(outlier_count < i & is_case == 0) %>% pull(sample) %>% unique %>% length()
      
      ratio <- (expan.test.sample.case / expan.test.sample.ctrl) / ((n_Case - expan.test.sample.case) / (n_Control - expan.test.sample.ctrl))
      result_df <- rbind(result_df, data.frame(threshold = paste0("<",i), ratio = ratio))
      
      # expan.test.sample.case = expan.test.sample %>% filter(outlier_count >= i & is_case == 1) %>% pull(sample) %>% unique %>% length()
      # expan.test.sample.ctrl = expan.test.sample %>% filter(outlier_count >= i & is_case == 0) %>% pull(sample) %>% unique %>% length()
      #  
      # ratio <- (expan.test.sample.case / expan.test.sample.ctrl) / ((n_Case - expan.test.sample.case) / (n_Control - expan.test.sample.ctrl))
      # result_df <- rbind(result_df, data.frame(threshold = paste0('≥',i), ratio = ratio))
      next
      
    }
    expan.test.sample.case = expan.test.sample %>% filter(outlier_count >= i & is_case == 1) %>% pull(sample) %>% unique %>% length()
    expan.test.sample.ctrl = expan.test.sample %>% filter(outlier_count >= i & is_case == 0) %>% pull(sample) %>% unique %>% length()
    
    ratio <- (expan.test.sample.case / expan.test.sample.ctrl) / ((n_Case - expan.test.sample.case) / (n_Control - expan.test.sample.ctrl))
    result_df <- rbind(result_df, data.frame(threshold = paste0('≥',i), ratio = ratio))
  }
  result_df

}


########## 03. Amyloid beta level comparision ########## 
{
  expan.test.sample <- expan.test.sample %>%
    # filter(is_case == 1) %>%
    mutate(
      over40 = ifelse(outlier_count >= 40, 'Over_40_STR_outliers',
                      ifelse(outlier_count < 10, 'Under_10_STR_outliers',NA))
    )  #%>%filter(!is.na(over40))
  expan.test.sample %>% dim
  
  expan.test.sample.target = expan.test.sample %>% filter(!is.na(over40))
  expan.test.sample.target$batch <- factor(expan.test.sample.target$batch)
  expan.test.sample.target$sex <- factor(expan.test.sample.target$sex)
  
  ###### RCL_global 
  summary(lm(  RCL_global ~ over40 + is_female + batch + age + s_avg_depth +
                     PC1 + PC2 + PC3 ,expan.test.sample.target))
  
  
  ###### RCL_global 
  summary(lm(  RCL_global ~ over40 + is_female + batch + age + s_avg_depth +
                     PC1 + PC2 + PC3 ,expan.test.sample.target %>% filter(DX == 'CU')))
  
  
}

########## 04. Go term ########## 
{
  over40_sample = expan.test.sample %>% filter(over40 == 'Over_40_STR_outliers') %>% pull(sample)
  
  over40_str = expan.test %>% filter(sample %in% over40_sample)
  
  over40_str_list = over40_str %>% pull(id)
  
  bed_over40 <- all_bed %>% filter(id %in% over40_str_list)
  
  gene_id <- bed_over40 %>%
    filter(abs(`Distance to nearest TSS`) < 500000) %>%
    pull(`Gene symbol`) 
  
  gene_list <- bitr(gene_id, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db, drop = TRUE) %>% pull(ENTREZID)
  
  
  ## Background STR
  gene_id_back <- dat.background %>%
    filter(abs(`Distance to nearest TSS`) < 500000) %>%
    filter(id %in% dat.annot.t$id) %>%
    pull(`Gene symbol`)
  
  gene_list_back <- bitr(gene_id_back, fromType = "SYMBOL",
                         toType = c("ENTREZID"),
                         OrgDb = org.Hs.eg.db, drop = TRUE) %>% pull(ENTREZID)
  
  GO = enrichGO( gene_list,
                 OrgDb = org.Hs.eg.db,
                 keyType="ENTREZID", ont='ALL', pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                 universe = gene_list_back
  )
  
  GO_df = GO@result %>% as.data.frame() %>% arrange(desc(p.adjust))
  GO_df$Description = factor(GO_df$Description , levels = (GO_df$Description) )
  
  write_tsv(GO_df , "~/Dropbox/SMC_AD_WGS_paper/Data/STR/eh_STR_over40_goterm.tsv")
  
  GO_df %>%  arrange(p.adjust) %>% head
  
}


