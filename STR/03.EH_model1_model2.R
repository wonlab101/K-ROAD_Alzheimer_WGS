#############################################################
## Short tandem repeat data analysis
## ExpansionHunter output file preprocessing
##
## This analysis includes the following steps:
## 01. STR model 1 length test for amyloid beta
## 02. STR model 1 length test for clinical diagnosis
#############################################################

######## 02. STR model 1 length test for amyloid beta #######
# Input file prepare
{
  library(readxl)
  library(tidyverse)
  library(GenomicRanges)
  library(RNOmni)
  
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
  pheno = pheno %>% select(sample,RCL40_visual, batch ,age,sex, DX,PC1,PC2,PC3,APOE )
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

########## 0.1. target str ########## 
{
  dat.annot.t = dat.annot %>% filter(id %in% target$id)
  dat.t = dat[dat.annot.t$index,]
  dat.site_cov.t = dat.site_cov[dat.annot.t$index,]
}

####### 0.2 phenotype & filtering ####### 
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

####### 0.3 make test file ####### 
{
  pheno$is_female <- factor(pheno$is_female)
  pheno$batch <- factor(pheno$batch)
  
  dat.t.test.1 = t(dat.t.test)
  dat.t.test.s = dat.t.test.1 %>% as_tibble() %>% mutate('sample' = rownames(dat.t.test.1))
  
  dat.site_cov.t.test.1 = t(dat.site_cov.t.test)
  
  # Sample coverage
  dat.site_cov.t.test.s = dat.site_cov.t.test.1 %>% as_tibble() %>% mutate('sample' = rownames(dat.site_cov.t.test.1))
  dat.site_cov.t.test.s = dat.site_cov.t.test.s %>% left_join(dat.sample_cov)
  
  rm(dat,dat.annot,dat.site_cov,dat.site_cov.t,
     dat.site_cov.t.test,dat.site_cov.t.test.1,
     dat.t,dat.t.test,dat.t.test.1)
  
  
}

####### 1. Model 1 Length Test ####### 
{
  ###Run single STR association analysis
  library(data.table)
  library(stringr)
  library(RNOmni)
  library(jsonlite)
  
  # Create APOE4 columns
  pheno <- pheno %>%
    mutate(
      E4_count = str_count(APOE, "E4"),  # Count the number of E4 alleles (0, 1, 2)
      has_E4 = ifelse(E4_count > 0, 1, 0)  # 1 if E4 is present, else 0
    )
  
  dat.annot.t %>% head()
  dat.test = dat.annot.t 
  
  str_list = dat.test$index
  
  dat.test$allele_count_logistic_p<-1
  dat.test$allele_count_effect<-0
  dat.test$allele_count_mean_control<-0
  dat.test$allele_count_mean_ad<-0
  dat.test$allele_count_sd<-0
  
  dat.test$allele_count_inv_norm_logistic_p<-1
  dat.test$allele_count_inv_norm_effect<-0
  
  
  pheno$is_female <- factor(pheno$is_female)
  pheno$batch <- factor(pheno$batch)
  
  
  #Run STR associations.
  for(i in 1:length(str_list)){
    if(i%%1000==0){print(i)}
    
    dat.temp <- dat.t.test.s[, c(i, 293752)]
    
    names(dat.temp) <- c("allele_count",'sample') 
    dat.cov.temp <- dat.site_cov.t.test.s[, c(i, 293752,293753)]
    names(dat.cov.temp) <- c("depth",'sample','s_avg_depth')
    dat.temp <- merge(dat.temp, dat.cov.temp, by = "sample")
    
    
    dat.temp <- merge(dat.temp, pheno, by = "sample")
    # Remove rows with missing allele counts
    dat.temp <- dat.temp[!is.na(dat.temp$allele_count), ]
    
    dat.temp$allele_count[is.na(dat.temp$allele_count)] <- 0
    dat.temp$allele_count_inv_norm<-RankNorm(as.numeric(dat.temp$allele_count))   #Inverse normal transform STR genotypes
    
    #Run Linear Model
    if(T){
      tryCatch(
        {
          allele_count_glm = coef(summary(glm(is_case ~ allele_count + is_female + batch + age + depth + s_avg_depth+
                                                PC1 + PC2 + PC3, data = dat.temp, family = "binomial")))
          
          dat.test[i,]$allele_count_logistic_p<-allele_count_glm[2,4]
          dat.test[i,]$allele_count_effect<-allele_count_glm[2,1]
          
          ### mean difference
          dat.test[i,]$allele_count_mean_ad<-mean(as.numeric(subset(dat.temp, is_case==1)$allele_count),na.rm = T)
          dat.test[i,]$allele_count_mean_control<-mean(as.numeric(subset(dat.temp, is_case==0)$allele_count),na.rm = T) 
          dat.test[i,]$allele_count_sd<-sd(as.numeric(dat.temp$allele_count),na.rm = T)  
         
          allele_count_inv_norm_glm = coef(summary(glm(is_case ~ allele_count_inv_norm + is_female + batch + age + depth + s_avg_depth +
                                                         has_E4 + # Comment out if you want to exclude APOE4
                                                         PC1 + PC2 + PC3, data = dat.temp, family = "binomial")))
          
          
          # allele_count_inv_norm의 p-value와 효과 추출 (변수가 존재할 때만)
          if ("allele_count_inv_norm" %in% rownames(allele_count_inv_norm_glm)) {
            dat.test[i,]$allele_count_inv_norm_logistic_p <- allele_count_inv_norm_glm["allele_count_inv_norm", "Pr(>|z|)"]
            dat.test[i,]$allele_count_inv_norm_effect <- allele_count_inv_norm_glm["allele_count_inv_norm", "Estimate"]
          } else {
            dat.test[i,]$allele_count_inv_norm_logistic_p <- NA
            dat.test[i,]$allele_count_inv_norm_effect <- NA
          }
          
        },error=function(e){
          print(paste(i,": error")) # Exclude STRs that cannot be analyzed.
        }
      )
      
    }
  }
  write.table(dat.test, "outputs/eh_RCL40_visual_single_str_test_inv_noPC_240419_pca_rm_APOE_with_orig.txt", row.names=F, col.names = T, sep="\t", quote=F) 
  
}


####### 2.1.  Model 2 Threshold Test  ####### 
{
  
  library(GenomicRanges)
  library(data.table)
  library(stringr)
  
  join_target <- target %>% select(id, ref_length) %>% distinct()
  dat.annot.t.m2 <- dat.annot.t %>% left_join(join_target, by = c('id' = 'id'))
  dat.list = dat.annot.t.m2
  
  # STR list
  str_list = dat.list$id
  
  # Fisher's exact test for each STR length thresholds (by parameter s)
  for(s in c(1,5,10,20)){
    
    print(paste0('Paramter s: ',s))
    
    dat.list = dat.list %>% dplyr::mutate(!!paste0('thr_', s,'_pval') := 1) %>% 
      dplyr::mutate(!!paste0('thr_', s,'_odds') := 0)
    
    
    for(i in 1:length(str_list)){
      if(i%%1000==0){print(i)}
      
      dat.temp <- dat.t.test.s[, c(i, 293752)]
      names(dat.temp) <- c("allele",'sample')
      
      dat.temp$allele_count = dat.temp$allele - dat.list[i,6]
      
      dat.temp <- merge(dat.temp, pheno, by = "sample")
      
      contingency_table <- table(dat.temp$allele_count >= s, dat.temp$is_case)
      # If a 2 by 2 table cannot be formed
      if (dim(contingency_table)[1] != 2 || dim(contingency_table)[2] != 2) {
        next
      }
      
      fisher_result <- fisher.test(contingency_table)
      
      dat.list[[paste0('thr_', s, '_pval')]][i] <- fisher_result$p.value
      dat.list[[paste0('thr_', s, '_odds')]][i] <- fisher_result$estimate
    }
    
    
  }
  
  write.table(dat.list, "outputs/eh_RCL40_visual_single_str_thr_test_ref_240419_PC_rm.txt", row.names=F, col.names = T, sep="\t", quote=F)
  
}

####### 2.2. Model 2 Threshold Test with frequency preparing   ####### 
{
  library(GenomicRanges)
  library(data.table)
  library(stringr)
  
  join_target <- target %>% select(id, ref_length) %>% distinct()
  dat.annot.t.m2 <- dat.annot.t %>% left_join(join_target, by = c('id' = 'id'))
  
  dat <- dat.t.test.s 
  
  dat.list<-dat.annot.t
  dat.list = dat.annot.t.m2
  
  # STR list
  str_list = dat.list$id
  
  #Perform Fisher's exact test at different STR tract length thresholds (parameter s)
  for(s in c(1,5,10,20)){
    
    print(paste0('Paramter s: ',s))
    
    dat.list = dat.list %>% mutate(!!paste0('thr_', s,'_freq') := NA)
    
    for(i in 1:length(str_list)){
      if(i%%1000==0){print(i)}
      
      dat.temp <- dat.t.test.s[, c(i, 293752)]
      names(dat.temp) <- c("allele",'sample')
      
      dat.temp$allele_count = dat.temp$allele - dat.list[i,6]
      
      
      # Calculate thr_s_freq and assign it directly to dat.list
      dat.list[[paste0('thr_', s, '_freq')]][i] <- sum(dat.temp$allele_count >= s, na.rm = T)
      
      
    }
    
    
  }
  
  write.table(dat.list, 'outputs/eh_RCL40_visual_single_str_thr_freq_ref.txt', row.names=F, col.names = T, sep="\t", quote=F) 
  
  
  #With reference
  join_target <- target %>% select(id, ref_length) %>% distinct()
  dat.annot.t.m2 <- dat.annot.t %>% left_join(join_target, by = c('id' = 'id'))
  
  
  dat.list = dat.annot.t.m2
  
  # STR list
  str_list = dat.list$id
  
  pheno.thr = pheno %>% select(sample,is_case)
  
  n <- length(str_list)
  dat.t.test.s.ref <- vector("list", n)
  
  
  for(i in 1:n){
    if(i%%1000==0){print(i)}
    
    dat.temp <- dat.t.test.s[, c(i, 293752)]
    names(dat.temp) <- c("allele_count",'sample')
    
    dat.temp$allele_count <- dat.temp$allele_count - dat.list[i,6]
    dat.temp$id <- rep(dat.list[i,4], nrow(dat.temp))
    dat.temp$index <- rep(dat.list[i,5], nrow(dat.temp))
    
    dat.temp <- merge(dat.temp, pheno.thr, by = "sample")
    
    dat.t.test.s.ref[[i]] <- dat.temp
  }
  
  
  ### data preprocessing
  {
    library(dplyr)
    library(data.table)
    
    
    # Define a function that performs a mutate operation on each element of a list
    mutate_list <- function(df) {
      df <- df %>%
        mutate(!!paste0('thr_', 1,'_freq') := ifelse(allele_count >= 1, 1, 0),
               !!paste0('thr_', 5,'_freq') := ifelse(allele_count >= 5, 1, 0),
               !!paste0('thr_', 10,'_freq') := ifelse(allele_count >= 10, 1, 0),
               !!paste0('thr_', 20,'_freq') := ifelse(allele_count >= 20, 1, 0))
      return(df)
    }
    
    # Create a modified list by applying a mutate operation to each element of the list
    dat.t.test.s.ref_list <- lapply(dat.t.test.s.ref, mutate_list)
    
    
    dat.t.test.s.ref_df <- rbindlist(dat.t.test.s.ref_list)

    data.table::fwrite(dat.t.test.s.ref_df, "eh.v5_RCL40_visual_allele_long_ref_pivot_longer_240419.txt.gz", compress = "gzip")
    dat.t.test.s.ref_df = data.table::fread('eh.v5_RCL40_visual_allele_long_ref_pivot_longer_240419.txt.gz')
    dat.t.test.s.ref_df <- dat.t.test.s.ref_df %>% filter(sample %in% pheno$sample)
    
    # Set up the given dataset (dat.t.test.s.ref) and the required threshold and frequency lists.
    threshold_list <- c(1, 5, 10, 20)
    frequency_list <- c(1, 5, 10, 100, Inf)  # Inf means No limit
    
    dat.thr <- as.data.table(dat.t.test.s.ref_df)
    
    # Create an empty dataframe to store the results
    result_df <- data.frame(threshold = integer(),
                            frequency = integer(),
                            p_value = numeric(),
                            fold_difference = numeric())
    
    dat.list = as.data.table(read_tsv('outputs/eh_RCL40_visual_single_str_thr_freq_ref.txt'))
    
    # Compute for each combination of threshold and frequency.
    for (thr in threshold_list) {
      print(paste0('thr: ',thr))
      for (freq in frequency_list) {
        print(paste0('freq: ',freq))
        if (freq == 1) {
          thr_freq_list = dat.list %>% filter(!!sym(paste0('thr_', thr,'_freq')) == freq) %>% pull(id)
        } else {
          thr_freq_list = dat.list %>% filter(!!sym(paste0('thr_', thr,'_freq')) <= freq) %>% pull(id)
        }
        
        dat.temp <- dat.thr[id %in% thr_freq_list]
        
        AD_group <- dat.temp[is_case == 1, .(thr_freq_count = sum(get(paste0('thr_', thr, '_freq')) > 0, na.rm = TRUE)), by = sample]$thr_freq_count
        ctrl_group <- dat.temp[is_case == 0, .(thr_freq_count = sum(get(paste0('thr_', thr, '_freq')) > 0, na.rm = TRUE)), by = sample]$thr_freq_count
        

        wilcox_test_result <- wilcox.test(AD_group, ctrl_group)
        
        # fold difference 
        mean_AD <- mean(AD_group)  
        mean_control <- mean(ctrl_group)  
        fold_difference <- mean_AD / mean_control
        
        # Add the results to a temporary dataframe
        temp_result <- data.frame(threshold = thr, frequency = freq, p_value = wilcox_test_result$p.value, fold_difference = fold_difference)
        
        # Add the temporary dataframe to the result dataframe
        result_df <- rbind(result_df, temp_result)
      }
    }
    
    
    
    write.table(result_df, 'outputs/eh_RCL40_visual_single_str_thr_freq_ref_test_PC_rm.txt', row.names=F, col.names = T, sep="\t", quote=F) 
    
    
  }
  
  
}

######## 02. STR model 1 length test for clinical diagnosis #######
# Input file prepare
{
  library(readxl)
  library(tidyverse)
  library(GenomicRanges)
  
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
  pheno = pheno %>% select(sample,RCL40_visual, batch ,age,sex, DX,PC1,PC2,PC3, APOE )
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
  
  all_bed %>% dim()
  
  # Genotype info
  dat<-data.table::fread("eh.v5_allele_long_240419.txt.gz", header=T, sep=",", data.table=F) 
  dat %>% dim()
  dat[1:5,1:5]
  
  # STR coverage
  dat.site_cov<-data.table::fread("eh.v5_coverage_240419.txt.gz", header=T, sep=",", data.table=F) 
  dat.site_cov %>% dim()
  
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

########## 0.1. target str ########## 
{
  dat.annot.t = dat.annot %>% filter(id %in% target$id)
  dat.annot.t %>% dim
  target %>% dim
  
  dat.t = dat[dat.annot.t$index,]
  dat.t %>% dim
  dat.t[1:5,1:5]
  
  dat.site_cov.t = dat.site_cov[dat.annot.t$index,]
  dat.site_cov.t %>% dim
  dat.site_cov.t[1:5,1:5]
}

####### 0.2. phenotype & filtering ####### 
{
  sample_list = pheno %>% filter(DX %in% c("DAT","CU")) %>% pull(sample)
  pheno = pheno %>% mutate( is_case = case_when(
    DX == "DAT" ~ 1,
    DX == "CU" ~ 0,
    T ~ NA
  ),
  is_female = case_when(
    sex == 2 ~ 1,
    sex == 1 ~ 0,
    T ~ NA
  ))
  pheno <- pheno %>% filter(!is.na(is_case))
  
  # Remove bad quality sample
  sample_list <- sample_list[sample_list != "WGS_1622"]
  
  dat.t.test = dat.t %>% select(all_of(sample_list))
  dat.site_cov.t.test = dat.site_cov.t %>% select(all_of(sample_list))
  
}

####### 0.3. make test file ####### 
{
  pheno$is_female <- factor(pheno$is_female)
  pheno$batch <- factor(pheno$batch)
  
  
  dat.t.test.1 = t(dat.t.test)
  dat.t.test.s = dat.t.test.1 %>% as_tibble() %>% mutate('sample' = rownames(dat.t.test.1))
  
  dat.site_cov.t.test.1 = t(dat.site_cov.t.test)
  
  dat.site_cov.t.test.s = dat.site_cov.t.test.1 %>% as_tibble() %>% mutate('sample' = rownames(dat.site_cov.t.test.1))
  dat.site_cov.t.test.s = dat.site_cov.t.test.s %>% left_join(dat.sample_cov)
  
  rm(dat,dat.annot,dat.site_cov,dat.site_cov.t,
     dat.site_cov.t.test,dat.site_cov.t.test.1,
     dat.t,dat.t.test,dat.t.test.1)
}

####### 1. Model 1 Length Test ####### 
{
  ###Run single STR association analysis
  library(data.table)
  library(stringr)
  library(RNOmni)
  library(jsonlite)
  
  # Create APOE4 columns
  pheno <- pheno %>%
    mutate(
      E4_count = str_count(APOE, "E4"),  # Count the number of E4 alleles (0, 1, 2)
      has_E4 = ifelse(E4_count > 0, 1, 0)  # 1 if E4 is present, else 0
    )
  
  dat.annot.t %>% head()
  dat.test = dat.annot.t 
  
  str_list = dat.test$index
  
  dat.test$allele_count_logistic_p<-1
  dat.test$allele_count_effect<-0
  dat.test$allele_count_mean_control<-0
  dat.test$allele_count_mean_ad<-0
  dat.test$allele_count_sd<-0
  
  dat.test$allele_count_inv_norm_logistic_p<-1
  dat.test$allele_count_inv_norm_effect<-0
  
  
  pheno$is_female <- factor(pheno$is_female)
  pheno$batch <- factor(pheno$batch)
  
  
  #Run STR associations.
  for(i in 1:length(str_list)){
    if(i%%1000==0){print(i)}
    
    dat.temp <- dat.t.test.s[, c(i, 293752)]
    
    names(dat.temp) <- c("allele_count",'sample') 
    dat.cov.temp <- dat.site_cov.t.test.s[, c(i, 293752,293753)]
    names(dat.cov.temp) <- c("depth",'sample','s_avg_depth')
    dat.temp <- merge(dat.temp, dat.cov.temp, by = "sample")
    
    
    dat.temp <- merge(dat.temp, pheno, by = "sample")
    # Remove rows with missing allele counts
    dat.temp <- dat.temp[!is.na(dat.temp$allele_count), ]
    
    dat.temp$allele_count[is.na(dat.temp$allele_count)] <- 0
    dat.temp$allele_count_inv_norm<-RankNorm(as.numeric(dat.temp$allele_count))   #Inverse normal transform STR genotypes
    
    #Run Linear Model
    if(T){
      tryCatch(
        {
          # allele_count_glm = coef(summary(glm(is_case ~ allele_count + is_female + batch + age + depth + s_avg_depth+
          #                                       PC1 + PC2 + PC3, data = dat.temp, family = "binomial")))
          # 
          # dat.test[i,]$allele_count_logistic_p<-allele_count_glm[2,4]
          # dat.test[i,]$allele_count_effect<-allele_count_glm[2,1]
          
          ### mean difference
          dat.test[i,]$allele_count_mean_ad<-mean(as.numeric(subset(dat.temp, is_case==1)$allele_count),na.rm = T)
          dat.test[i,]$allele_count_mean_control<-mean(as.numeric(subset(dat.temp, is_case==0)$allele_count),na.rm = T) 
          dat.test[i,]$allele_count_sd<-sd(as.numeric(dat.temp$allele_count),na.rm = T)  
          
          allele_count_inv_norm_glm = coef(summary(glm(is_case ~ allele_count_inv_norm + is_female + batch + age + depth + s_avg_depth +
                                                         has_E4 + # Comment out if you want to exclude APOE4
                                                         PC1 + PC2 + PC3, data = dat.temp, family = "binomial")))
          
          
          # allele_count_inv_norm의 p-value와 효과 추출 (변수가 존재할 때만)
          if ("allele_count_inv_norm" %in% rownames(allele_count_inv_norm_glm)) {
            dat.test[i,]$allele_count_inv_norm_logistic_p <- allele_count_inv_norm_glm["allele_count_inv_norm", "Pr(>|z|)"]
            dat.test[i,]$allele_count_inv_norm_effect <- allele_count_inv_norm_glm["allele_count_inv_norm", "Estimate"]
          } else {
            dat.test[i,]$allele_count_inv_norm_logistic_p <- NA
            dat.test[i,]$allele_count_inv_norm_effect <- NA
          }
          
        },error=function(e){
          print(paste(i,": error")) # Exclude STRs that cannot be analyzed.
        }
      )
      
    }
  }
  write.table(dat.test, "outputs/eh_DX_single_str_test_inv_noPC_240419_pca_rm_APOE.txt", row.names=F, col.names = T, sep="\t", quote=F) 
  
}

####### 2.1. Model 2 Threshold Test  ####### 
{
  
  library(GenomicRanges)
  library(data.table)
  library(stringr)
  
  join_target <- target %>% select(id, ref_length) %>% distinct()
  dat.annot.t.m2 <- dat.annot.t %>% left_join(join_target, by = c('id' = 'id'))
  dat.list = dat.annot.t.m2
  
  # STR list
  str_list = dat.list$id
  
  # Fisher's exact test for each STR length thresholds (by parameter s)
  for(s in c(1,5,10,20)){
    
    print(paste0('Paramter s: ',s))
    
    dat.list = dat.list %>% dplyr::mutate(!!paste0('thr_', s,'_pval') := 1) %>% 
      dplyr::mutate(!!paste0('thr_', s,'_odds') := 0)
    
    for(i in 1:length(str_list)){
      if(i%%1000==0){print(i)}

      dat.temp <- dat.t.test.s[, c(i, 293752)]
      names(dat.temp) <- c("allele",'sample')
      
      dat.temp$allele_count = dat.temp$allele - dat.list[i,6]
      
      
      dat.temp <- merge(dat.temp, pheno, by = "sample")
      
      contingency_table <- table(dat.temp$allele_count >= s, dat.temp$is_case)
      
      # If a 2 by 2 table cannot be formed
      if (dim(contingency_table)[1] != 2 || dim(contingency_table)[2] != 2) {
        next
      }
   
      fisher_result <- fisher.test(contingency_table)
      
      dat.list[[paste0('thr_', s, '_pval')]][i] <- fisher_result$p.value
      dat.list[[paste0('thr_', s, '_odds')]][i] <- fisher_result$estimate
    }
    
    
  }
  
  write.table(dat.list, "outputs/eh_DX_single_str_thr_test_ref_240419_PC_rm.txt", row.names=F, col.names = T, sep="\t", quote=F) 
  
}

####### 2.2. Model 2 Threshold Test with frequency preparing ####### 
{
  library(GenomicRanges)
  library(data.table)
  library(stringr)
  
  join_target <- target %>% select(id, ref_length) %>% distinct()
  dat.annot.t.m2 <- dat.annot.t %>% left_join(join_target, by = c('id' = 'id'))
  
  dat <- dat.t.test.s 
  
  dat.list<-dat.annot.t
  
  dat.list = dat.annot.t.m2
  
  # STR list
  str_list = dat.list$id
  
  #Perform Fisher's exact test at different STR tract length thresholds (parameter s)
  for(s in c(1,5,10,20)){
    
    print(paste0('Paramter s: ',s))
    
    dat.list = dat.list %>% mutate(!!paste0('thr_', s,'_freq') := NA)
    
    for(i in 1:length(str_list)){
      if(i%%1000==0){print(i)}
      
      dat.temp <- dat.t.test.s[, c(i, 293752)]
      names(dat.temp) <- c("allele",'sample')
      
      dat.temp$allele_count = dat.temp$allele - dat.list[i,6]
      
      # Calculate thr_s_freq and assign it directly to dat.list
      dat.list[[paste0('thr_', s, '_freq')]][i] <- sum(dat.temp$allele_count >= s, na.rm = T)
      
    }
    
    
  }
  
  write.table(dat.list, 'outputs/eh_DX_single_str_thr_freq_ref_240419_PC_rm.txt', row.names=F, col.names = T, sep="\t", quote=F) #Write output file
  

  #With reference
  join_target <- target %>% select(id, ref_length) %>% distinct()
  dat.annot.t.m2 <- dat.annot.t %>% left_join(join_target, by = c('id' = 'id'))
  
  dat.list = dat.annot.t.m2
  
  # STR list
  str_list = dat.list$id
  
  pheno.thr = pheno %>% select(sample,is_case)
  
  n <- length(str_list)
  dat.t.test.s.ref <- vector("list", n)
  
  
  for(i in 1:n){
    if(i%%1000==0){print(i)}
    
    dat.temp <- dat.t.test.s[, c(i, 293752)]
    names(dat.temp) <- c("allele_count",'sample')
    
    dat.temp$allele_count <- dat.temp$allele_count - dat.list[i,6]
    dat.temp$id <- rep(dat.list[i,4], nrow(dat.temp))
    dat.temp$index <- rep(dat.list[i,5], nrow(dat.temp))
    
    dat.temp <- merge(dat.temp, pheno.thr, by = "sample")
    
    dat.t.test.s.ref[[i]] <- dat.temp
  }
  
  
  ### data preprocessing
  {
    library(dplyr)
    library(data.table)
    
    
    # Define a function that performs a mutate operation on each element of a list
    mutate_list <- function(df) {
      df <- df %>%
        mutate(!!paste0('thr_', 1,'_freq') := ifelse(allele_count >= 1, 1, 0),
               !!paste0('thr_', 5,'_freq') := ifelse(allele_count >= 5, 1, 0),
               !!paste0('thr_', 10,'_freq') := ifelse(allele_count >= 10, 1, 0),
               !!paste0('thr_', 20,'_freq') := ifelse(allele_count >= 20, 1, 0))
      return(df)
    }
    
    # Create a modified list by applying a mutate operation to each element of the list
    dat.t.test.s.ref_list <- lapply(dat.t.test.s.ref, mutate_list)
    
    
    dat.t.test.s.ref_df <- rbindlist(dat.t.test.s.ref_list)
  }
    
    data.table::fwrite(dat.t.test.s.ref_df, "eh.v5_DX_allele_long_ref_pivot_longer_240419.txt.gz", compress = "gzip")
    dat.t.test.s.ref_df = data.table::fread('eh_DX_allele_long_ref_pivot_longer_240419.txt.gz')
    dat.t.test.s.ref_df <- dat.t.test.s.ref_df %>% filter(sample %in% pheno$sample)
    
    # Set up the given dataset (dat.t.test.s.ref) and the required threshold and frequency lists.
    threshold_list <- c(1, 5, 10, 20)
    frequency_list <- c(1, 5, 10, 100, Inf)  # Inf means No limit
    
    
    dat.thr <- as.data.table(dat.t.test.s.ref_df)
    
    # Create an empty dataframe to store the results
    result_df <- data.frame(threshold = integer(),
                            frequency = integer(),
                            p_value = numeric(),
                            fold_difference = numeric())
    
    dat.list = as.data.table(read_tsv('outputs/eh_DX_single_str_thr_freq_ref.txt'))
    
    # Compute for each combination of threshold and frequency.
    for (thr in threshold_list) {
      print(paste0('thr: ',thr))
      for (freq in frequency_list) {
        print(paste0('freq: ',freq))
        if (freq == 1) {
          thr_freq_list = dat.list %>% filter(!!sym(paste0('thr_', thr,'_freq')) == freq) %>% pull(id)
        } else {
          thr_freq_list = dat.list %>% filter(!!sym(paste0('thr_', thr,'_freq')) <= freq) %>% pull(id)
        }
        
        dat.temp <- dat.thr[id %in% thr_freq_list]
        
        AD_group <- dat.temp[is_case == 1, .(thr_freq_count = sum(get(paste0('thr_', thr, '_freq')) > 0, na.rm = TRUE)), by = sample]$thr_freq_count
        ctrl_group <- dat.temp[is_case == 0, .(thr_freq_count = sum(get(paste0('thr_', thr, '_freq')) > 0, na.rm = TRUE)), by = sample]$thr_freq_count
        
        wilcox_test_result <- wilcox.test(AD_group, ctrl_group)
        
        # fold difference 
        mean_AD <- mean(AD_group) 
        mean_control <- mean(ctrl_group)  
        fold_difference <- mean_AD / mean_control
        
        # Add the results to a temporary dataframe
        temp_result <- data.frame(threshold = thr, frequency = freq, p_value = wilcox_test_result$p.value, fold_difference = fold_difference)
        
        # Add the temporary dataframe to the result dataframe
        result_df <- rbind(result_df, temp_result)
      }
    }
    
    write.table(result_df, 'outputs/eh_DX_single_str_thr_freq_ref_test_PC_rm.txt', row.names=F, col.names = T, sep="\t", quote=F)
}



