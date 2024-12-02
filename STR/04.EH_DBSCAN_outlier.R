#############################################################
## Short tandem repeat data analysis
## ExpansionHunter output file preprocessing
##
## This analysis includes the following steps:
## 01. Outlier calling by DBSCAN
#############################################################

########## 01. Outlier calling by DBSCAN ########## 
library(data.table)
library(stringr)
library(RNOmni)
library(dbscan)
# set test files

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
  
  
  
}

#Set default dbscan parameters
minpts<-2
eps<-2


#Read in test annotation file
join_target <- target %>% select(id, ref_length) %>% distinct()
dat.annot <- dat.annot.t %>% left_join(join_target, by = c('id' = 'id'))  #read in bed file of STRs

#Read in merged STR genotypes 
dat<-dat.t.test %>% mutate(index = rownames(dat.t.test)) %>% mutate(index = as.integer(index))
# dat[1:5,1550:1558]
dat<-left_join(dat.annot,dat, by = c('index' = 'index')) #여기 살짝쿵 애매
dat[1:5,1:10]
dat %>% dim()
dat[1:5,1515:1521]

#read in coverage file
dat.site_cov<-dat.site_cov.t.test  %>% mutate(index = rownames(dat.site_cov.t)) %>% mutate(index = as.integer(index))
dat.site_cov[1:5,1:10]
# dat.site_cov<-cbind(dat.annot, dat.site_cov)
dat.site_cov<-left_join(dat.annot,dat.site_cov, by = c('index' = 'index'))
dat.site_cov %>% dim()
dat.site_cov[1:5,1:10]


#Read in phenotype file
dat.pheno<-pheno
dat.pheno<-subset(dat.pheno, sample%in%names(dat))
dat.pheno<- dat.pheno %>% left_join(dat.sample_cov)



#Generate case and control sample lists
case_list<-subset(dat.pheno, is_case==1 & sample%in%names(dat))$sample
control_list<-subset(dat.pheno, is_case==0 & sample%in%names(dat))$sample
dat<-dat[,c(1:6, which(names(dat)%in%c(case_list, control_list)))]


#Process STR site coverage file
dat.site_cov<-dat.site_cov[,c(names(dat))]
dat.site_cov[1:10,1:10]

# dat.site_cov$id<-paste0(dat.site_cov$chr, "_", dat.site_cov$pos_start,"_",dat.site_cov$pos_end)

#Helper function to calculate mode of STR tract lengths
Mode <- function(x){ 
  a = table(as.numeric(x)) # x is a vector
  return(as.numeric(names(a[which.max(a)])))
}

#Helper function to calculate number of non-reference STR alleles
non_ref_count<-function(temp_vector){ 
  temp_vector<-as.numeric(temp_vector)
  length(temp_vector[!is.na(temp_vector) & temp_vector!=Mode(temp_vector)])
}

#Make output file
dat.test<-dat[,c(1:6)]
dat.test$num_allele<-0
dat.test$alleles_is_case<-NA
dat.test$alleles_control<-NA
dat.test$outliers<-NA
dat.test$outlier_size<-NA

# setDT(dat)
# setDT(dat.site_cov)
# setDT(pheno)

dat = dat %>% arrange(index)
dat.site_cov = dat.site_cov %>% arrange(index)
dat%>% dim
dat.site_cov %>% dim

error_list = c()

#Run DBSCAN for each STR
for(i in 1:nrow(dat.test)){
  if(i%%1000==0){print(i)}
  # i = which(str_list == 9032) 
  
  
  dat.temp<-data.frame(t(dat[i,c(7:ncol(dat))]), stringsAsFactors = F) #Generate temporary dataframe of genotypes for STR
  str<-paste0(dat[i,c(1:3)], collapse="-")
  names(dat.temp)[1]<-"str"
  dat.temp$sample<-rownames(dat.temp)
  
  dat.site_cov_temp<-dat.site_cov[i,]
  dat.site_cov_temp<-data.frame(t(dat.site_cov[i,c(7:ncol(dat.site_cov_temp))]), stringsAsFactors = F)
  dat.temp$site_cov<-dat.site_cov_temp[,1]  
  
  dat.temp<-merge(dat.temp, dat.pheno, by="sample", all.x=F, all.y=F)
  
  dat.temp$str<-as.numeric(dat.temp$str)
  dat.temp <- dat.temp[!is.na(dat.temp$str), ]
  
  non_ref<-non_ref_count(dat.temp$str) #calculate number of individuals with non-reference tract lengths
  dat.test[i,]$num_allele<-length(unique(dat.temp$str)) #Calculate number of unique str alleles
  
  
  if(non_ref > 0 & nrow(dat.temp) > 0){
    tryCatch({
      is_case_table <- data.frame(table(subset(dat.temp, is_case == 1)$str)) 
      dat.test[i,]$alleles_is_case <- paste(paste(is_case_table[,1], is_case_table[,2], sep = ","), collapse = "|") #Generate list of case STR alleles
      
      control_table <- data.frame(table(subset(dat.temp, is_case == 0)$str)) #Generate list of control STR alleles
      dat.test[i,]$alleles_control <- paste(paste(control_table[,1], control_table[,2], sep = ","), collapse = "|") 
      
      dat.temp$residuals <- residuals(lm(str ~ is_female + eh.coverage.mean + batch + site_cov + PC1 + PC2 + PC3, data = dat.temp)) #Calculate residuals of STR tract lengths after correcting for covariants
      
      ref <- as.numeric(names(which.max(table(dat.temp$str)))) #Set reference STR tract length for dbscan as mode tract length in the cohort
      range <- max(as.numeric(eps) * ref, quantile(dat.temp$residuals, 0.95, na.rm = TRUE) - quantile(dat.temp$residuals, 0.05, na.rm = TRUE)) #Calculate range for dbscan as done in Trost et al
      
      scan <- dbscan::dbscan(matrix(dat.temp$residuals), eps = range, minPts = ceiling(log2(minpts * nrow(dat.temp)))) #Run dbscan
      
      # Calculate number of clusters in dbscan and identify cutoffs for outliers
      if(length(unique(scan$cluster)) == 1 | sum(scan$cluster == 0) == 0){
        cutoff <- Inf
      } else {
        cutoff <- max(dat.temp[scan$cluster != 0,]$residuals)
        cutoff <- ifelse(cutoff < 2, 2, cutoff)
      }
      
      # Output outliers from dbscan
      if(nrow(subset(dat.temp, residuals > cutoff)) > 0){
        dat.test[i,]$outliers <- paste(dat.temp[dat.temp$residuals > cutoff,]$sample, collapse = ";") 
        dat.test[i,]$outlier_size <- paste(dat.temp[dat.temp$residuals > cutoff,]$str, collapse = ";")
      }
    }, error = function(e) {
      print(paste0('error: ',dat.test[i,4]))
    })
  }
  
}

write_tsv(dat.test,"~/Dropbox/ADWGS/STR/EH/outputs/eh_DBSCAN_outlier.tsv")

dat.dbscan<-read.delim(paste0("~/Dropbox/ADWGS/STR/EH/outputs/eh_DBSCAN_outlier.tsv"), stringsAsFactors = F, sep="\t", header=T)
dat.dbscan<-subset(dat.dbscan, !is.na(outliers))

#Summaryize dbscan output
dat.out<-data.frame(chr=character(), pos_start=numeric(), pos_end=numeric(), motif=character(), sample=character(), outlier_length=numeric(), stringsAsFactors = F)
for(i in c(1:nrow(dat.dbscan))){
  if(i%%1000==0){print(i)}
  dat.temp<-data.frame(chr=dat.dbscan[i,]$chr, pos_start=dat.dbscan[i,]$pos_start, 
                       # pos_end=dat.dbscan[i,]$pos_end, 
                       motif=dat.dbscan[i,]$RU,
                       ref=dat.dbscan[i,]$ref_length,
                       sample=as.vector(str_split_fixed(dat.dbscan[i,]$outliers, "\\;", Inf)),outlier_length=as.vector(as.numeric(str_split_fixed(dat.dbscan[i,]$outlier_size, "\\;", Inf))), 
                       stringsAsFactors = F)
  dat.out<-rbind(dat.temp, dat.out)
}
dat.out$sample %>% unique() %>% length

dat.out.save = dat.out %>% left_join(pheno)

write_tsv(dat.out.save, "outputs/eh_DBSCAN_outlier.output.tsv")
