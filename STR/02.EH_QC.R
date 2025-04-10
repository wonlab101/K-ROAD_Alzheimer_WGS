#############################################################
## Short tandem repeat data analysis
## ExpansionHunter output file preprocessing
##
## This analysis includes the following steps:
## 01. STR coverage PCA
## 02. Coverage comparision
#############################################################

######## Load files #######
# Input file prepare
{
  library(readxl)
  library(tidyverse)
  library(GenomicRanges)
  library(RNOmni)
  
  rm(list=ls())
  sessionInfo()
    
  rename = dplyr::rename
  select = dplyr::select
  mutate = dplyr::mutate
  
  
  # Phenotype data
  pheno = read_tsv("phenotype_table.tsv") %>% dplyr::rename(sample = IID)
  pheno = pheno %>% select(sample,RCL40_visual, batch ,age,sex, DX,PC1,PC2,PC3,APOE )
  pheno = pheno %>% filter(sample != 'WGS_1622')
  
  # Target chormosomes
  chr_list = paste0("chr",1:22)
  
  # bed file
  dat.annot<-read.delim("eh.v5_w_gangstr.v13.polymorphic.bed", sep="\t", stringsAsFactors = F) #read in bed file of STRs
  dat.annot = dat.annot %>% filter(chr %in% chr_list)
  dat.annot %>% dim()
  
  # BED file from Guo et al. 2023 supplementary table 
  # https://www.medrxiv.org/content/10.1101/2023.11.16.23298623v1
  all_bed = readxl::read_xlsx('supp_1.xlsx', sheet = 1) %>% rename('start' = "STR start position (GRCh38)", 'end' =  "STR end position (GRCh38)", 'ref_length' = "Reference tract length" ) %>% mutate(id = paste(chr ,   start  , motif,sep = '-'))
  
  # Genotype info
  dat<-data.table::fread("eh.v5_allele_long.txt.gz", header=T, sep=",", data.table=F) 
  
  # STR coverage
  dat.site_cov<-data.table::fread("eh.v5_coverage.txt.gz", header=T, sep=",", data.table=F) 
  
  # Sample mean coverage
  dat.sample_cov = data.table::fread("eh.v5.sample_coverage.txt.gz")
  
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
  pheno = pheno %>% mutate( is_female = case_when(
    sex == 2 ~ 1,
    sex == 1 ~ 0,
    T ~ NA
  ))

  pheno$is_female <- factor(pheno$is_female)
  pheno$batch <- factor(pheno$batch)
}


######## 01. STR coverage PCA ######## 

eh.coverage.table = dat.site_cov.t

pcsToUse = 1:2

eh.coverage.table_NA_rm <- na.omit(eh.coverage.table)
eh.coverage.table_NA_rm[1:5,1:5] 
eh.coverage.table_NA_rm = eh.coverage.table_NA_rm
eh.coverage.table_NA_rm = eh.coverage.table_NA_rm[paste(c(ratio_1)),] 
eh.coverage.table_NA_rm %>% dim()

pca <- prcomp(t(eh.coverage.table_NA_rm)
              # ,center = T
              ,scale. = T
)

pcs <- paste0("PC", pcsToUse)
d <- data.frame(V1=pca$x[,pcsToUse[1]],
                V2=pca$x[,pcsToUse[2]], 
                sample = rownames(pca$x)
)
d<- d %>% left_join(pheno, by = c('sample' = 'sample'))
# colnames(d)[1:2] <- pcs

pca$x[1:5,1:5]

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
percentVar <- round(100 * percentVar)

d$batch = as.character(d$batch)

sd_pc1 <- sd(d$V1)
sd_pc2 <- sd(d$V2)

ggplot(d, aes(V1, V2, color = batch)) +
  geom_point(size=3,alpha = 0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle('PCA from EH calling coverage') +
  scale_shape_manual(values = c(17))+
  # theme_step1() +
  geom_vline(xintercept = 4*sd_pc1, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 4*sd_pc2, color = "red", linetype = "dashed") +
  geom_vline(xintercept = -4*sd_pc1, color = "red", linetype = "dashed") +
  geom_hline(yintercept = -4*sd_pc2, color = "red", linetype = "dashed") +
  geom_hline(yintercept = -220, color = "grey", linetype = "dashed")

bad = d %>% filter(V1 > 4*sd_pc1 | V1 < -4*sd_pc1 |  V2 > 4*sd_pc2|  V2 < -220) %>% pull(sample)

d = d %>% mutate(pc_remove = ifelse(sample %in% bad, 1, 0))

write_tsv(d, 'eh.v5_coverage_1558_PCA.txt')


######## 02. Coverage comparision ######## 

eh.coverage.mean = data.table::fread( "eh.v5.sample_coverage.txt.gz")
eh.coverage.mean = eh.coverage.mean %>% filter(!sample %in%pca_rm_samples)

eh.coverage.mean %>% dim()
pheno = read_tsv("phenotype_table.tsv") %>% dplyr::rename(sample = IID)
eh.coverage.mean = eh.coverage.mean %>% left_join(pheno)

eh.coverage.mean$RCL40_visual <- as.factor(eh.coverage.mean$RCL40_visual)

eh.coverage.mean$eh.coverage.mean = as.numeric(eh.coverage.mean$eh.coverage.mean)

eh.coverage.mean$batch = as.factor(eh.coverage.mean$batch)

# Amyloid beta 
d = glm(eh.coverage.mean ~ RCL40_visual + batch, data=eh.coverage.mean)
summary(d)

# Clinical diagnosis
d = glm(eh.coverage.mean ~ DX + batch, data=subset(eh.coverage.mean, DX != 'MCI'))
summary(d)
  
  
