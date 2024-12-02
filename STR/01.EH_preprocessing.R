#############################################################
## Short tandem repeat data analysis
## ExpansionHunter output file preprocessing
##
## This analysis includes the following steps:
## 01. Get STR genotype
## 02. Get STR coverage
## 03. Genotype merge
## 04. Longer genotype merge
## 05. Coverage merge
#############################################################


######## 01. Get STR genotype ######## 

library(stringr)
setwd("~/Dropbox/ADWGS/STR/EH/")
sample_list<-list.files("vcf/")

dat.bed<-read.delim("eh.v5_w_gangstr.v13.polymorphic_240419.bed", header=T, sep="\t", stringsAsFactors = F)
dat.bed %>% head()
dat.bed$pos_start <- as.integer(dat.bed$pos_start)
dat.bed$index <- as.integer(dat.bed$index)

for(s in sample_list){
  
  # genotypes from VCF
  filename<-paste0("vcf/",s)
  dat.temp<-read.delim(filename, comment.char = "#", header=F, sep="\t")
  names(dat.temp)[c(1:2)]<-c("chr", "pos_start")
  dat.temp$pos_start = as.numeric(dat.temp$pos_start)
  dat.temp$RU <- str_extract(str_split_fixed(dat.temp$V8, "\\;", 8)[,4], "(?<=RU=)[^;]+")
  dat.temp$gt<-str_split_fixed(dat.temp$V10, "\\:", 8)[,3]
  dat.temp<-subset(dat.temp, select=c(chr, pos_start,RU, gt))
  dat.temp$chr_num <- gsub("chr", "", dat.temp$chr)
  dat.temp$chr_num <- as.integer( dat.temp$chr_num)
  
  # For STRs with multiple entries, select the one containing the longest allele
  dat.temp$long<-as.numeric(str_split_fixed(dat.temp$gt, "\\/", 2)[,2])
  dat.temp<-dat.temp[order(dat.temp$chr_num,dat.temp$pos_start,dat.temp$RU, -dat.temp$long),]
  dat.temp$id <- paste(dat.temp$chr, dat.temp$pos_start,dat.temp$RU, sep = '-')
  dat.temp<-dat.temp[!duplicated(dat.temp$id),]
  dat.vcf<-dat.temp
  
  # Merge the genotypes with the BED file and then reorder them
  dat.bed_tmp<-merge(dat.bed, dat.vcf, by=c("id"), all.x=T, all.y=F)
  dat.bed_tmp<-dat.bed_tmp[order(dat.bed_tmp$index),]
  
  # output data with genotype
  dat.gt_tmp<-subset(dat.bed_tmp, select=c(gt))
  dat.long_tmp<-subset(dat.bed_tmp, select=c(long))
  
  pattern <- "WGS_\\d+"
  s <- str_extract(s, pattern)
  
  names(dat.gt_tmp)<-s
  names(dat.long_tmp)<-s
  
  write.table(dat.gt_tmp, paste0( 'gt_240419/',s, ".eh-v5.genotype.txt"), row.names = F, col.names = T, sep="\t", quote=F)
  write.table(dat.long_tmp, paste0( 'gt_long_240419/',s, ".eh-v5.genotype_long.txt"), row.names = F, col.names = T, sep="\t", quote=F)
  
}



######## 02. Get STR coverage ########


library(stringr)
library("rjson")
setwd("~/Dropbox/ADWGS/STR/EH/")

dat.bed<-read.delim("eh.v5_w_gangstr.v13.polymorphic_240419.bed", header=T, sep="\t", stringsAsFactors = F)
dat.bed %>% head()
dat.bed$pos_start <- as.integer(dat.bed$pos_start)
dat.bed$index <- as.integer(dat.bed$index)


sample_list<-list.files("vcf")

for(s in sample_list){
  
  # Read in vcf and obtain coverage
  filename<-paste0("vcf/",s)
  dat.temp<-read.delim(filename, comment.char = "#", header=F, sep="\t")
  names(dat.temp)[c(1:2)]<-c("chr", "pos_start")
  dat.temp$pos_start = as.numeric(dat.temp$pos_start)
  dat.temp$RU <- str_extract(str_split_fixed(dat.temp$V8, "\\;", 8)[,4], "(?<=RU=)[^;]+")
  dat.temp$gt<-str_split_fixed(dat.temp$V10, "\\:", 8)[,3]
  dat.temp$coverage<-str_split_fixed(dat.temp$V10, "\\:", 8)[,8]
  dat.temp<-subset(dat.temp, select=c(chr, pos_start, RU,gt, coverage))
  dat.temp$chr_num <- gsub("chr", "", dat.temp$chr)
  dat.temp$chr_num <- as.integer( dat.temp$chr_num)
  
  # For STRs with multiple entries, select the one containing the longest allele with coverage
  dat.temp$long<-as.numeric(str_split_fixed(dat.temp$gt, "\\/", 2)[,2])
  dat.temp<-dat.temp[order(dat.temp$chr_num,dat.temp$pos_start,dat.temp$RU, -dat.temp$long),]
  dat.temp$id <- paste(dat.temp$chr, dat.temp$pos_start,dat.temp$RU, sep = '-')
  dat.temp<-dat.temp[!duplicated(dat.temp$id),]
  dat.vcf<-dat.temp
  
  # Merge the genotypes and coverages with the BED file and then reorder them
  dat.bed_tmp<-merge(dat.bed, dat.vcf, by=c("id"), all.x=T, all.y=F)
  dat.bed_tmp<-dat.bed_tmp[order(dat.bed_tmp$index),]
  
  # select coverage
  dat.bed_tmp<-subset(dat.bed_tmp, select=c(coverage))
  
  pattern <- "WGS_\\d+"
  s <- str_extract(s, pattern)
  
  names(dat.bed_tmp)<-s
  write.table(dat.bed_tmp, paste0( 'coverage_240419/',s, ".eh-v5.coverage.txt"), row.names = F, col.names = T, sep="\t", quote=F)
}



######## 03. Genotype merge  ########
# File directory
folder_path <- "~/Dropbox/ADWGS/STR/EH/gt_240419/"

# Get file list from folder
file_list <- list.files(folder_path, full.names = TRUE)
unique_names <- gsub(".*/|\\.eh-v5\\.genotype\\.txt", "", file_list)

# Combine files
combined_data <- do.call(cbind, sapply(file_list,data.table::fread, simplify = FALSE))
combined_data %>% dim

names(combined_data) <- unique_names
combined_data[1:5,1:5]

eh.gt.table = combined_data %>% as_tibble()

data.table::fwrite(eh.gt.table, "~/Dropbox/ADWGS/STR/EH/eh.v5_allele_240419.txt.gz", compress = "gzip")

######## 04. Longer genotype merge  ########
folder_path <- "~/Dropbox/ADWGS/STR/EH/gt_long_240419/"

file_list <- list.files(folder_path, full.names = TRUE)
unique_names <- gsub(".*/|\\.eh-v5\\.genotype_long\\.txt", "", file_list)

# Combine files
combined_data <- do.call(cbind, sapply(file_list,data.table::fread, simplify = FALSE))
combined_data %>% dim

names(combined_data) <- unique_names
combined_data[1:5,1:5]

eh.gt.table = combined_data %>% as_tibble()

data.table::fwrite(eh.gt.table, "~/Dropbox/ADWGS/STR/EH/eh.v5_allele_long_240419.txt.gz", compress = "gzip")

######## 05. Coverage merge ######## 
folder_path <- "~/Dropbox/ADWGS/STR/EH/coverage_240419/"

file_list <- list.files(folder_path, full.names = TRUE)
unique_names <- gsub(".*/|\\.eh-v5\\.coverage\\.txt", "", file_list)

# Combine files
combined_data <- do.call(cbind, sapply(file_list,data.table::fread, simplify = FALSE))

names(combined_data) <- unique_names
combined_data[1:5,1:5]

data.table::fwrite(combined_data, "~/Dropbox/ADWGS/STR/EH/eh.v5_coverage_240419.txt.gz", compress = "gzip")

eh.coverage.table = data.table::fread("~/Dropbox/ADWGS/STR/EH/eh.v5_coverage.txt.gz")




