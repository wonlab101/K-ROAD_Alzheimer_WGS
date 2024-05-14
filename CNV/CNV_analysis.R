#library
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(stringr)
library(logistf)
library(readxl)
library(gridExtra)
library(ggmanh)

setwd("/Users/idea_macbookpro/Desktop/CNV/batch1-7/AD_MCI_NC/")
del1<-read.delim("del1_filtering.tsv", sep="\t", header=T, na.strings = c("", " ", NA), stringsAsFactors = FALSE, row.names=NULL)
del2<-read.delim("del2_filtering.tsv", sep="\t", header=T, na.strings = c("", " ", NA), stringsAsFactors = FALSE, row.names=NULL)
dup<-read.delim("dup_filtering.tsv", sep="\t", header=T, na.strings = c("", " ", NA), stringsAsFactors = FALSE, row.names=NULL)
#ins<-read.delim("ins_filtering.tsv", sep="\t", header=T, na.strings = c("", " ", NA), stringsAsFactors = FALSE, row.names=NULL)
data<-rbind(del1,del2,dup)
data<-distinct(data, AnnotSV_ID, .keep_all = T)
master<-read.csv("WGS_master_sheet.csv", na.strings = c(""," ", NA), stringsAsFactor = FALSE)
master<-master %>% rename("status"="DX")
master_final <- master %>% select(FID, sex, age, status, PC1, PC2, PC3, RCL40)
master_final<-arrange(master_final, FID)
#apoe<-read.delim("apoe", sep="", header=T, na.strings = c("", " ", NA), stringsAsFactors = FALSE, row.names=NULL)
#master_final<-merge(master_final, apoe, by.x="FID", by.y="ID")
sample_number=nrow(master_final)

test<-data[,c(15:(sample_number+14))]
test_t=t(test)
a<-substr(test_t, 1,3)
b<-gsub("1/1", "2", a)
c<-gsub("0/1", "1", b)
d<-gsub("./.", "0", c)
d<-as.data.frame(d)
info<-data %>% select(SV_type,CytoBand, SV_chrom, SV_start, SV_end, SV_length, AnnotSV_ID)
save<-d
#DAT vs CN 용
set<-master_final %>% select(FID, status)
d<-cbind(set,d)
d<-d %>% filter(status=="DAT" | status=="CU")
d<-d[,-1:-2]
master_final<-master_final %>% filter(status=="DAT" | status=="CU")
sample_number=nrow(master_final)
gene<-data.frame(data$Gene_name)
gene<-as.vector(gene[,1])
names(d)<-gene
for(i in 1:ncol(d)){
  d[,i] = as.numeric(d[,i])
}
#for(i in 1:ncol(d)){
#  d[,i] = d[,i]-1
#}

for(i in 1:ncol(d)){
  d[(sample_number + 1),i]<-sum(d[,i], na.rm = T)
}

df<-d
df_t<-t(df)
df_t<-data.frame(df_t)
df_t<-cbind(df_t,info)
#for DAT/CU
frq<-df_t %>% filter(X969 > 0)
frq<-frq[,-ncol(frq)+7]
df_tt<-data.frame(t(frq))
new_info<-df_tt[969:975,]
df_tt<-df_tt[-(969:975),]
###
frq<-df_t %>% filter(X1556 > 0)
frq<-frq[,-ncol(frq)+7]
df_tt<-data.frame(t(frq))
new_info<-df_tt[1556:1562,]
df_tt<-df_tt[-(1556:1562),]

#f=ncol(df_tt)
#for(i in 1:f){
#  df_tt[,i] = as.character(df_tt[,i])
#}
r<-lapply(df_tt, function(df_tt) gsub("2", "1", df_tt))
r<-data.frame(r)

for(i in 1:ncol(r)){
  r[,i] = as.numeric(r[,i])
}

#for(i in 1:ncol(r)){
#  r[,i] = r[,i] -1
#}

for(i in 1:ncol(r)){
  r[(sample_number + 1),i]<-sum(r[,i], na.rm = T)
}
#r<-rbind(r,new_info)
r_t<-t(r)
r_t<-data.frame(r_t)
r_t[,(sample_number + 1)] = as.numeric(r_t[,(sample_number + 1)])
##for DAT, CU
r_t<-r_t %>% filter(X969 < sample_number)
r<-t(r_t)
r<-data.frame(r)
r<-r[-(969:976),]
####
r_t<-r_t %>% filter(X1556 < sample_number) %>% filter(X1556 > 1)
r<-t(r_t)
r<-data.frame(r)
r<-r[-1556,]

#variant 2, 1, 0
#r<-df_tt
#for(i in 1:ncol(r)){
#  r[,i] = as.numeric(r[,i])
#}

list<-master_final %>% select(FID)
w<-cbind(list,r)
cov<-master_final %>% select(FID, sex, age, PC1, PC2, PC3, RCL40,status)
df <- merge(w, cov, by="FID")
#df <- df %>% drop_na(apoe_binary)
mod_summaries <- list()
number= ncol(df)-ncol(cov)
for(i in 2:number) {
  df[,i] = as.numeric(df[,i])
}

for(i in 2:number) {
  test <- df %>% select(c(status, sex, age, PC1, PC2, PC3, RCL40))
  a<-data.frame(df[,i])
  c<-cbind(test,a)
  names(c)[8]<-c("gene")
  mod_summaries[[i - 1]] <- summary(t<-glm(formula=as.numeric(RCL40) ~ gene + age + sex + PC1 + PC2 + PC3, family = "binomial", data = c, na.action = na.omit))$coefficients[2,4]
}
regression<-data.frame((mod_summaries))
gene<-names(df)
up=ncol(df); down=ncol(df)-ncol(cov)+1;
gene<-gene[c(-1,-up:-down)]
names(regression)<-gene
regression<-data.frame(t(regression))
regression$gene<-gene
names(regression)[1]<-c("p")
new_info_t<-as.data.frame(t(new_info))
new_info_t$gene<-dimnames(new_info_t)[[1]]
regression<-merge(regression, new_info_t, by="gene")
FDR<-p.adjust(regression$p, "fdr")
BON<-p.adjust(regression$p, "bonferroni")
BON_cytoband<-regression$p*307
res=data.frame(gene=regression$gene, pvals=regression$p, BONF=round(BON,3), FDR=round(FDR,3), BON_c=round(BON_cytoband,8))
res<-merge(new_info_t, res, by="gene")

#HPSE2
HPSE2<-subset(data,AnnotSV_ID=="10_98736468_98737614_DEL_1")
HPSE2_status<-df %>% select(FID, RCL40, HPSE2.2)
nrow(HPSE2_status %>% filter(RCL40 == 1 & HPSE2.2==1))
#SCHLAP1
SCHLAP1<-subset(data,AnnotSV_ID=="2_180701175_180701900_DEL_1")
SCHLAP1_status<-df %>% select(FID, RCL40, SCHLAP1)
nrow(SCHLAP1_status %>% filter(RCL40 == 0 & SCHLAP1==1))
#ZNF7
ZNF7<-subset(data,AnnotSV_ID=="8_144831620_144832494_DEL_1")
ZNF7_status<-df %>% select(FID, RCL40, ZNF7)
nrow(ZNF7_status %>% filter(RCL40 == 1 & ZNF7==1))
#VAPA
VAPA<-subset(data,AnnotSV_ID=="18_9922853_9924273_DEL_1")
VAPA_status<-df %>% select(FID, RCL40, VAPA)
nrow(VAPA_status %>% filter(RCL40 == 1 & VAPA==1))
#122,167,780,980
#calculate odd ratio
t<-glm(formula=RCL40 ~ VAPA + sex + age + PC1 + PC2 + PC3, family = "binomial", data = df, na.action = na.omit)
summary(t)
odd<-exp(cbind(coef(t), confint.default(t)))
exp(coef(t))
exp(confint(t))

#BON_c 
data<-read.table("exclude_intergenic_variant_1_0.txt", header=T)
for(i in 1:nrow(data)){
  if (data[i,12]>=1){
    data[i,12]=1
  }
}

#Function definition
theme_step1 <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  theme(title = element_text(family = 'Arial', size = 18, color = 'black'), text = element_text(family = 'Arial', size = 16, color = 'black'),
        axis.title = element_text(family = 'Arial', size = 18, color = 'black'), axis.text = element_text(family = 'Arial', size = 16, color = 'black'), 
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = NA), axis.line = element_line(colour = "black", size = rel(1)),
        legend.background = element_rect(color = 'black'), legend.title = element_text(family = 'Arial', size = 16),
        legend.text = element_text(family = 'Arial', size = 14),
        legend.direction = "vertical", 
        legend.box = c("horizontal", "vertical"),
        legend.spacing.x = unit(0.1, 'cm'),
        plot.margin = unit(c(0.25, 1, 1, 0.5), 'cm'),
        axis.title.y = element_text(margin = margin(r = 10, unit = "pt")))
}

#Plot
Chr<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
Del<-c(1649,1558, 1211, 1102, 1133, 1215, 1407, 1035, 856, 1053, 981, 1057, 618, 542, 618, 760, 932, 507, 926, 548, 309, 498)
Dup<-c(344, 262, 230, 207, 207, 292, 273, 186, 139, 176, 189, 201, 109, 119, 120, 145, 198, 88, 247, 125, 57, 87)
#Ins<-c(653, 611, 533, 442, 377, 469, 473, 385, 350, 415, 371, 433, 270, 194, 263, 249, 278, 203, 277, 239, 136, 179)
Del<-as.data.frame(cbind(Chr,Del)) %>% rename(SV=Del); Del$SV_type<-"DEL"
Dup<-as.data.frame(cbind(Chr,Dup)) %>% rename(SV=Dup); Dup$SV_type<-"DUP"
#Ins<-as.data.frame(cbind(Chr,Ins)) %>% rename(SV=Ins); Ins$SV_type<-"INS"
plot_df2<-rbind(Del,Dup)
#plot_df2$Chr<-factor(plot_df2$Chr, levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) 
plot_df2$SV<-as.numeric(plot_df2$SV)
plot_df2$CHR<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
plot_df2$CHR<-factor(plot_df2$CHR, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")) 
#Figure 1
fig1<-ggplot(plot_df2,aes(x=CHR,y=SV,fill=SV_type))+
  geom_bar(stat="identity") + scale_fill_manual(labels=c("DEL","DUP"), values=c("#D83F31","#219C90")) + ggtitle("Number of CNV per chromosome") + xlab("Chromosome") + ylab("Number of CNV") + theme_step1()
ggsave(file='Figure1.png', plot=fig1, width = 8, height = 4.51)

#Figure 2
data$length<-log(abs(data$SV_length), base=10)
fig2<-ggplot(data, aes(length, fill=SV_type, color=SV_type)) + scale_fill_manual(labels=c("DEL","DUP"), values=c("#D83F31","#219C90")) + geom_density(alpha=0.3) + scale_x_continuous(breaks = c(2, 3, 4, 5), labels = c("100bp", "1kb", "10kb", "100kb")) + theme_step1()
ggsave(file='Figure2.png', plot=fig2)

#figure 3
res<-fread("logistic_regression_result")
res<-fread("logistic_regression_result_DAT_CU")
res<-subset(res, P > 0.00000001)
res<-res %>% rename(P=pvals)
res$SV_start<-as.numeric(res$SV_start)
fig3 <- manhattan_plot(x = res, pval.colname = "P", chr.colname = "SV_chrom", pos.colname = "SV_start", ylim = 10, signif = c(2.4e-06, 1.6e-04), chr.col=c("black", "lightgrey"), label.font.size =6) + theme(plot.title = element_text(hjust = 0.5))  + theme_step1() 
fig3 <- fig3 + theme(plot.title = element_text(hjust = 0.5))  + theme_step1() 
ggsave(file='Figure3.png', plot=fig3)
ggsave(file='manhattan_dash.png', plot=fig3)


pdf("manhattan.pdf", width = 10, height = 10)
manhattan_plot(x = res, pval.colname = "P", chr.colname = "SV_chrom", pos.colname = "SV_start", ylim = 10, signif = c(2.4e-06, 1.6e-04), chr.col=c("black", "lightgrey"), label.font.size =6) + theme(plot.title = element_text(hjust = 0.5))  + theme_step1() 


#figure merge
plot1<-grid.arrange(fig1, fig2, fig3, ncol=2, nrow=2)
ggsave(file='Figure1.png', plot=plot1)


res<-("")
res<-res %>% rename(P=pvals)
res$SV_start<-as.numeric(res$SV_start)
p1.1 <- manhattan_plot(x = res, pval.colname = "P", chr.colname = "SV_chrom", pos.colname = "SV_start", ylim = 10, signif = c(5e-08, 1e-06),chr.col=c("black", "lightgrey"), label.font.size =6) 
p1.1 <-  p1.1 + labs(tabs = "c") + theme(plot.title = element_text(hjust = 0.5))  + theme_step1() 

p1.6 <- ggplot(df_DEG, aes(x = cluster_id, y = gene, fill = logFC)) +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = significance_label, fontface = "bold", angle = 90), color = "black", size = 5, hjust = 0.5,) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_gradient2(low = "seagreen", mid = "white", high = "firebrick2", limits = c(-0.7, 0.7)) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 5)) +
  facet_grid(rsID ~ major_celltype, scales = "free", space = "free") +  
  theme_step1() 

#figure4
CNV<-c("DEL","DEL","DUP","DUP")
Number<-c("20515","98805","4001","24602")
Data<-c("This study","Wang et al., 2023","This study","Wang et al., 2023")
a<-cbind(CNV,Number,Data)
a<-data.frame(a)
a$CNV<-as.character(a$CNV)
a$Data<-as.character(a$Data)
a$Number<-as.numeric(a$Number)
a[1,2]=20515/1555; a[2,2]=231385/16905; a[3,2]=4001/1555; a[4,2]=45839/16905;
plot5<-ggplot(a, aes(x=CNV, y=Number, fill=Data)) + geom_bar(position=position_dodge(0.5),stat="identity", size=.2, width = 0.4) + scale_fill_manual(labels=c("This study","Wang 2023"), values=c("#D83F31","#219C90")) + ylab("Mean of CNV") + theme_step1() 
ggsave(
  filename = 'Mean_CNV.pdf',
  plot = plot5,
  width = 7,
  height = 7,
  device = "pdf"
)









fig3<-manhattan(res, chr="SV_chrom", bp="SV_start", snp="AnnotSV_ID", p="pvals") + labs(tag='C')

grid.arrange(fig1, fig2, fig3, ncol=2)
ggsave(file='Figure3.png', plot=plot2)
fig3<- manhattanly(res, snp = "AnnotSV_ID",p="pvals", chr = "SV_chrom", bp = "SV_start", gene="gene")
manhattan(res, chr="SV_chrom", bp="SV_start", snp="AnnotSV_ID", p="Frequnecy", logp = FALSE)
manhattan(res, chr="SV_chrom", bp="SV_start", snp="AnnotSV_ID", p="Bon")
#Figure 3
don <- res %>% group_by(SV_chrom) %>% summarise(chr_len=max(SV_start)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% 
  select(-chr_len) %>%  left_join(res, ., by=c("SV_chrom"="SV_chrom")) %>% arrange(SV_chrom, SV_start) %>% mutate(BPcum=SV_start + tot) %>% rename(P=pvals, chromosome=BPcum)

axisdf = don %>%
  group_by(SV_chrom) %>%
  summarize(center=( max(chromosome) + min(chromosome) ) / 2 )

fig3<-ggplot(don, aes(x=chromosome, y=-log10(P))) + geom_point( aes(color=as.factor(SV_chrom)), alpha=0.8, size=1.3) + scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) + scale_x_continuous( label = axisdf$SV_chrom, breaks= axisdf$center ) + scale_y_continuous(limits = c(0,8)) +  theme_bw() +
  theme(legend.position="none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + labs(tag='C')