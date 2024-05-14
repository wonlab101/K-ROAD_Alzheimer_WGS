import os
import time
import hail as hl
from pprint import pprint

from hail.plot import show
from gnomad_qc.v2.resources.sample_qc import *
from gnomad_qc.v2.resources import get_gnomad_data

import argparse
import sys
import logging
import csv
import pandas as pd
import numpy as np
import glob
import time
from tqdm import tqdm


from gnomad.utils.filtering import filter_to_autosomes
from pprint import pprint

from bokeh.io import output_notebook, show
from bokeh.layouts import gridplot
from bokeh.models import Span
from bokeh.transform import factor_cmap, factor_mark
from bokeh.plotting import figure, output_file, save

os.chdir("/data1/mycho/WGS_AD_2011/")

vcf_dir = "/data2/Oneomics/VCF_hg38.1-6/batch1-7/"
fn_vcf = "SMC_WGS.batch1-7.vcf.gz"

#dir_output = "/data2/Oneomics/VCF_hg38.1-6/batch1-7/data"
dir_output = "/data1/mycho/WGS_AD_2011/4.VCF.QCed.set1-7/"
prefix = "SMC_batch1-7_hg38"
fn_LCR="/data1/mycho/Source/annotation/LCR-hs38.auto_XY.bed"
fn_meta="/data1/mycho/WGS_AD_2011/1.data/phenotype/WGS_1824.pheno_update.230726"

config = { 'spark.driver.memory':'200g', 'spark.local.dir' : '/data1/mycho/WGS_AD_2011/temp'}
hl.init(master='local[60]', spark_conf=config,  default_reference='GRCh38')

###01.1 load matrix file
mt = hl.read_matrix_table(dir_output +'/data/' + prefix + '.VQ1.GQ.mt')

###02. sampleQC
###02.1 sampleQC1
#use hail sampleQC & variantQC function
mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                (mt.variant_qc.AF[1] > 0.001) &
                (mt.variant_qc.call_rate > 0.95))

print('After filter, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))

mt=hl.sample_qc(mt, name='sample_qc')

#vq_ht = mt.variant_qc.call_rate
mean_ht = mt.sample_qc.dp_stats.mean
call_ht = mt.sample_qc.call_rate
n_call = mt.sample_qc.n_called
n_het = mt.sample_qc.n_het


print(dir_output + 'data_qc/' + prefix + '_sample_qc_mean_dp')
mean_ht.export(dir_output + 'data_qc/' + prefix + '_sample_qc_mean_dp')
print(dir_output + 'data_qc/' + prefix  +'_sample_qc_call_rate')
call_ht.export(dir_output + 'data_qc/' + prefix  +'_sample_qc_call_rate')
print(dir_output + 'data_qc/' + prefix + '_sample_qc_n_called')
n_call.export(dir_output + 'data_qc/' + prefix + '_sample_qc_n_called')
print(dir_output + 'data_qc/' + prefix + '_sample_qc_n_het')
n_het.export(dir_output + 'data_qc/' + prefix + '_sample_qc_n_het')

prefix_sample_call_rate = prefix + '_sample_qc_call_rate'
df_sample_qc_call_rate = pd.read_csv(dir_output + 'data_qc/' + prefix_sample_call_rate, sep='\t')
df_sample_qc_call_rate.columns = ['s', 'call_rate']

prefix_sample_dp = prefix  +'_sample_qc_mean_dp'
df_sample_qc_mean_dp = pd.read_csv(dir_output + 'data_qc/' + prefix_sample_dp , sep='\t')
df_sample_qc_mean_dp.columns = ['s', 'mean_dp']

prefix_sample_n_called = prefix  +'_sample_qc_n_called'
df_sample_qc_n_called = pd.read_csv(dir_output + 'data_qc/' + prefix_sample_n_called , sep='\t')
df_sample_qc_n_called.columns = ['s', 'n_called']

prefix_sample_n_het = prefix  +'_sample_qc_n_het'
df_sample_qc_n_het = pd.read_csv(dir_output + 'data_qc/' + prefix_sample_n_het , sep='\t')
df_sample_qc_n_het.columns = ['s', 'n_het']

df_sample_qc = pd.merge(df_sample_qc_mean_dp, df_sample_qc_call_rate, how='inner', on='s')
df_sample_qc = pd.merge(df_sample_qc, df_sample_qc_n_called, how='inner', on='s')
df_sample_qc = pd.merge(df_sample_qc, df_sample_qc_n_het, how='inner', on='s')

df_sample_qc.head()

df_sample_qc['missingness_rate'] = 1- df_sample_qc['call_rate']
df_sample_qc['het_rate'] = df_sample_qc['n_het']/df_sample_qc['n_called']

df_sample_qc.to_csv(dir_output + 'data_qc/' + prefix + '_sample_qc', sep='\t', index=False)

##Call rate & mean dp 
sp_qc = dir_output + 'data_qc/' + prefix + '_sample_qc'
sp_qc_table= (hl.import_table(sp_qc, impute=True).key_by('s'))

sp_qc_table.show(5)

### 1) sample level call rate
p = hl.plot.histogram(sp_qc_table.call_rate, range=(.80,1),
                      legend='Call Rate', title = 'Call Rate')

p.renderers.extend(
    [Span(location=0.90 , dimension='height', line_color='red', line_width=1)
     ])

show(p)

output_file(dir_output + 'plot/' + prefix + '_Sample_qc_call_rate.html')
save(p)


### 2) low coverage (dp_mean)

p = hl.plot.histogram(sp_qc_table.mean_dp, range=(10,70),
                      legend='Mean Sample DP', title = 'Mean Sample DP')

p.renderers.extend(
    [Span(location=4 , dimension='height', line_color='red', line_width=1)
     ])

show(p)


output_file(dir_output + 'plot/' + prefix + '_Sample_qc_mean_dp.html')
save(p)


### 3) scatter plot
'''
p = hl.plot.scatter(sp_qc_table.mean_dp,sp_qc_table.call_rate, 
                    xlabel='Mean DP', ylabel='Call Rate', title='Dp-Call rate scatter plot')

p.renderers.extend(
    [Span(location=0.9 , dimension='width', line_color='red', line_width=1),
     Span(location=15, dimension='height', line_color='red', line_width=1),
    ])


show(p)

output_file(dir_output + 'plot/' + prefix + '_Sample_qc_scatter_plot.html')
save(p)

med_het = sp_qc_table.aggregate(hl.agg.approx_median(sp_qc_table['het_rate']))
sd_het = sp_qc_table.aggregate(hl.agg.stats(sp_qc_table['het_rate'])).stdev

med_miss = sp_qc_table.aggregate(hl.agg.approx_median(sp_qc_table['missingness_rate']))
sd_miss = sp_qc_table.aggregate(hl.agg.stats(sp_qc_table['missingness_rate'])).stdev

p = hl.plot.scatter(sp_qc_table.het_rate,sp_qc_table.missingness_rate,
                    xlabel='Het Rate', ylabel='missingness rate', title='Missing-Het rate scatter plot')

p.renderers.extend(
    [Span(location=med_het+4*sd_het, dimension='height', line_color='red', line_width=1),
     Span(location=med_het-4*sd_het, dimension='height', line_color='red', line_width=1),
     Span(location=med_miss-4*sd_miss, dimension='width', line_color='red', line_width=1),
     Span(location=med_miss-4*sd_miss, dimension='width', line_color='red', line_width=1),
         ])


show(p)

output_file(dir_output + 'plot/' + prefix + '_Sample_qc_missing_het_scatter_plot.html')
save(p)
'''

### 4) filtering

mt_filtered_dp_callrate = sp_qc_table.filter((sp_qc_table.mean_dp >= 15)
                                             & (sp_qc_table.call_rate >= 0.9),
                                             keep=False)

mt_filtered_dp_callrate.export(dir_output + 'filtered_sample/' + prefix + 'filtered_dp_callrate')

## Sex check
mt = hl.read_matrix_table(dir_output +'/data/' + prefix + '.VQ1.GQ.mt')
### 01) run imputed sex

mt_sampleQC1 = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                (mt.variant_qc.AF[1] > 0.001) &
                (mt.variant_qc.call_rate > 0.95))

mt_sampleQC1 = hl.sample_qc(mt_sampleQC1, name='sample_qc')
mt_X = hl.filter_intervals(mt_sampleQC1, [hl.parse_locus_interval('chrX')])


#add split multi hts : based on hail doc 
#(https://hail.is/docs/0.2/methods/genetics.html#hail.methods.impute_sex)
x_nonpar=mt_X.locus.in_x_nonpar()
mt_X=mt_X.annotate_rows(x_nonpar=x_nonpar)
mt_X = mt_X.filter_rows(mt_X.x_nonpar==True)

#female : F < 0.2, male : F > 0.8, aaf_threshold=0.05
sex_ht = hl.impute_sex(mt_X.GT, female_threshold = 0.2, male_threshold = 0.8, aaf_threshold=0.05)

sex_colnames = ['f_stat', 'is_female']
sex_ht = sex_ht.select(*sex_colnames)

mt_X = mt_X.annotate_cols(**sex_ht[mt_X.col_key])

mt_X.f_stat.export(dir_output + 'data_qc/' + prefix + '_f_stat')

X_cov=(dir_output + '/data_qc/' + prefix + '_X_coverage')
mt_X=hl.sample_qc(mt_X, name='sample_qcX')

mt_X.sample_qcX.dp_stats.mean.export(X_cov)

#print('After filter, %d samples and %d variants remain.' % (mt_X.count_cols(), mt_X.count_rows()))

# chr 20 mean coverage()
mt_chr20 = hl.filter_intervals(mt_sampleQC1, [hl.parse_locus_interval('chr20')])
#print('After filter, %d samples and %d variants remain.' % (mt_chr20.count_cols(), mt_chr20.count_rows()))

chr20_cov=(dir_output + '/data_qc/' + prefix + '_chr20_coverage')

mt_chr20=hl.sample_qc(mt_chr20, name='sample_qc20')

mt_chr20.sample_qc20.dp_stats.mean.export(chr20_cov)

#nomalized Y coverage by chr Y
mt_Y = hl.filter_intervals(mt, [hl.parse_locus_interval('chrY')])

print('After filter, %d samples and %d variants remain.' % (mt_Y.count_cols(), mt_Y.count_rows()))

y_nonpar=mt_Y.locus.in_y_nonpar()
mt_Y=mt_Y.annotate_rows(y_nonpar=y_nonpar)
mt_Y = mt_Y.filter_rows(mt_Y.y_nonpar==True)


Y_cov=(dir_output + '/data_qc/' + prefix + '_Y_coverage')

mt_Y = hl.sample_qc(mt_Y, name='sample_qcY')

mt_Y.sample_qcY.dp_stats.mean.export(Y_cov)

chr20_cov=(dir_output + 'data_qc/' + prefix + '_chr20_coverage')
df_chr20 = pd.read_csv(chr20_cov, sep='\t')
df_chr20.columns = ['s', 'chr20_cov']
df_chr20.to_csv(chr20_cov, sep='\t', index=False)
df_chr20.head()

X_cov=(dir_output + 'data_qc/' + prefix + '_X_coverage')
df_chrX = pd.read_csv(X_cov, sep='\t')
df_chrX.columns = ['s', 'X_cov']
df_chrX.to_csv(X_cov, sep='\t', index=False)
df_chrX.head()

Y_cov=(dir_output + 'data_qc/' + prefix + '_Y_coverage')
df_chrY = pd.read_csv(Y_cov, sep='\t')
df_chrY.columns = ['s', 'Y_cov']
df_chrY.to_csv(Y_cov, sep='\t', index=False)
df_chrY.head()

### normalized X coverage
df_chrX_chr20 = pd.merge(df_chrX, df_chr20, how='inner', on='s')
df_chrX_chr20['normalized_X_coverage'] = df_chrX_chr20['X_cov'] /df_chrX_chr20['chr20_cov']
df_normalized_X_coverage = df_chrX_chr20[['s', 'normalized_X_coverage']]
normalized_X_coverage = (dir_output + 'data_qc/' + prefix + '_normalized_X_coverage')
df_normalized_X_coverage.to_csv(normalized_X_coverage, sep='\t', index=False)

### normalized Y coverage
df_chrY_chr20 = pd.merge(df_chrY, df_chr20, how='inner', on='s')
df_chrY_chr20['normalized_Y_coverage'] = df_chrY_chr20['Y_cov'] /df_chrY_chr20['chr20_cov']
df_normalized_Y_coverage = df_chrY_chr20[['s', 'normalized_Y_coverage']]
normalized_Y_coverage = (dir_output + 'data_qc/' + prefix + '_normalized_Y_coverage')
df_normalized_Y_coverage.to_csv(normalized_Y_coverage, sep='\t', index=False)

normalized_X_coverage = (dir_output + 'data_qc/' + prefix + '_normalized_X_coverage')
normalized_Y_coverage = (dir_output + 'data_qc/' + prefix + '_normalized_Y_coverage')
f_stat = (dir_output + 'data_qc/' + prefix + '_f_stat')
impute_sex = (dir_output + 'data_qc/' + prefix + '_impute_sex')

df_normal_X = pd.read_csv(normalized_X_coverage, sep='\t')
df_normal_Y = pd.read_csv(normalized_Y_coverage, sep='\t')
df_f_stat = pd.read_csv(f_stat, sep='\t')

df_impute_sex = pd.merge(df_f_stat, df_normal_X, how='inner', on='s')
df_impute_sex = pd.merge(df_impute_sex, df_normal_Y, how='inner', on='s')

df_impute_sex.to_csv(impute_sex,sep='\t', index=False)
df_impute_sex.head()

print(df_impute_sex)

conditionlist = [
    (df_impute_sex['f_stat'] <= 0.5) &  (df_impute_sex['normalized_Y_coverage'] <= 0.1),
    (df_impute_sex['f_stat'] >= 0.6) &  (df_impute_sex['normalized_Y_coverage'] >= 0.1),
    ((df_impute_sex['f_stat'] < 0.4) & (df_impute_sex['normalized_Y_coverage'] > 0.1)),
    ((df_impute_sex['f_stat'] >= 0.5) & (df_impute_sex['normalized_Y_coverage'] <= 0.1)) |
    ((df_impute_sex['f_stat'] >= 0.4) & (df_impute_sex['f_stat'] <= 0.6)
     & (df_impute_sex['normalized_Y_coverage'] > 0.1))]
choicelist = ['Female', 'Male', 'sex_aneuploidy', 'ambiguous_sex']
df_impute_sex['impute_sex'] = np.select(conditionlist, choicelist, default='NA')


impute_sex = (dir_output + 'data_qc/' + prefix + '_impute_sex')
df_impute_sex=pd.read_csv(impute_sex, sep='\t')

conditionlist = [
    (df_impute_sex['f_stat'] <= 0.6),
    (df_impute_sex['f_stat'] >= 0.8)]
choicelist = ['F', 'M' ]
df_impute_sex['impute_sex'] = np.select(conditionlist, choicelist, default='NA')


df_impute_sex.to_csv(impute_sex,sep='\t', index=False)
df_impute_sex.head()

impute_sex_ht = (hl.import_table(impute_sex, impute=True).key_by('s'))

#, 'normalized_X_coverage': hl.tfloat64, 'normalized_Y_coverage':hl.tfloat64

impute_sex_ht.show(5)

### 3) ambiguous sex

impute_sex = (dir_output + 'data_qc/' + prefix + '_impute_sex')
df_impute_sex=pd.read_csv(impute_sex, sep='\t')
df_impute_sex.head()

df_sex = pd.read_csv(fn_meta, sep='\t')
df_sex = df_sex[['FID', 'sex']]
df_sex.columns = ['s', 'sex']
df_sex = pd.merge(df_impute_sex, df_sex, how='inner', on = 's')
df_sex.head()

df_ambiguous_sex_lower = df_sex['f_stat'] > 0.3
df_ambiguous_sex_upper = df_sex['f_stat'] < 0.8

df_ambiguous_sex = df_sex[ df_ambiguous_sex_lower & df_ambiguous_sex_upper ]

df_ambiguous_sex

df_ambiguous_sex.to_csv(dir_output + '/filtered_sample/' + prefix + '_ambiguous_sex', sep='\t', index=False)

### 4) unmatched sex
impute_sex = (dir_output + 'data_qc/' + prefix + '_impute_sex')
df_impute_sex=pd.read_csv(impute_sex, sep='\t')

df_sex = pd.read_csv(fn_meta, sep='\t')
df_sex = df_sex[['FID', 'sex']]
df_sex['SEX'] = np.where(df_sex['sex']==2, 'F','M')
df_sex = df_sex[['FID','SEX']]
df_sex.columns = ['s', 'SEX']
df_sex.head()

df_mis = pd.merge(df_impute_sex, df_sex, how='inner', on = 's')
df_mis = df_mis[df_mis['impute_sex'] != df_mis['SEX']]
print(df_mis)

df_mis.to_csv(dir_output + '/filtered_sample/' + prefix + '_unmatched_sex', sep='\t', index=False)

###save file
####save 
mt_sampleQC1.write(dir_output +'/data/' + prefix +'.sampleQC1.mt', overwrite=True)



###02.2 sampleQC2
## 03.02.01 KING relatedness / PCA

# read VQSR filtered matrix
mt_sampleQC2 = hl.read_matrix_table(dir_output +'/data/' + prefix + '.VQ1.GQ.mt')
print('%d samples and %d varianxts are uploaded.' % (mt_sampleQC2.count_cols(), mt_sampleQC2.count_rows()))

mt_sampleQC2=hl.variant_qc(mt_sampleQC2, name='variant_qc')

mt_sampleQC2 = mt_sampleQC2.filter_rows((hl.len(mt_sampleQC2.alleles) == 2) &
                                        hl.is_snp(mt_sampleQC2.alleles[0], mt_sampleQC2.alleles[1]) &
                                        (mt_sampleQC2.variant_qc.AF[1] > 0.001) &
                                        (mt_sampleQC2.variant_qc.call_rate > 0.95))

INBREEDING_COEFF_HARD_CUTOFF = -0.25

mt_sampleQC2 = filter_to_autosomes(mt_sampleQC2)

mt_sampleQC2_AT5 = mt_sampleQC2.filter_rows(mt_sampleQC2.variant_qc.AF[1] > 0.05)

mt_sampleQC2.write(dir_output +'/data/' + prefix +'.sampleQC2.VQ1.mt', overwrite=True)

mt_sampleQC2 = hl.read_matrix_table(dir_output +'/data/' + prefix +'.sampleQC2.VQ1.mt')

print('After filter, %d samples and %d variants remain.'
      % (mt_sampleQC2.count_cols(), mt_sampleQC2.count_rows()))

mt_sampleQC2_AT5.write(dir_output +'/data/' + prefix +'.sampleQC2_AT5.VQ1.mt', overwrite=True)

mt_sampleQC2_AT5 = hl.read_matrix_table(dir_output +'/data/' + prefix +'.sampleQC2_AT5.VQ1.mt')

print('After filter, %d samples and %d variants remain.'
      % (mt_sampleQC2_AT5.count_cols(), mt_sampleQC2_AT5.count_rows()))

mt_sampleQC2_AT5 = hl.read_matrix_table(dir_output +'/data/' + prefix +'.sampleQC2_AT5.VQ1.mt')
hl.export_vcf(mt_sampleQC2_AT5, dir_output +'/data/' + prefix +'.sampleQC2_AT5.VQ1_230705.vcf.bgz')


## 03.02.02 Genotype Concordance


#Read SampleQC2_AT5
mt_sampleQC2_AT5 = hl.read_matrix_table(dir_output +'/data/' + prefix +'.sampleQC2_AT5.VQ1.mt')

#Read GWAS data
fn_gwas="/data1/mycho/WGS_AD_2011/4.VCF.QCed.set1-7/data/SMC_merged.maf01.rmdup.completed.anno.TOPMed_hg38_concordance.vcf"
hl.import_vcf(fn_gwas, force_bgz=True, reference_genome='GRCh38', contig_recoding={'1': 'chr1', '2':'chr2','3':'chr3','4':'chr4','5':'chr5','6':'chr6','7':'chr7','8':'chr8','9':'chr9','10':'chr10','11': 'chr11', '12':'chr12','13':'chr13','14':'chr14','15':'chr15','16':'chr16','17':'chr17','18':'chr18','19':'chr19','20':'chr20','21':'chr21','22':'chr22'}).write(dir_output +'/data/' + prefix +'.gwas_chip_data.mt', overwrite=True)

mt_GWAS = hl.read_matrix_table(dir_output +'/data/' + prefix +'.gwas_chip_data.mt')

#### 2) QC
mt_GWAS = mt_GWAS.key_rows_by(**hl.min_rep(mt_GWAS.locus, mt_GWAS.alleles))
print('%d samples and %d variants are uploaded.' % (mt_GWAS.count_cols(), mt_GWAS.count_rows()))
mt_GWAS = hl.split_multi_hts(mt_GWAS)
print('After split multi hts, %d samples and %d variants remain.' % (mt_GWAS.count_cols(), mt_GWAS.count_rows()))

#hail variantQC & sampleQC
mt_GWAS=hl.sample_qc(mt_GWAS)
mt_GWAS=hl.variant_qc(mt_GWAS)

mt_GWAS = mt_GWAS.filter_rows((hl.len(mt_GWAS.alleles) == 2) & 
                hl.is_snp(mt_GWAS.alleles[0], mt_GWAS.alleles[1]) &
                (mt_GWAS.variant_qc.call_rate > 0.95))
print('After filter, %d samples and %d variants remain.' % (mt_GWAS.count_cols(), mt_GWAS.count_rows()))

mt_GWAS = filter_to_autosomes(mt_GWAS)
print('After filter, %d samples and %d variants remain.' % (mt_GWAS.count_cols(), mt_GWAS.count_rows()))

mt_GWAS.write(dir_output +'/data/' + prefix +'.gwas_chip_data.QCed.mt', overwrite=True)

#### 3) calculate concordance
summary, samples, variants = hl.concordance(mt_sampleQC2_AT5, mt_GWAS)
samples.write( dir_output +'/data/' + prefix + '.concordacne.samples', overwrite=True)
variants.write( dir_output +'/data/' + prefix + '.concordacne.variants', overwrite=True)

with open(dir_output +'/data/' + prefix + '.concordacne.summary', 'w', newline='') as f: 
    writer = csv.writer(f) 
    writer.writerow(summary) 

samples = hl.read_table(dir_output +'/data/' + prefix + '.concordacne.samples')
variants = hl.read_table(dir_output +'/data/' + prefix + '.concordacne.variants')

# non-ref disconcordance
sample_non_ref_con = (samples.concordance[3][3]+samples.concordance[4][4]) / (samples.concordance[1][3]+samples.concordance[1][4]+
                                                                             samples.concordance[2][3]+samples.concordance[2][4]+
                                                                             samples.concordance[3][3]+samples.concordance[3][4]+
                                                                             samples.concordance[4][3]+samples.concordance[4][4])

samples = samples.annotate(non_ref_concordance = sample_non_ref_con)

samples_filtered_non_ref_concordance = samples.filter(samples.non_ref_concordance < 0.9)

print('After filter, %d samples filtered.' % (samples_filtered_non_ref_concordance.count()))
print(samples_filtered_non_ref_concordance.s.show())

samples_filtered_non_ref_concordance.export(dir_output + '/filtered_sample/' + prefix + '_non_ref_concordance', delimiter='\t')




