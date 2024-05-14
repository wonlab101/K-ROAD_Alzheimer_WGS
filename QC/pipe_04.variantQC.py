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
fn_sample_list="/data1/mycho/WGS_AD_2011/4.VCF.QCed.set1-7/data/sample.list"

config = { 'spark.driver.memory':'200g', 'spark.local.dir' : '/data1/mycho/WGS_AD_2011/temp'}
hl.init(master='local[60]', spark_conf=config,  default_reference='GRCh38')
###01.1 load matrix file
mt_variantQC1 = hl.read_matrix_table(dir_output +'/data/' + prefix + '.VQ1.GQ.mt')

print('%d samples and %d varianxts are uploaded.' % (mt_variantQC1.count_cols(), mt_variantQC1.count_rows()))

samples_to_remove = {"WGS_1810","WGS_0179","WGS_0298","WGS_0468","WGS_0808","WGS_0959","WGS_1013","WGS_1064","WGS_1108","WGS_1109","WGS_1366","WGS_1370","WGS_1637","WGS_1663","WGS_1700","WGS_1722","WGS_1810","WGS_0345","WGS_0972","WGS_0463","WGS_0477","WGS_0496","WGS_0515","WGS_0519","WGS_0542","WGS_0560","WGS_0806","WGS_1133","WGS_0546","WGS_0897","WGS_0130","WGS_0221","WGS_0666","WGS_0919","WGS_1209","WGS_1375","WGS_0984","WGS_0955","WGS_0973","WGS_0997","WGS_1042","WGS_0941","WGS_1006","WGS_1014","WGS_1314","WGS_1010","WGS_1028","WGS_1015","WGS_0407","WGS_0455","WGS_1049","WGS_1385","WGS_1297","WGS_0004","WGS_0028","WGS_0298","WGS_1064","WGS_1637","WGS_1663","WGS_1722","WGS_0569"}

set_to_remove = hl.literal(samples_to_remove)
mt_variantQC1 = mt_variantQC1.filter_cols(~set_to_remove.contains(mt_variantQC1['s']))

print('After Sample QC, %d samples and %d variants remain.' % (mt_variantQC1.count_cols(), mt_variantQC1.count_rows()))

#hail variantQC & sampleQC
mt_variantQC1=hl.variant_qc(mt_variantQC1)
mt_variantQC1=hl.sample_qc(mt_variantQC1)

HARD_CUTOFF = 0

mt_variantQC1 = mt_variantQC1.filter_rows(mt_variantQC1.variant_qc.AC[1] > HARD_CUTOFF)

mt_variantQC1.write(dir_output +'/data/' + prefix +'.variantQC1.mt', overwrite=True)

mt_variantQC1 = hl.read_matrix_table(dir_output +'/data/' + prefix +'.variantQC1.mt')

print('After AC filter, %d samples and %d variants remain.'
      % (mt_variantQC1.count_cols(), mt_variantQC1.count_rows()))


# 02 VaraintQC 1
## 02.01. inbreeding coefficient


p = hl.plot.histogram(mt_variantQC1.info.InbreedingCoeff,
                      legend='InbreedingCoeff', title='Inbreeding coefficient')

p.renderers.extend(
                          [Span(location=-0.3 , dimension='height', line_color='red', line_width=1),
                               ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC1_InbreedingCoeff.html')
save(p)

INBREEDING_COEFF_HARD_CUTOFF = -0.3
mt_variantQC1 = mt_variantQC1.filter_rows(mt_variantQC1.info.InbreedingCoeff >= INBREEDING_COEFF_HARD_CUTOFF)


print('After inbreeding coeff filter, %d samples and %d variants remain.' % (mt_variantQC1.count_cols(), mt_variantQC1.count_rows()))

## 02.03. Missingness rate

p = hl.plot.histogram(mt_variantQC1.variant_qc.call_rate,
                      legend='call_rate', title='Missingness rate')

p.renderers.extend(
                          [Span(location=0.9 , dimension='height', line_color='red', line_width=1),
                               ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC1_Missingness.html')
save(p)

Call_rate_HARD_CUTOFF = 0.9

mt_variantQC1 = mt_variantQC1.filter_rows(mt_variantQC1.variant_qc.call_rate >= Call_rate_HARD_CUTOFF)

print('After callrate filter, %d samples and %d variants remain.' % (mt_variantQC1.count_cols(), mt_variantQC1.count_rows()))

## 02.04. HWE (only control sample)
### only control samples

mt_variantQC1_NC = mt_variantQC1.filter_cols(mt_variantQC1.pheno.DX == "CU" ) 
mt_variantQC1_NC.write(dir_output +'/data/' + prefix +'.variantqc.NC.mt', overwrite=True)
mt_variantQC1_NC = hl.read_matrix_table(dir_output +'/data/' + prefix +'.variantqc.NC.mt')
print('After CU filter, %d samples filtered %d variants filtered.'
      % (mt_variantQC1_NC.count_cols(), mt_variantQC1_NC.count_rows()))

#p_value_hwe (float64) – p-value from test of Hardy-Weinberg equilibrium. 
#See functions.hardy_weinberg_test() for details.
p = hl.plot.histogram(mt_variantQC1_NC.variant_qc.p_value_hwe, 
                      legend='p_value_hwe', title='P value, HWE')

p.renderers.extend(
    [Span(location=1e-09 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC1_P_value_HWE.html')
save(p)

#het_freq_hwe (float64) – Expected frequency of heterozygous samples under Hardy-Weinberg equilibrium. 
#See functions.hardy_weinberg_test() for details.
p = hl.plot.histogram(mt_variantQC1_NC.variant_qc.het_freq_hwe, 
                      legend='het_freq_hwe', title='Het freq HWE')

output_file(dir_output + 'plot/' + prefix + '_VariantQC1_Het_freq_HWE.html')
save(p)

# cutoff from MIGEN (Prof.Won)
HWE_HARD_CUTOFF = 1e-09

mt_variantQC1_HWE = mt_variantQC1_NC.filter_rows(mt_variantQC1_NC.variant_qc.p_value_hwe 
                                                      <= HWE_HARD_CUTOFF)

print('After HWE filter,  %d variants filtered.' % (mt_variantQC1_HWE.count_rows()))


mt_variantQC1 = mt_variantQC1.filter_rows(hl.is_defined(mt_variantQC1_HWE.rows()[mt_variantQC1.row_key]), keep=False) 

print('After HWE merge all data filter, %d samples and %d variants remain.' 
      % (mt_variantQC1.count_cols(), mt_variantQC1.count_rows()))


mt_variantQC1.write(dir_output +'/data/' + prefix +'.variantQC2.mt', overwrite=True)


# 03. Variant QC 2
## 03.01. split SNP / indel
mt_variantQC2 = hl.read_matrix_table(dir_output +'/data/' + prefix +'.variantQC2.mt')
mt_variantQC2 = hl.variant_qc(mt_variantQC2)
mt_variantQC2 = hl.sample_qc(mt_variantQC2)

print(dir_output +'/data/' + prefix +'.variantQC2.mt')
print('After filter, %d samples and %d variaznts remain.' 
      % (mt_variantQC2.count_cols(), mt_variantQC2.count_rows()))

mt_variantQC2_SNP = mt_variantQC2.filter_rows(hl.is_snp(mt_variantQC2.alleles[0], mt_variantQC2.alleles[1]))
print('After SNP filter, %d samples and %d variants remain.' 
      % (mt_variantQC2_SNP.count_cols(), mt_variantQC2_SNP.count_rows()))


mt_variantQC2_INDEL = mt_variantQC2.filter_rows(hl.is_indel(mt_variantQC2.alleles[0], mt_variantQC2.alleles[1]))
print('After INDEL filter, %d samples and %d variants remain.' 
      % (mt_variantQC2_INDEL.count_cols(), mt_variantQC2_INDEL.count_rows()))

## 03.02. SNP filtering
#### 1) QD
p = hl.plot.histogram(mt_variantQC2_SNP.info.QD, 
                      legend='QD', title='SNP QD')

p.renderers.extend(
    [Span(location=2 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_SNP_QD.html')
save(p)

HARD_CUTOFF = 2

mt_variantQC2_SNP = mt_variantQC2_SNP.filter_rows(mt_variantQC2_SNP.info.QD >= HARD_CUTOFF)

print('After QD filter, %d samples and %d variants remain.' 
      % (mt_variantQC2_SNP.count_cols(), mt_variantQC2_SNP.count_rows()))

p = hl.plot.histogram(mt_variantQC2_SNP.info.QD, 
                      legend='QD', title='SNP QD, After QC')

p.renderers.extend(
    [Span(location=2 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_SNP_QD_after.html')
save(p)

### 2) SOR
p = hl.plot.histogram(mt_variantQC2_SNP.info.SOR, 
                      legend='SOR', title = 'SNP SOR')

p.renderers.extend(
    [Span(location=3 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_SNP_SOR.html')
save(p)

HARD_CUTOFF = 3

mt_variantQC2_SNP = mt_variantQC2_SNP.filter_rows(mt_variantQC2_SNP.info.SOR <= HARD_CUTOFF)

print('After SOR filter, %d samples and %d variants remain.' % (mt_variantQC2_SNP.count_cols(), mt_variantQC2_SNP.count_rows()))

p = hl.plot.histogram(mt_variantQC2_SNP.info.SOR, 
                      legend='SOR', title ='SNP SOR, After QC')

p.renderers.extend(
    [Span(location=3 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_SNP_SOR_after.html')
save(p)

### 3) FS
p = hl.plot.histogram(mt_variantQC2_SNP.info.FS, 
                      legend='FS', title='SNP FS')

p.renderers.extend(
    [Span(location=60 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_SNP_FS.html')
save(p)

HARD_CUTOFF = 60

mt_variantQC2_SNP = mt_variantQC2_SNP.filter_rows(mt_variantQC2_SNP.info.FS <= HARD_CUTOFF)

print('After FS filter, %d samples and %d variants remain.' 
      % (mt_variantQC2_SNP.count_cols(), mt_variantQC2_SNP.count_rows()))

p = hl.plot.histogram(mt_variantQC2_SNP.info.FS, 
                      legend='FS', title='SNP FS, After QC')

p.renderers.extend(
    [Span(location=60 , dimension='height', line_color='red', line_width=1),
     ])
output_file(dir_output + 'plot/' + prefix + '_VariantQC2_SNP_FS_after.html')
save(p)

### 4) MQ
p = hl.plot.histogram(mt_variantQC2_SNP.info.MQ, 
                      legend='MQ', title = 'SNP MQ')


p.renderers.extend(
    [Span(location=50 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_SNP_MQ.html')
save(p)

HARD_CUTOFF = 50


mt_variantQC2_SNP = mt_variantQC2_SNP.filter_rows(mt_variantQC2_SNP.info.MQ >= HARD_CUTOFF)

print('After MQ filter, %d samples and %d variants remain.' 
      % (mt_variantQC2_SNP.count_cols(), mt_variantQC2_SNP.count_rows()))

p = hl.plot.histogram(mt_variantQC2_SNP.info.MQ, 
                      legend='MQ', title = 'SNP MQ, After QC')

p.renderers.extend(
    [Span(location=50 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_SNP_MQ_after.html')
save(p)

### 5) MQRankSum
p = hl.plot.histogram(mt_variantQC2_SNP.info.MQRankSum, 
                      legend='MQRankSum', title='SNP MQRankSum')

p.renderers.extend(
    [Span(location=-12.5 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_SNP_MQRankSum.html')
save(p)

HARD_CUTOFF = -12.5

mt_variantQC2_SNP = mt_variantQC2_SNP.filter_rows(mt_variantQC2_SNP.info.MQRankSum >= -12.5)

print('After MQRankSum filter, %d samples and %d variants remain.' % (mt_variantQC2_SNP.count_cols(), mt_variantQC2_SNP.count_rows()))

p = hl.plot.histogram(mt_variantQC2_SNP.info.MQRankSum, 
                      legend='MQRankSum', title = 'SNP MQRankSum, After QC')

p.renderers.extend(
    [Span(location=-12.5 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_SNP_MQRankSum_after.html')
save(p)

### 6) ReadPosRankSum

p = hl.plot.histogram(mt_variantQC2_SNP.info.ReadPosRankSum, 
                      legend='ReadPosRankSum', title='SNP ReadPosRankSum')

p.renderers.extend(
    [Span(location=-8.0 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_SNP_ReadPosRankSum.html')
save(p)



HARD_CUTOFF = -8.0

mt_variantQC2_SNP = mt_variantQC2_SNP.filter_rows(mt_variantQC2_SNP.info.ReadPosRankSum >= HARD_CUTOFF)

print('After ReadPosRankSum filter, %d samples and %d variants remain.' % (mt_variantQC2_SNP.count_cols(), mt_variantQC2_SNP.count_rows()))

p = hl.plot.histogram(mt_variantQC2_SNP.info.ReadPosRankSum, 
                      legend='ReadPosRankSum', title='SNP ReadPosRankSum, After QC')

p.renderers.extend(
    [Span(location=-8.0 , dimension='height', line_color='red', line_width=1),
     ])


output_file(dir_output + 'plot/' + prefix + '_VariantQC2_SNP_ReadPosRankSum_after.html')
save(p)


## 03.03. INDEL filtering
#### 1) QD

p = hl.plot.histogram(mt_variantQC2_INDEL.info.QD, 
                      legend='QD', title = 'INDEL QD')

p.renderers.extend(
    [Span(location=2 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_INDEL_QD.html')
save(p)

HARD_CUTOFF = 2

mt_variantQC2_INDEL = mt_variantQC2_INDEL.filter_rows(mt_variantQC2_INDEL.info.QD >= HARD_CUTOFF)

print('After INDEL QD filter, %d samples and %d variants remain.' 
      % (mt_variantQC2_INDEL.count_cols(), mt_variantQC2_INDEL.count_rows()))

p = hl.plot.histogram(mt_variantQC2_INDEL.info.QD, 
                      legend='QD', title='INDEL QD, After QC')

p.renderers.extend(
    [Span(location=2 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_INDEL_QD_after.html')
save(p)

### 2) FS
p = hl.plot.histogram(mt_variantQC2_INDEL.info.FS, 
                      legend='FS', title = 'INDEL FS')

p.renderers.extend(
    [Span(location=200 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_INDEL_FS.html')
save(p)

HARD_CUTOFF = 200

mt_variantQC2_INDEL = mt_variantQC2_INDEL.filter_rows(mt_variantQC2_INDEL.info.FS <= HARD_CUTOFF)

print('After INDEL FS filter, %d samples and %d variants remain.' 
      % (mt_variantQC2_INDEL.count_cols(), mt_variantQC2_INDEL.count_rows()))

p = hl.plot.histogram(mt_variantQC2_INDEL.info.FS, 
                      legend='FS', title='INDEL FS, After QC')

p.renderers.extend(
    [Span(location=200 , dimension='height', line_color='red', line_width=1),
     ])
output_file(dir_output + 'plot/' + prefix + '_VariantQC2_INDEL_FS_after.html')
save(p)

### 3) ReadPosRankSum
p = hl.plot.histogram(mt_variantQC2_INDEL.info.ReadPosRankSum, 
                      legend='ReadPosRankSum', title='INDEL ReadPosRankSum')

p.renderers.extend(
    [Span(location=-20 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_INDEL_QD_after.html')
save(p)

HARD_CUTOFF = -20

mt_variantQC2_INDEL = mt_variantQC2_INDEL.filter_rows(mt_variantQC2_INDEL.info.ReadPosRankSum >= HARD_CUTOFF)

print('After ReadPosRankSum filter, %d samples and %d variants remain.' % (mt_variantQC2_INDEL.count_cols(), mt_variantQC2_INDEL.count_rows()))

p = hl.plot.histogram(mt_variantQC2_INDEL.info.ReadPosRankSum, 
                      legend='ReadPosRankSum', title = 'INDEL ReadPosRankSum, After QC')

p.renderers.extend(
    [Span(location=-20 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_INDEL_ReadPosRankSum_after.html')
save(p)

## 4) MQ
p = hl.plot.histogram(mt_variantQC2_INDEL.info.MQ, 
                      legend='MQ', title = 'INDEL MQ')

p.renderers.extend(
    [Span(location=50 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_INDEL_MQ.html')
save(p)

HARD_CUTOFF = 50

mt_variantQC2_INDEL = mt_variantQC2_INDEL.filter_rows(mt_variantQC2_INDEL.info.MQ >= HARD_CUTOFF)

print('After MQ filter, %d samples and %d variants remain.' % (mt_variantQC2_INDEL.count_cols(), mt_variantQC2_INDEL.count_rows()))

p = hl.plot.histogram(mt_variantQC2_INDEL.info.MQ, 
                      legend='MQ', title = 'INDEL MQ_after QC')


p.renderers.extend(
    [Span(location=50 , dimension='height', line_color='red', line_width=1),
     ])


output_file(dir_output + 'plot/' + prefix + '_VariantQC2_INDEL_MQ_after.html')
save(p)

## 5) MQRankSum
p = hl.plot.histogram(mt_variantQC2_INDEL.info.MQRankSum, 
                      legend='MQRankSum', title = 'INDEL MQRankSum')


p.renderers.extend(
    [Span(location=-12.5 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_INDEL_MQRanksum.html')
save(p)

HARD_CUTOFF = -12.5

mt_variantQC2_INDEL = mt_variantQC2_INDEL.filter_rows(mt_variantQC2_INDEL.info.MQRankSum >= HARD_CUTOFF)

print('After MQRankSum filter, %d samples and %d variants remain.' % (mt_variantQC2_INDEL.count_cols(), mt_variantQC2_INDEL.count_rows()))

p = hl.plot.histogram(mt_variantQC2_INDEL.info.MQRankSum, 
                      legend='MQRankSum', title = 'INDEL MQRankSum_After QC')



p.renderers.extend(
    [Span(location=-12.5 , dimension='height', line_color='red', line_width=1),
     ])

output_file(dir_output + 'plot/' + prefix + '_VariantQC2_INDEL_MQRanksum_after.html')
save(p)

# Export data
print('After filter, %d samples and %d variants remain.' 
      % (mt_variantQC2_SNP.count_cols(), mt_variantQC2_SNP.count_rows()))

print('After filter, %d samples and %d variants remain.' 
      % (mt_variantQC2_INDEL.count_cols(), mt_variantQC2_INDEL.count_rows()))


## merge SNP & INDEL matrix
mt_variantQC2_merged = mt_variantQC2_SNP.union_rows(mt_variantQC2_INDEL)
print('After filter, %d samples and %d variants remain.' 
      % (mt_variantQC2_merged.count_cols(), mt_variantQC2_merged.count_rows()))


#hail variantQC & sampleQC
mt_variantQC2_merged=hl.variant_qc(mt_variantQC2_merged)
mt_variantQC2_merged=hl.sample_qc(mt_variantQC2_merged)

mt_variantQC2_merged.write(dir_output +'/data/' + prefix +'.variantQC2.final.mt', overwrite=True)

mt = hl.read_matrix_table(dir_output +'/data/' + prefix +'.variantQC2.final.mt')
print(dir_output +'/data/' + prefix +'.variantQC2.final.mt')


# SampleQC3
mt=hl.variant_qc(mt, name='variant_qc')

mt = mt.filter_rows((hl.len(mt.alleles) == 2) & 
                                        hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                                        (mt.variant_qc.AF[1] > 0.001) &
                                        (mt.variant_qc.call_rate > 0.99))

mt = filter_to_autosomes(mt)


mt.write(dir_output +'/data/' + prefix +'.variantQC2.final.VQ1.mt', overwrite=True)

mt_sampleQC2 = hl.read_matrix_table(dir_output +'/data/' + prefix +'.variantQC2.final.VQ1.mt')

print('After filter, %d samples and %d variants remain.' 
      % (mt_sampleQC2.count_cols(), mt_sampleQC2.count_rows()))

qc_list = ['n_snp', 'n_singleton', 
           'r_insertion_deletion', 'n_insertion', 'n_deletion',
          'r_het_hom_var', 'n_het', 'n_hom_var',
          'r_ti_tv', 'n_transition', 'n_transversion']

df_result = pd.read_csv(fn_sample_list, sep='\t', header=None)
df_result.columns = ['s']

for qc_name in qc_list:
    tmp_ht = mt_sampleQC2.sample_qc[qc_name]
    tmp_ht.export(dir_output + 'data_qc.sampleQC3/' + prefix + '_sample_qc_' + qc_name)
    df_tmp = pd.read_csv(dir_output + 'data_qc.sampleQC3/' + prefix + '_sample_qc_'+ qc_name, sep='\t')
    df_tmp.columns = ['s', qc_name]
    df_result = pd.merge(df_result, df_tmp, how='inner', on='s')

df_result.to_csv(dir_output + 'data_qc.sampleQC3/' + prefix + '_sample_qc2', sep='\t', index =False)

qc_ht = (hl.import_table(dir_output + 'data_qc.sampleQC3/' + prefix + '_sample_qc2', impute=True).key_by('s'))

qc_ht.show(5)

qc_list2 = ['n_snp', 
           'r_insertion_deletion', 'n_insertion', 'n_deletion',
          'r_het_hom_var', 'n_het', 'n_hom_var',
          'r_ti_tv', 'n_transition', 'n_transversion']

for qc_name in qc_list2:
    median = qc_ht.aggregate(hl.agg.approx_median(qc_ht[qc_name])) 
    sd = qc_ht.aggregate(hl.agg.stats(qc_ht[qc_name])).stdev
    
    mt_sampleQC2_lower = qc_ht.filter(qc_ht[qc_name] <= (median-5*sd))
    mt_sampleQC2_upper = qc_ht.filter(qc_ht[qc_name] >= (median+5*sd))
    
    mt_sampleQC2_total = mt_sampleQC2_lower.union(mt_sampleQC2_upper)
    
 
      
    print('After filter, %d samples filtered.' % (mt_sampleQC2_total.count()))
    
    mt_sampleQC2_total.s.export(dir_output + '/filtered_sample.sampleQC3/' + prefix + '_sampleQC2_' + qc_name)
    
    p = hl.plot.histogram(qc_ht[qc_name], legend=qc_name,  title = qc_name)

    p.renderers.extend(
        [Span(location=median+5*sd, dimension='height', line_color='red', line_width=1),
         Span(location=median-5*sd, dimension='height', line_color='red', line_width=1),
          ])

    output_file(dir_output + 'filtered_sample.sampleQC3/' + prefix + '_SampleQC2_'+ qc_name +'.html')
    save(p)


qc_name='n_singleton'

median = qc_ht.aggregate(hl.agg.approx_median(qc_ht[qc_name])) 
sd = qc_ht.aggregate(hl.agg.stats(qc_ht[qc_name])).stdev
    
mt_sampleQC2_lower = qc_ht.filter(qc_ht[qc_name] <= (median-8*sd))
mt_sampleQC2_upper = qc_ht.filter(qc_ht[qc_name] >= (median+8*sd))
    
mt_sampleQC2_total = mt_sampleQC2_lower.union(mt_sampleQC2_upper)

      
print('After filter, %d samples filtered.' % (mt_sampleQC2_total.count()))
    
mt_sampleQC2_total.s.export(dir_output + '/filtered_sample.sampleQC3/' + prefix + '_sampleQC2_' + qc_name)
    
print(mt_sampleQC2_total.s.show())
    
p = hl.plot.histogram(qc_ht[qc_name], legend=qc_name,  title = qc_name)

p.renderers.extend(
        [Span(location=median+8*sd, dimension='height', line_color='black', line_width=1),
         Span(location=median-8*sd, dimension='height', line_color='black', line_width=1),
         ])


output_file(dir_output + 'filtered_sample.sampleQC3/' + prefix + '_SampleQC2_'+ qc_name +'.html')
save(p)

# n_transition, n_transversion
#qc_name='ntitv'
#p = hl.plot.scatter(qc_ht.n_transition, qc_ht.n_transversion, 
#                    xlabel='n_transition', ylabel='n_transversion')

#output_file(dir_output + 'filtered_sample.sampleQC3/' + prefix + '_SampleQC2_'+ qc_name +'.html')
#save(p)

df_result=pd.DataFrame({'s':[]})


for qc_name in qc_list:
    df_tmp=pd.read_csv(dir_output + '/filtered_sample.sampleQC3/' + prefix + '_sampleQC2_' + qc_name)
    df_tmp[qc_name]=1
    print(qc_name)
    print(df_tmp.empty)
    if df_tmp.empty==False:
        df_result = pd.merge(df_result, df_tmp, how='outer', on='s')


print(len(df_result))
df_result

df_result.s.to_csv(dir_output + '/filtered_sample.sampleQC3/' + prefix + '_sampleQC2', sep='\t', index=False)

mt = hl.read_matrix_table(dir_output +'/data/' + prefix +'.variantQC2.final.mt')
print(dir_output +'/data/' + prefix +'.variantQC2.final.mt')

print('After filter, %d samples and %d variants remain.' 
      % (mt.count_cols(), mt.count_rows()))

###adding final sampleQC2 and make final vcf
