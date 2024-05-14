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
fn_meta="/data1/mycho/WGS_AD_2011/1.data/phenotype/WGS_1824.pheno_update.230706"
config = { 'spark.driver.memory':'100g', 'spark.local.dir' : '/data1/mycho/WGS_AD_2011/temp'}
hl.init(master='local[40]', spark_conf=config,  default_reference='GRCh38')


### 3) read data

#Read SampleQC2_AT5
mt_sampleQC2_AT5 = hl.read_matrix_table(dir_output +'/data/' + prefix +'.sampleQC2_AT5.VQ1.mt')

#Read GWAS data
fn_gwas="/data1/mycho/WGS_AD_2011/4.VCF.QCed.set1-7/data/SMC_merged.maf01.rmdup.completed.anno.TOPMed_hg38_concordance.vcf.temp_concordance.m2"
hl.import_vcf(fn_gwas, force_bgz=True, reference_genome='GRCh38', contig_recoding={'1': 'chr1', '2':'chr2','3':'chr3','4':'chr4','5':'chr5','6':'chr6','7':'chr7','8':'chr8','9':'chr9','10':'chr10','11': 'chr11', '12':'chr12','13':'chr13','14':'chr14','15':'chr15','16':'chr16','17':'chr17','18':'chr18','19':'chr19','20':'chr20','21':'chr21','22':'chr22'}).write(dir_output +'/data/' + prefix +'.gwas_chip_data.temp_concordance.m2.mt', overwrite=True)

mt_GWAS = hl.read_matrix_table(dir_output +'/data/' + prefix +'.gwas_chip_data.temp_concordance.m2.mt')

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

mt_GWAS.write(dir_output +'/data/' + prefix +'.gwas_chip_data.QCed.temp.concor_matching_m2.mt', overwrite=True)

#### 3) calculate concordance
summary, samples, variants = hl.concordance(mt_sampleQC2_AT5, mt_GWAS)
samples.write( dir_output +'/data/' + prefix + '.concordacne.samples.concor_matching_m2', overwrite=True)
variants.write( dir_output +'/data/' + prefix + '.concordacne.variants.concor_matching_m2', overwrite=True)

with open(dir_output +'/data/' + prefix + '.concordacne.summary.concor_matching_m2', 'w', newline='') as f: 
    writer = csv.writer(f) 
    writer.writerow(summary) 

samples = hl.read_table(dir_output +'/data/' + prefix + '.concordacne.samples.concor_matching_m2')
variants = hl.read_table(dir_output +'/data/' + prefix + '.concordacne.variants.concor_matching_m2')

# non-ref disconcordance
sample_non_ref_con = (samples.concordance[3][3]+samples.concordance[4][4]) / (samples.concordance[1][3]+samples.concordance[1][4]+
                                                                             samples.concordance[2][3]+samples.concordance[2][4]+
                                                                             samples.concordance[3][3]+samples.concordance[3][4]+
                                                                             samples.concordance[4][3]+samples.concordance[4][4])

samples = samples.annotate(non_ref_concordance = sample_non_ref_con)

samples_filtered_non_ref_concordance = samples.filter(samples.non_ref_concordance < 0.9)

print('After filter, %d samples filtered.' % (samples_filtered_non_ref_concordance.count()))
print(samples_filtered_non_ref_concordance.s.show())

samples_filtered_non_ref_concordance.export(dir_output + '/filtered_sample/' + prefix + '_non_ref_concordance.concor_matching_m2', delimiter='\t')




