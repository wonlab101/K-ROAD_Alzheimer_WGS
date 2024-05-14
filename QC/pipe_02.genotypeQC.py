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

config = { 'spark.driver.memory':'200g', 'spark.local.dir' : '/data1/mycho/WGS_AD_2011/temp'}
hl.init(master='local[40]', spark_conf=config,  default_reference='GRCh38')

###01.1 load matrix file
mt = hl.read_matrix_table(dir_output +'/data/' + prefix +'.VQ1.mt')
print('import_VQ1 matrix : After filter, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))

###01.2 pre-processing
#add meta data
table= (hl.import_table(fn_meta, impute=True).key_by('FID'))
table.describe(widget=True)
table.show(5)

mt = mt.annotate_cols(pheno = table[mt.s])
print('add meta data : After Sample QC, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))

samples_to_remove = {"WGS_0345"}
set_to_remove = hl.literal(samples_to_remove)
mt = mt.filter_cols(~set_to_remove.contains(mt['s']))
print('remove WGS_0345 : After Sample QC, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))


###01.3 Genotype QC
mt = mt.annotate_entries(AB = mt.AD[1] / hl.sum(mt.AD))

#GQ
mt = mt.filter_entries((mt.GQ>=20))
print('GQ>=20 : After filter, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))

#AB
mt = mt.filter_entries((mt.GT.is_hom_ref()) |
                       (mt.GT.is_het()&(mt.AB >= 0.2))& (mt.GT.is_het()&(mt.AB <= 0.8)) |
                       (mt.GT.is_hom_var()&(mt.AB >= 0.9)), keep=True)
print('AB : After filter, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))

#DP
mt_auto = filter_to_autosomes(mt)
mt_auto = mt_auto.filter_entries((mt_auto.DP>=10) & (mt_auto.DP<=200))
print('DP_AUTO : After filter, %d samples and %d variants remain.' % (mt_auto.count_cols(), mt_auto.count_rows()))

mt_X = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX')])
print('DP_X :After filter, %d samples and %d variants remain.' % (mt_X.count_cols(), mt_X.count_rows()))

mt_X_female = mt_X.filter_cols(mt_X.pheno.sex == 2 )
mt_X_female = mt_X_female.filter_entries((mt_X_female.DP>=10) & (mt_X_female.DP<=200))
print('DP_X_female : After filter, %d samples and %d variants remain.' % (mt_X_female.count_cols(), mt_X_female.count_rows()))

mt_X_male = mt_X.filter_cols(mt_X.pheno.sex == 1 )
mt_X_male = mt_X_male.filter_entries((mt_X_male.DP>=5) & (mt_X_male.DP<=200))
print('DP_X_male : After filter, %d samples and %d variants remain.' % (mt_X_male.count_cols(), mt_X_male.count_rows()))

mt_X_merged = mt_X_female.union_cols(mt_X_male)
print('DP_X : After filter, %d samples and %d variants remain.' % (mt_X_merged.count_cols(), mt_X_merged.count_rows()))

def align_mt2_cols_to_mt1(mt1, mt2):
    mt1 = mt1.add_col_index()
    mt2 = mt2.add_col_index()
    new_col_order = mt2.index_cols(mt1.col_key).col_idx.collect()
    new_col_order_seq = list(new_col_order)  # Convert the list to a tuple
    return mt2.choose_cols(new_col_order_seq)

mt_X_merged_sort = align_mt2_cols_to_mt1(mt_auto, mt_X_merged)

mt_Y = hl.filter_intervals(mt, [hl.parse_locus_interval('chrY')])
print('DP_Y : After filter, %d samples and %d variants remain.' % (mt_Y.count_cols(), mt_Y.count_rows()))

mt_Y = mt_Y.filter_entries((mt_Y.DP>=5) & (mt_Y.DP<=200))

#Merge
all_datasets = [mt_auto, mt_X_merged_sort, mt_Y]
mt_genotypeQCed = hl.MatrixTable.union_rows(*all_datasets)
print('DP_MERGE : After filter, %d samples and %d variants remain.' % (mt_genotypeQCed.count_cols(), mt_genotypeQCed.count_rows()))

#hail variantQC & sampleQC
mt_genotypeQCed=hl.sample_qc(mt_genotypeQCed, name='sample_qc')
mt_genotypeQCed=hl.variant_qc(mt_genotypeQCed, name='variant_qc')

####save 
mt.write(dir_output +'/data/' + prefix +'.VQ1.GQ.mt', overwrite=True)


##PLOT_beforeQC
mt = hl.read_matrix_table(dir_output +'/data/' + prefix +'.VQ1.mt')
print('After filter, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))

mt=hl.sample_qc(mt)
mt=hl.variant_qc(mt)

call_ht = mt.sample_qc.n_called

print(dir_output + 'data_qc/' + prefix  +'_raw_qc_n_call')
call_ht.export(dir_output + 'data_qc/' + prefix  +'_raw_qc_n_call')

df_sum = pd.read_csv(dir_output + 'data_qc/' + prefix  +'_raw_qc_n_call', sep='\t')

df_sum.loc['total']= df_sum.sum(axis=0)

df_total = df_sum.loc['total']
df_total

mt = hl.read_matrix_table(dir_output +'/data/' + prefix +'.VQ1.mt')
print('After filter, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))

p = hl.plot.histogram(mt.GQ, title = 'GQ pre', )

p.renderers.extend(
        [Span(location=20, dimension='height', line_color='red', line_width=1),
         ])




output_file(dir_output + 'plot/' + prefix + '_GQ_pre.html')
save(p)

p = hl.plot.histogram(mt.DP, title = 'DP pre 1',)

p.renderers.extend(
        [Span(location=10, dimension='height', line_color='red', line_width=1),
         Span(location=200, dimension='height', line_color='red', line_width=1),
         ])



output_file(dir_output + 'plot/' + prefix + '_DP_pre_1.html')
save(p)

p = hl.plot.histogram(mt.DP, title = 'DP pre 2', range=(0,300),)

p.renderers.extend(
        [Span(location=10, dimension='height', line_color='red', line_width=1),
         Span(location=200, dimension='height', line_color='red', line_width=1),
         ])



output_file(dir_output + 'plot/' + prefix + '_DP_pre_2.html')
save(p)

mt = mt.annotate_entries(AB = mt.AD[1] / hl.sum(mt.AD))

mt_hom_ref = mt.filter_entries((mt.GT.is_hom_ref()), keep=True)

p = hl.plot.histogram(mt_hom_ref.AB, title = 'AB - hom_ref pre', range=(0,1),)

p.renderers.extend(
        [Span(location=0.2, dimension='height', line_color='red', line_width=1),
         ])



output_file(dir_output + 'plot/' + prefix + '_AB_hom_ref_pre.html')
save(p)

mt_het = mt.filter_entries((mt.GT.is_het()), keep=True)

p = hl.plot.histogram(mt_het.AB, title = 'AB - het pre', range=(0,1),)

p.renderers.extend(
        [Span(location=0.2, dimension='height', line_color='red', line_width=1),
         ])
p.renderers.extend(
        [Span(location=0.8, dimension='height', line_color='red', line_width=1),
         ])


output_file(dir_output + 'plot/' + prefix + '_AB_het_pre.html')
save(p)

mt_hom_var = mt.filter_entries((mt.GT.is_hom_var()), keep=True)

p = hl.plot.histogram(mt_hom_var.AB, title = 'AB - hom_var pre', range=(0,1),)

p.renderers.extend(
        [Span(location=0.9, dimension='height', line_color='red', line_width=1),
         ])



output_file(dir_output + 'plot/' + prefix + '_AB_hom_var_pre.html')
save(p)













##PLOT_afterQC
mt = hl.read_matrix_table(dir_output +'/data/' + prefix + '.VQ1.GQ.mt')
print('After filter, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))

call_ht = mt.sample_qc.n_called

print(dir_output + 'data_qc/' + prefix  +'_genotype_qc_n_call')
call_ht.export(dir_output + 'data_qc/' + prefix  +'_genotype_qc_n_call')

df_sum = pd.read_csv(dir_output + 'data_qc/' + prefix  +'_genotype_qc_n_call', sep='\t')
df_sum.loc['total']= df_sum.sum(axis=0)





df_total = df_sum.loc['total']
df_total

p = hl.plot.histogram(mt.GQ, title = 'GQ after', )

p.renderers.extend(
        [Span(location=20, dimension='height', line_color='red', line_width=1),
         ])



output_file(dir_output + 'plot/' + prefix + '_GQ_after.html')
save(p)

p = hl.plot.histogram(mt.DP, title = 'DP after', range=(0,300),)

p.renderers.extend(
        [Span(location=10, dimension='height', line_color='red', line_width=1),
         Span(location=200, dimension='height', line_color='red', line_width=1),
         ])



output_file(dir_output + 'plot/' + prefix + '_DP_after_2.html')
save(p)

mt_hom_ref = mt.filter_entries((mt.GT.is_hom_ref()), keep=True)

p = hl.plot.histogram(mt_hom_ref.AB, title = 'AB - hom_ref after', range=(0,1),)

p.renderers.extend(
        [Span(location=0.2, dimension='height', line_color='red', line_width=1),
         ])



output_file(dir_output + 'plot/' + prefix + '_AB_hom_ref_after.html')
save(p)

mt_het = mt.filter_entries((mt.GT.is_het()), keep=True)

p = hl.plot.histogram(mt_het.AB, title = 'AB - het after', range=(0,1),)

p.renderers.extend(
        [Span(location=0.2, dimension='height', line_color='red', line_width=1),
         ])
p.renderers.extend(
        [Span(location=0.8, dimension='height', line_color='red', line_width=1),
         ])



output_file(dir_output + 'plot/' + prefix + '_AB_het_after.html')
save(p)

mt_hom_var = mt.filter_entries((mt.GT.is_hom_var()), keep=True)

p = hl.plot.histogram(mt_hom_var.AB, title = 'AB - hom_var after', range=(0,1),)

p.renderers.extend(
        [Span(location=0.9, dimension='height', line_color='red', line_width=1),
         ])



output_file(dir_output + 'plot/' + prefix + '_AB_hom_var_after.html')
save(p)

