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

config = { 'spark.driver.memory':'100g', 'spark.local.dir' : '/data1/mycho/WGS_AD_2011/temp'}
hl.init(master='local[16]', spark_conf=config,  default_reference='GRCh38')

###01.1 load VCF file

hl.import_vcf(vcf_dir + fn_vcf, force_bgz=True, reference_genome='GRCh38').write(dir_output +'/data/' + prefix +'.mt', overwrite=True)
print(dir_output +'/data/' + prefix +'.mt')

mt = hl.read_matrix_table(dir_output +'/data/' + prefix +'.mt')
print('After filter, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))

#samples_to_remove = {"WGS_0345"}


#set_to_remove = hl.literal(samples_to_remove)
#mt = mt.filter_cols(~set_to_remove.contains(mt['s']))

#print('After Sample QC, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))



## 02-01) Split mult, VQSR, LCR

mt = hl.split_multi_hts(mt, permit_shuffle=True)
mt = mt.key_rows_by(**hl.min_rep(mt.locus, mt.alleles))

print('Multi_hts : After filter, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))


mt = mt.filter_rows(hl.len(mt.filters) == 0)
print('VQSR : After filter, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))


bed_LCR = hl.import_bed(fn_LCR)
mt = mt.annotate_rows(LCR_region = bed_LCR[mt.locus])
mt = mt.filter_rows(hl.is_defined(mt.LCR_region), keep=False)
print('LCR : After filter, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))


#hail variantQC & sampleQC
mt=hl.sample_qc(mt)
mt=hl.variant_qc(mt)

## For Removing the variants that count Zero(0) becasue of sample removing
HARD_CUTOFF = 0

mt = mt.filter_rows(mt.variant_qc.AC[1] > HARD_CUTOFF)
print('remove Zero : After filter, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))

mt.write(dir_output +'/data/' + prefix +'.VQ1.mt', overwrite=True)


