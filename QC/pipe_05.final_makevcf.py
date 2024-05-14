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
hl.init(master='local[40]', spark_conf=config,  default_reference='GRCh38')
###01.1 load matrix file
mt = hl.read_matrix_table(dir_output +'/data/' + prefix +'.final.230808.mt')
print(dir_output +'/data/' + prefix +'.variantQC2.final.mt')

#hail variantQC & sampleQC
mt=hl.variant_qc(mt)
mt=hl.sample_qc(mt)

print('After Sample QC, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))
hl.export_vcf(mt, dir_output +'/data/' + prefix +'.QCed.final_230808.vcf.bgz')
hl.export_plink(mt, dir_output +'/data/' + prefix +'.QCed.final_230808', ind_id = mt.s, fam_id= mt.s)


#mt = hl.read_matrix_table(dir_output +'/data/' + prefix +'.final.230808.mt')

#for chromosome in range(1,23):
#    chromosome_mt = mt.filter_rows(mt.locus.contig == "chr" + str(chromosome))
#    chromosome_mt.write(dir_output +'/data/' + prefix +'.chr'+str(chromosome) +'.final.230808.mt', overwrite=True)
#    hl.export_vcf(chromosome_mt, dir_output +'/data/' + prefix + '.chr'+str(chromosome) +'.QCed.final_230808.vcf.bgz')


