#!/bin/sh

echo $(date)

cwas permutation_test -i CWAS_input.categorization_result.zarr \
-o_dir /output \
-s phenotype_table.tsv \
-p 45 \
-n 1000 \
-b --use_n_carrier

echo $(date)
echo "Done"
