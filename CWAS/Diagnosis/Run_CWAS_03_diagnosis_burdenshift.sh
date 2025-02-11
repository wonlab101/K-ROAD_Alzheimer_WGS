#!/bin/sh

echo $(date)

cwas burden_shift \
-i CWAS_input.burden_test.txt \
-b CWAS_input.binom_pvals.txt.gz \
-c_info CWAS_input.category_info.txt \
-c_count CWAS_input.category_counts.txt \
-o_dir /output \
-c_cutoff 9


echo $(date)
echo "Done"
