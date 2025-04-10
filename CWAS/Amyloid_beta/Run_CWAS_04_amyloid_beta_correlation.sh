#!/bin/sh

echo $(date)

input_dir=/input
output_dir=/output

cwas correlation \
-i $output_dir/categorization/CWAS_input.categorization_result.zarr \
-cm sample \
-o_dir $output_dir/correlation \
-p 38 \
-c_info $output_dir/CWAS_input.category_info.txt

echo $(date)
echo "Done"

