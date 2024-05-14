
import sys
import pandas as pd
import numpy as np

# Access command-line arguments
args = sys.argv
fn1 = args[1]
fn2 = args[2]
fn3 = args[3]
output1 = args[4]
output2 = args[5]
output3 = args[6]

# Read the DataFrame from file
df1 = pd.read_csv(fn1, delimiter='\t')
df2 = pd.read_csv(fn2, delimiter='\t')
df3 = pd.read_csv(fn3, delimiter='\t')
print("____",fn1,":  ",df1.shape)
print("____",fn2,":  ",df2.shape)
print("____",fn3,":  ",df3.shape)

# Identify the ID column in each data frame
id_column = 'ID'

# Extract the common IDs
common_ids = set(df1[id_column]).intersection(set(df2[id_column])).intersection(set(df3[id_column]))
num_common_ids = len(common_ids)
print("____Number of common IDs:", num_common_ids)

# Filter data frames to include only common IDs
df1_common = df1[df1[id_column].isin(common_ids)]
df2_common = df2[df2[id_column].isin(common_ids)]
df3_common = df3[df3[id_column].isin(common_ids)]

#save each file separately
df1_common.to_csv(output1, sep='\t', index=False)
df2_common.to_csv(output2, sep='\t', index=False)
df3_common.to_csv(output3, sep='\t', index=False)


