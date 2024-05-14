import sys
import pandas as pd
import numpy as np


# Access command-line arguments
args = sys.argv
filename = args[1]
OR=args[2]
column_names = args[3:]

# Read the DataFrame from file
df = pd.read_csv(filename, delimiter='\t')

print(args)
print(column_names)
print(df.columns)

# Get specific columns by names
columns = df.loc[:, column_names]

# Rename the column names
new_column_names = ['CHR', 'POS', 'EFF', 'NONEFF', 'PVAL', 'OR', 'SE', 'N']
columns = columns.rename(columns=dict(zip(column_names, new_column_names)))

if OR == "OR":
        columns['BETA'] = np.log(columns['OR'])
else:
        columns['BETA'] = columns['OR']

print(columns)

columns['ID'] = columns.apply(lambda row: ':'.join([str(row['CHR']), str(row['POS']), ':'.join(sorted([str(row['EFF']), str(row['NONEFF'])]))]), axis=1)


columns.to_csv(filename +".metalinput" , sep='\t', index=False, na_rep='NA')
