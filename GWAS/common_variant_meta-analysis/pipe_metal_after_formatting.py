import pandas as pd
import sys

args = sys.argv
filename = args[1]

df = pd.read_csv(filename , delimiter='\t')

df[['CHR','POS','A1','A2']]=df['MarkerName'].str.split(':', expand=True)


df.to_csv(filename + ".formatting", sep='\t', index=False, na_rep='NA')

