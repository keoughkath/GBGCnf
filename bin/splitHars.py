import pandas as pd
import sys, os

inhars = sys.argv[1]
mafdir = sys.argv[2]
outdir = sys.argv[3]

hars = pd.read_table(inhars,
                    header=None, names=['chrom','start','end','pvalue']).sort_values(by='pvalue').reset_index(drop=True)
hars['name'] = hars.index
hars = hars[['chrom','start','end','name']]

# these are split into 10Mb chunks, low end is specified in file name

for maf in os.listdir(mafdir):
    chrom = maf.split('.')[0]
    low = int(maf.split('.')[1])
    high = low + 10000000
    relevant_hars = hars.query('chrom == @chrom and start >= @low and end <= @high')
    if len(relevant_hars) > 0:
        relevant_hars.to_csv(f'{outdir}/{maf.split(".maf")[0]}.bed',
                            header=None, index=None, sep='\t')