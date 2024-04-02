import sys
import os
import pandas as pd

bioproject=sys.argv[1]
filenames =sys.argv[2:]

counts = pd.read_csv(filenames[0], sep='\t', index_col=0)["NumReads"]

for f in filenames[1:]:
    counts = pd.concat([counts, pd.read_csv(f, sep='\t', index_col=0)["NumReads"]], axis=1)
counts = round(counts).astype(int)
counts.columns = list(map(lambda f: f.split("/")[1].rstrip("_quant"), filenames[1:]))
counts.to_csv(f"../data/datasets/{bioproject}/counts.tsv", sep='\t')
