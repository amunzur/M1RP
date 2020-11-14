import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Choose patient

pt = 'ID10'

# Import mutation data.
# Change to phylo vgs input. Documentation here: https://github.com/morrislab/phylowgs

muts = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/pyclone/%s_loci.tsv' % pt, sep='\t')

# Create first columns of input dataframe
unimuts = muts['mutation_id'].unique().tolist()

input = pd.DataFrame({'ID': ["S" + str(i) for i in np.arange(len(unimuts))], 'gene': unimuts})
