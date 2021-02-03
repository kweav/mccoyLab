#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys

abundances_df = pd.read_csv("txi.kallisto.abundances.csv", index_col=0)
new_abundances_df = pd.DataFrame(index=abundances_df.index)
full_column_names = np.array(abundances_df.columns)
column_names = np.array([exp for exp, ct in (name.split('_') for name in abundances_df.columns)])

LA_lab_exps = np.array(['SRX669482', 'SRX669483', 'SRX669484', 'SRX669485', 'SRX669486',
               'SRX669487', 'SRX669488', 'SRX669489', 'SRX669490', 'SRX669491',
               'SRX669492', 'SRX669493', 'SRX669494', 'SRX669495', 'SRX669496',
               'SRX669497'])

H_lab_exps = np.array(['SRX3010287', 'SRX3010288', 'SRX3009974', 'SRX3009975', 'SRX3010022',
              'SRX3010023', 'SRX3010290', 'SRX3010291', 'SRX3010176', 'SRX3010177',
              'SRX3009997', 'SRX3009998', 'SRX3010228', 'SRX3010229', 'SRX3010199',
              'SRX3010200', 'SRX7517165', 'SRX7517166', 'SRX7517167', 'SRX7517168',
              'SRX3010044', 'SRX3010045', 'SRX3010068', 'SRX3010069'])

intersected, comm1_LA, comm2 = np.intersect1d(column_names, LA_lab_exps, assume_unique=True, return_indices=True)
intersected, comm1_H, comm2 = np.intersect1d(column_names, H_lab_exps, assume_unique=True, return_indices=True)

data_lab_LA = abundances_df.iloc[:, comm1_LA]
data_lab_H = abundances_df.iloc[:, comm1_H]

by_gene_variance_LA = data_lab_LA.var(axis=1)
by_gene_variance_H = data_lab_H.var(axis=1)

by_sample_variance_LA = data_lab_LA.var(axis=0)
by_sample_variance_H = data_lab_H.var(axis=0)

total_variance_LA = np.var(data_lab_LA.to_numpy().flatten())
total_variance_H = np.var(data_lab_H.to_numpy().flatten())
