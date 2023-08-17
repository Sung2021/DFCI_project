# scrublet application
# 23.08.17

# import library
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr

# define function
def run_scrublet(input_file, output_file):
    df = pd.read_csv(input_file, header=0, index_col=0)
    adata = sc.AnnData(df)
    sc.external.pp.scrublet(adata)
    adata.obs.to_csv(output_file)

# usage function
input_file = 'Takeda.NK.30293.count.mtx.csv'
output_file = 'Takeda.NK.30293.scr.csv'
run_scrublet(input_file, output_file)

input_file = 'Takeda.NK.30298.count.mtx.csv'
output_file = 'Takeda.NK.30298.scr.csv'
run_scrublet(input_file, output_file)

