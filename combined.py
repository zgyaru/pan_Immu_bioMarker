import scanpy as sc
import pandas as pd
import numpy as np
import sys
print(sys.version)


BCC = sc.read("./PAN_Cancer/anndata/pre/2019_NatMed_BCC_pre.h5ad",cache=True)
LBCL_2020 = sc.read("./PAN_Cancer/anndata/pre/2020_NatMed_LBCL_pre.h5ad",cache=True)
ccRCC = sc.read("./PAN_Cancer/anndata/pre/2021_CancerCell_ccRCC_pre.h5ad",cache=True)
TNBC = sc.read("./PAN_Cancer/anndata/pre/2021_CancerCell_TNBC_pre.h5ad",cache=True)
NSCLC_179994 = sc.read("./PAN_Cancer/anndata/pre/2021_natCancer_NSCLC_179994_pre.h5ad",cache=True)
NSCLC = sc.read("./PAN_Cancer/anndata/pre/2021_NatCancer_NSCLC_pre.h5ad",cache=True)
LBCL_2022 = sc.read("./PAN_Cancer/anndata/pre/2022_NatMed_LBCL_pre.h5ad",cache=True)
   


adata = BCC.concatenate([LBCL_2020,ccRCC,TNBC,NSCLC_179994,NSCLC,LBCL_2022],
                                  batch_categories=['2019_NatMed_BCC','2020_NatMed_LBCL','2021_CancerCell_ccRCC',
                                                    '2021_CancerCell_TNBC','2021_natCancer_NSCLC_179994','2021_NatCancer_NSCLC','2022_NatMed_LBCL'],
                                                    batch_key='batch',
                                                    join = 'outer',
                                                    fill_value=0)
#print(adata)
adata.obs = adata.obs.dropna(axis=1)
adata.var = adata.var.dropna(axis=1)







