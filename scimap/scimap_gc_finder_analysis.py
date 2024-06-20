import scimap as sm
import pandas as pd
import anndata as ad
import os

# Process and Phenotype mcmicro Output
## displays spatial plot of cells colored by option (phenotype, ROI, etc.)
def plot(adata, color):
   sm.pl.spatial_scatterPlot(adata, colorBy=color, figsize=(7,7), s=0.7, fontsize=5, catCmap='Set1')

## set working directory
os.chdir('/home/austin/mcmicro/L2-Section2')

## import mcmicro output to scimap
feature_table_path = ["quantification/L2-Section2--unmicst_cell.csv"]

adata = sm.pp.mcmicro_to_scimap(feature_table_path)

## set auto gates
adata = sm.pp.rescale(adata)

## phenotype cells
phenotype = pd.read_csv('scimap/phenotype_workflow.csv')
phenotype.style.format(na_rep='')

adata = sm.tl.phenotype_cells(adata, phenotype=phenotype, label="phenotype")

## save anndata object as h5ad
adata.write('scimap/L2-Section2_phenotyped.h5ad')

# Perform Spatial Clustering
bdata = adata.copy()
bdata = sm.tl.spatial_count(bdata)
bdata = sm.tl.spatial_cluster(bdata, k=9)

sm.pl.spatial_scatterPlot(bdata, colorBy='spatial_kmeans', figsize=(7,7), s=0.5, fontsize=5, fileName='scimap/spatial_kmeans.png')

# Find GC ROIs
## cluster 6 is GC visually
gc_cluster = bdata.obs[bdata.obs['spatial_kmeans'] == '6']
