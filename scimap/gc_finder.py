import scimap as sm
import pandas as pd
import anndata as ad
import math

# change cell phenotypes based on list of cells
# returns modified anndata object
def reclassify_cells(zdata, cell_list):
   for index in zdata.obs.index:
      if index in cell_list:
         zdata.obs.loc[index, 'phenotype'] = 'False GC'
   return zdata

def cell_dist(zdata, cell1, cell2):
   cell1_x = zdata.obs.loc[cell1, 'X_centroid']
   cell1_y = zdata.obs.loc[cell1, 'Y_centroid']
   cell2_x = zdata.obs.loc[cell2, 'X_centroid']
   cell2_y = zdata.obs.loc[cell2, 'Y_centroid']
   return math.dist([cell1_x, cell1_y], [cell2_x, cell2_y])

# determine if individual GC cell is within a B cell follicle or not
# returns True if in follice, False if not in follicle
def cell_in_follicle(zdata, cell):
   # pick a radius within which to evaluate neighboring cells
   radius = 15
   # generate list of cells within radius of the cell of interest 
   neighborhood = list()
   for index in zdata.obs.index:
      if cell_dist(zdata, cell, index) <= radius:
         neighborhood.append(index)
   # determine percentage of cells within radius that are not of valid type for GC definition
   invalid_cell_counter = 0
   valid_phenotypes = ['B Cell', 'GC']
   for neighbor in neighborhood:
      if zdata.obs.loc[neighbor, 'phenotype'] not in valid_phenotypes:
         invalid_cell_counter += 1
   percent_invalid = invalid_cell_counter / len(neighborhood)
   # if neighboring cells within the radius are not >90% GC or B cell, not valid - excludes things totally outside follicle
   # if neighboring cells within the radius are not >40% GC cell, not valid - excludes "lines"
   # else it is a valid GC cell

# verify whether all GC cells in the anndata object are valid or not
# returns modified anndata object
def gc_finder(zdata):
   invalid_cells = list()
   for index in zdata.obs.index:
      if not cell_in_follicle(zdata, index):
         invalid_cells.append(index)
   return reclassify_cells(zdata, invalid_cells)
