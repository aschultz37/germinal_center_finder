import scimap as sm
import pandas as pd
import anndata as ad
from matplotlib.path import Path
import math
from datetime import datetime                #DEBUG

# change cell phenotypes based on list of cells
# returns modified anndata object
def reclassify_cells(zdata: ad.AnnData, cell_list: list, new_cell_type: str) -> ad.AnnData:
   zdata.obs.phenotype = zdata.obs.phenotype.cat.add_categories(new_cell_type)
   for index in zdata.obs.index:
      if index in cell_list:
         zdata.obs.loc[index, 'phenotype'] = new_cell_type
   return zdata

# calculates the Euclidian distance between two cells
# returns distance as a float
def cell_distance(zdata: ad.AnnData, cell1: str, cell2: str) -> float:
   cell1_x = zdata.obs.loc[cell1, 'X_centroid']
   cell1_y = zdata.obs.loc[cell1, 'Y_centroid']
   cell2_x = zdata.obs.loc[cell2, 'X_centroid']
   cell2_y = zdata.obs.loc[cell2, 'Y_centroid']
   return math.dist([cell1_x, cell1_y], [cell2_x, cell2_y])

# generates points defining n equal arcs given a center point (x, y) tuple
# returns a list of the points
def generate_arcs(center: tuple, radius: float, num_arcs: int) -> list:
   arcs = list()
   angles = list(math.radians(45)*i for i in range(0, num_arcs))
   for angle in angles:
      x_coord = center[0] + (radius * math.cos(angle))
      y_coord = center[1] + (radius * math.sin(angle))
      arcs.append((x_coord, y_coord))
   return arcs

# checks if a point falls within the arc of a circle, estimated as a triangle polygon
# returns True if within arc, False if not
def check_within_arc(center: ad.AnnData, start_point: tuple, end_point: tuple, point_of_interest: tuple) -> bool:
   triangle = Path([center, start_point, end_point])
   return triangle.contains_point(point_of_interest)

# determines the distribution of cells in a list around a central cell within a given radius
# returns ... something that represents that, TBD, maybe just a Boolean for if valid or not
def cell_distribution(zdata: ad.AnnData, center_cell: str, invalid_cell_list: list, valid_cell_list: list, radius: float) -> bool:
   num_arcs = 8
   consecutive_arcs = 3
   validity_cutoff = 0.50
   center = (zdata.obs.loc[center_cell, 'X_centroid'], zdata.obs.loc[center_cell, 'Y_centroid'])
   arcs = generate_arcs(center, radius, num_arcs)
   arcs_invalid_counter = [0] * num_arcs
   arcs_valid_counter = [0] * num_arcs
   # determine # invalid cells in each section
   for invalid_cell in invalid_cell_list:
      cell_coords = (zdata.obs.loc[invalid_cell, 'X_centroid'], zdata.obs.loc[invalid_cell, 'Y_centroid'])
      for i in range(0, num_arcs):
         if check_within_arc(center, arcs[i], arcs[(i+1) % num_arcs], cell_coords):
            arcs_invalid_counter[i] += 1
   # determine # valid cells in each section
   for valid_cell in valid_cell_list:
      cell_coords = (zdata.obs.loc[valid_cell, 'X_centroid'], zdata.obs.loc[valid_cell, 'Y_centroid'])
      for i in range(0, num_arcs):
         if check_within_arc(center, arcs[i], arcs[(i+1) % num_arcs], cell_coords):
            arcs_valid_counter[i] += 1
   # determine % valid in each section
   arcs_total_counter = [arcs_valid_counter[i] + arcs_invalid_counter[i] for i in range(0, num_arcs)]
   for i in range(0, len(arcs_total_counter)):
      if arcs_total_counter[i] == 0:
         arcs_total_counter[i] = -1
   arcs_pct_valid = [arcs_valid_counter[i] / arcs_total_counter[i] for i in range(0, num_arcs)]
   # if certain number of consecutive arcs are primarily valid, it is OK, probably edge of a GC
   arcs_valid_bool = [arcs_pct_valid[i] > validity_cutoff for i in range(0, num_arcs)]
   for i in range(0, num_arcs):
      consec_valid = 0
      for j in range(i, i+consecutive_arcs):
         consec_valid += arcs_valid_bool[i%num_arcs]
      if consec_valid >= consecutive_arcs:
         return True
   return False

# determine if individual GC cell is within a B cell follicle or not
# returns True if in follice, False if not in follicle
def cell_in_follicle(zdata: ad.AnnData, radius: int, follicle_cutoff_pct: float, edge_cutoff_pct: float, cell: str) -> bool:
   # pick a radius within which to evaluate neighboring cells and cutoffs for follicle composition
   #radius = 15
   #follicle_cutoff_pct = 0.20
   #edge_cutoff_pct = 0.20
   # generate list of cells within radius of the cell of interest 
   neighborhood = list()
   # NOTE: redo this by first making a new col. of avg coord, sorting on it, then only checking cells that have avg coord w/in r of the center
   #for index in zdata.obs.index:
   #   if cell_distance(zdata, cell, index) <= radius:
   #      neighborhood.append(index)
   zdata.obs['avg_coord'] = zdata.obs[['X_centroid', 'Y_centroid']].mean(axis=1)
   ac_within_range_list = zdata.obs.index[abs(zdata.obs['avg_coord'] - zdata.obs.loc[cell, 'avg_coord']) < (0.8 * radius)].tolist()
   for index in ac_within_range_list:
      if cell_distance(zdata, cell, index) <= radius:
         neighborhood.append(index)
   # determine percentage of cells within radius that are not of valid type for GC definition
   invalid_cells = list()
   valid_cells = list()
   valid_phenotypes = ['B Cell', 'GC']
   for neighbor in neighborhood:
      if zdata.obs.loc[neighbor, 'phenotype'] not in valid_phenotypes:
         invalid_cells.append(neighbor)
      else:
         valid_cells.append(neighbor)
   percent_invalid = len(invalid_cells) / len(neighborhood)
   # if neighboring cells within the radius are not >90% GC or B cell, not valid
   # excludes things totally outside follicle
   if percent_invalid > follicle_cutoff_pct:
      return False
   # if neighboring cells within the radius are not >40% GC cell, not valid - excludes "lines"
      # actually, this needs to check across what arcs the invalid cells are distributed
   if not cell_distribution(zdata, cell, invalid_cells, valid_cells, radius):
      return False
   # else it meets criteria and is a valid GC cell
   return True

# verify whether all GC cells in the anndata object are valid, changes phenotype if not
# returns modified anndata object
def gc_finder(zdata: ad.AnnData, new_cell_type: str, follicle_cutoff_pct: float, edge_cutoff_pct:float, radius: int) -> ad.AnnData:
   wdata = zdata.copy()
   gc_type_cell = ['GC']
   invalid_cells = list()
   cycle_counter = 0                                                                   #DEBUG
   for index in wdata.obs.index:
      if cycle_counter % 5000 == 0:                                                    #DEBUG
         print(str(datetime.now()) + ': On cell ' + str(cycle_counter) + " of " + str(len(wdata.obs.index)))   #DEBUG
      if wdata.obs.loc[index, 'phenotype'] in gc_type_cell:
         if not cell_in_follicle(wdata, radius, follicle_cutoff_pct, edge_cutoff_pct, index):
            invalid_cells.append(index)
      cycle_counter += 1                                                               #DEBUG
   return reclassify_cells(wdata, invalid_cells, new_cell_type)
