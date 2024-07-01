# MCMICRO Setup  
Full documentation available at: https://mcmicro.org/  
  
## In Brief  
Set up a conda environment (such as by using Anaconda) called `mcmicro`.  
Activate the environment and install the required packages and programs for the Nextflow pipeline.  
The MCMICRO website (above) provides detailed instructions under Tutotial --> Installation.  
  
# Running MCMICRO
Run the bash script in this repository, but modify the paths to fit your filesystem.  
The script contains comments explaining directory structure, as does the tutorial on the mcmicro website.  
The relevant output will be in subdirectories such as `registration/`, `qc/s3seg/`, and `quantification/`.  
- 'registration/' contains the overlaid (registered) image.
- 'qc/s3seg/' contains the cell masks from segmentation, which should be overlaid and viewed using FIJI/ImageJ.
- 'quantification/' contains a .csv describing the morphology of each cell and the quantified amount of marker in each cell.
  
# SCIMAP Setup
Full documentatino available at: https://scimap.xyz/  
  
## In Brief
Set up a conda environment (such as by using Anaconda) called `scimap`.  
Activate the environment and install the requires packages, which are specified in the installation instructions on the SCIMAP website.  
  
# Running SCIMAP
Relevant scripts are in the scimap/ directory in this repository.  
Edit file paths as appropriate for your filesystem.  
Some parameters might need to be changed. See comments in code for parameters when using gc_finder.py as well.  
**Strongly recommend** using scimap_gc_finder_analysis.py instead of the separate `scimap_analysis.py` and `gc_finder.py`.  
N.B. `phenotype_workflow.csv` is project-specific and needs to be edited for given markers and phenotypes of interest.  

# Citations
## MCMICRO
Schapiro, D., Sokolov, A., Yapp, C. et al. MCMICRO: a scalable, modular image-processing pipeline for multiplexed tissue imaging. Nat Methods 19, 311–315 (2022). https://doi.org/10.1038/s41592-021-01308-y  

## SCIMAP
Nirmal et al., (2024). SCIMAP: A Python Toolkit for Integrated Spatial Analysis of Multiplexed Imaging Data. Journal of Open Source Software, 9(97), 6604, https://doi.org/10.21105/joss.06604
