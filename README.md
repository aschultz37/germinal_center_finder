# mcmicro Setup  
Full documentation available at: https://mcmicro.org/  
  
## In Brief  
Set up a conda environment (such as by using Anaconda) for mcmicro.  
Activate the environment and install the required packages and programs for the Nextflow pipeline.  
The mcmicro website (above) provides detailed instructions under Tutotial --> Installation.  
  
# Running mcmicro
Run the bash script in this repository, but modify the paths to fit your filesystem.  
The script contains comments explaining directory structure, as does the tutorial on the mcmicro website.  
The relevant output will be in subdirectories such as registration/, qc/s3seg/, and quantification/.  
- Registration contains the overlaid (registered) image.
- S3seg contains the cell masks from segmentation, which should be overlaid and viewed using FIJI/ImageJ.
- Quantification contains a .csv describing the morphology of each cell and the quantified amount of marker in each cell.

# Running scimap
Relevant scripts are in the scimap/ directory in this repository.  
Edit file paths as appropriate for your filesystem.  
Some parameters might need to be changed. See comments in code for parameters when using gc_finder.py as well.
