# Probabilistic Coseismic Displacement Hazard Model
Framework for calculating coastal coseismic displacement using national seismic hazard model solutions 
in the Wellington Region. This repository accompanies the manuscript "A probabilistic model for coseismic vertical 
displacement hazard in coastal settings" for publication in Geosphere by Delano et al. (in prep). 

## Scripts
Scripts are all in the `scripts` directory. 
Data are in the `data` directory.

NOTE: each script has a series of "USER INPUTS" at the top of each script. 
Check the user inputs before running and follow comment instructions. 
The full range of results/figures requires the user to change some of these variables before running the scripts. 
Examples for changing variables: 
- desired calculation extents/locations
- figures (True/False; saves disk space and run time), 
- grid size (file size/run time)
- output file directory (to avoid overwriting alternative trials)

The main scripts should be run in the following 
order:
A. Prep crustal fault meshes and greens functions. The results are put into a new folder in the crustal directory
1. `01_crustal_subset_fault_sections_by_mesh.py` This script extracts the fault sections from the NSHM that are 
   located along meshes in the Wellington area. 
2. `02_crustal_discretize_remeshed.py` Reads in fault traces from the NSHM, turns them into rectangular 
   patches with metadata, and matches each patch to triangles in the finer mesh. 
3. choose one or more Green's function scripts:
   1. `03_crustal_discretized_grid_gfs.py` Calculates greens functions from the mesh on a 
      regional grid.
   2. `03_crustal_discretized_coast_gfs.py` Calculates greens functions from the mesh 
      on points along the coastline.
   3. `03_crustal_discretized_sites_gfs.py` Calculates greens functions from the mesh at 
      specific coordinates.
   
B. subduction zone mesh and greens function file prep. 
4. `04_sz_discretize.py` Reads in fault traces from the NSHM, turns them into 
   rectangular patches with metadata, and matches each patch to triangles in the mesh. 
5. Choose one or more greens functions for displacements:
   1. `05_sz_discretized_gfs_grid.py` Calculates greens functions on a grid
   2. `05_sz_discretized_gfs_coast.py` Calculates greens functions along the coastline shapefile
   3. `05_sz_discretized_gfs_sites.py` Calculates greens functions at specific coordinates

C. calculate scenario displacements and branch probabilities
6. `06_run_displacements_and_probabilities.py` Calculates displacements for each scenario and writes 
   pickle to results directory. Can run multiple branches from a fault model at the same time for same gf type or a 
   single branch (to save time). There are many possible user inputs, read prompts before choosing variables. 

D. calculate weighted mean probabilities for all branches in a fault model
7. `07_run_aggregate_weighted_branches.py` Calculates weighted mean probabilities for all branches in a 
   fault model (e.g., SZ unmodified dip) and writes pickle to results directory. Only works for sites (not grid or coast). 
   Can aggregate a subduction and crustal fault model together (i.e. paired fault model) because it just adds all 
   displacements/probabilities together. Many user inputs, read comments to decide how to set up. 

E. Make additional plots (optional, useful for comparisons)
8. `08_run_compare_fault_models.py` Additional plotting scripts to compare different branches or fault models.
   Many user input toggle options based on what you want to plot together. 
   

