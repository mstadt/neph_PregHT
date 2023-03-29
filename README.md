# About
This version of the nephron model was used for the paper 
[Stadt et al. Effect of pregnancy and hypertension on kidney function in female rats: Modeling and functional implications](https://www.biorxiv.org/content/10.1101/2022.12.15.520674v1.abstract)

# Directory information
### model_code/
Model code files used to generate simulations. Simulation results given in **plot_results/data_files/**

### plot_results/
**data_files/** contains data files created from model simulation results

# Instructions for model simulations
### This is for running simulations using files under **model_code/**. Move model output to directory with files from **plot_results/** to make figures.
To run the parallel simulation code use command: **python3 parallel_simulate.py --sex [option] --species [option] --type [option] --diabetes [option] --inhibition [option] --pregnant [option] --obese [option]**

The options here are:

sex: **Male, Female** (required);

species: **human, rat, mouse** (required);

type: **superficial, multiple** (required);

diabetes: **Severe, Moderate, Non** (optional, default: Non);

pregnant: **mid, late** (optional, default: non, only for female rat);

obese: **Y, N** (optional, default: N);

inhibition: **ACE, SGLT2, NHE3-50, NHE3-80, NKCC2-70, NKCC2-100, NCC-70, NCC-100, ENaC-70, ENaC-100, SNB-70, SNB-100, HKA-100, HKApreg-100** (optional, default: None).

unx: **N, Y** (optional, default: N)

Notes:
* Human only have ACE and SGLT2 inhibition cases. The others are for rats.
* pregnancy: only has been characterized for normal pregnant rat nephron models at this time (i.e., not done for humans and for diabetes)
* mouse model is not finished for any case

### Understanding output

All the output files' names are in following structure: 'sex_species_segment_concentration/flow_of_solute_in_compartment.txt'. 

Here is an example: female_rat_ccd_con_of_Cl_in_Bath.txt. It contains interstitial concentration of Chloride along cortical collecting duct in female rat.

Another example: male_hum_pt_flow_of_Na_in_Lumen.txt. It contains luminal flow of Sodium along proximal convolute tubule in male human.

These results are scaled per nephron.

The unit of concentration from outputs is **mmol/L (mM)**.

The unit of volume is **nl/min**.

The unit of flow is **pmol/min**.
