# About
This is the version of sex-specific mathematical models for epithelial transport along the nephron implemented in Python 3 used in [Stadt et al. Effect of pregnancy and hypertension on kidney function in female rats: Modeling and functional implications, under review](https://www.biorxiv.org/content/10.1101/2022.12.15.520674v1.abstract). Related research papers are listed below. Please cite appropriately. To improve computation time, parallel computation of different types of nephron is implemented for the multiple nephron model. 

# Instructions
To run the parallel simulation code use command: **python3 parallel_simulate.py --sex [option] --species [option] --type [option] --diabetes [option] --inhibition [option] --pregnant [option]**

The options here are:

sex: **Male, Female** (required);

species: **human, rat, mouse** (required);

type: **superficial, multiple** (required);

diabetes: **Severe, Moderate, Non** (optional, default: Non);

pregnant: **mid, late** (optional, default: non, only for female rat);

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

**/plot/** contains scripts for plotting output


