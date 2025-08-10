# Afterlip-stress-change-inversion-
This repository contains MATLAB programs by A. Meneses-Gutierrez and T. Saito for shear stress change inversion in a specified direction, developed for the study “Linking Coseismic Slip and Afterslip in Intraplate Earthquakes: A Case Study of the 2016 Central Tottori Earthquake, Japan” (2025JB031677).

The inversion requires as input (*.mat files):
Observed surface displacements
The geometry of the fault plane for which stress changes are estimated
A set of basis functions to parameterize the stress drop distribution
 
This software requires the code for calculating displacements, 
strains, and stresses associated with triangular dislocations (TDs)
in an elastic half-space, developed by Nikkhoo and Walter (2015).
We downloaded the code from:
https://volcanodeformation.com/software
If you use the software in your work, please also cite the original reference:
[Nikkhoo, M., & Walter, T. R. (2015). Triangular dislocation: an analytical,
 artefact-free solution. Geophysical Journal International, 201(2), 1117–1139.
  https://doi.org/10.1093/gji/ggv035]

This code is provided as is and may contain errors or unexpected behavior.
Although efforts have been made to ensure accuracy and functionality,
 it has not undergone comprehensive testing.
Use it at your own risk. Users should verify the outputs and adapt the
 code as needed for their particular requirements.

######################################################################
Input files (.mat):
Displacement observations:read_2016Tottori
Basis functions information: area_pos_B01
Fault geometry and mesh information: faultP_par01

How to run:
slip_vector_b01.m estimates normal vectors and slip vectors
surf_dis_b01.m estimates the surface displacement at the GNSS stations
due to 1 m of slip in a direction consistent with the uniform coseismic
slip
meshresp_b01.m estimate traction due to unit slip for each subfault
stress_base_b01.m estimate the slip distribution corresponding to the basis function
strs_inv_b01.m Perform inversion and output products into folder "results"

Output files (.sd and figures):
Traction in the given direction (Traction_inv_after_drop.sd)
Displacement at the surface at observations location (Displacement_after_drop.sd)
Slip distribution (slip_BF_vectors_after_drop.sd; gmt format: BF_after_drop.sd)
