# Description

This folder contains the necessary files to run the R-C++ verions of ForestGALES/fgr.  
The two C++ files containing the functions from the fgr R package to calculate Critical Wind Speeds for the tree and stand-level versions of ForestGALES. Both files can easily be modified to account for new country-specific allometric or dendrometric equations as required. 
1. The FGdataprep_Norway_cpp.cpp file is a data preparation step and contains functions to check the supplied input data and compute all other variables required for the ForestGALES model. The functions in this file are dependent on the structure of the two ForestGALES parameter files - fgr_constants and species_parameters.txt, found in the folder NorForestGALES/3.CWS_calculations.
2. The FGfunctions_Norway_cpp.cpp file contains the functions to perform the core calculations of the ForestGALES model outputting the critical wind speed for breakage, overturning,
the minimum CWS for damage and the damage type.  

For simplifying the use of both files in the R environment, they were compiled in a bare-minimum R package - fgRcppNorway_1.0

