Date: 23/01/2025

# NorForestGALES
Adapting ForestGALES for its application across country-level high resolution forest rasters.

This repository is part of the Climate Smart Forestry Norway project funded by the Norwegian Research Council (NFR project 302701, 2020-2024). The goal is to calculate the critical wind speed (CWS) for breakage and overturning for each individual pixel in a nationwide forest raster using the ForestGales model described in the R package 'fgr' (Locatelli & al, 2022). The code is associated with the scientific publication accepted for publication in Forest Ecosystems: Merlin, M., Locatelli, T., Gardiner, B., Astrup, R. 2025. Large-scale modelling wind damage vulnerability through combination of high-resolution forest resources maps and ForestGALES. The example here is with the Norwegian forest attribute raster collection SR16 (https://kilden.nibio.no). The parameters and some equations of the model were modified to better represent the Norwegian forests. 

## Contributors
The main work was done by Morgane Merlin, with contributions from Barry Gardiner (IEFC) and Tom Locatelli (UK Forest Research) for the ForestGales model and fgr functions, Nicolas Cattaneo (NIBIO), Nils Egil Søvde (NIBIO).

## Description
This repository contains the different scripts used in the routine implemented for calculating CWS for SR16 pixels under a summer scenario on relatively flat terrain (i.e. no snow loading of tree crown or effect of slope are considered here).
The repository contains:
1. a folder with the scripts and software for calculating the distance to the forest edge and size of the upwind forest gap in four cardinal directions. The code was written by Nils Egil Søvde (NIBIO). Currently, the edge of the forest is defined as forest to non-forest transitions, and the clearcuts. Clear limits to what is considered forest are used, based on the NFI common thresholds. Namely, mean forest height > 10 m and number of trees > 400.
3. a folder with the scripts invoking the fgr R package functions translated into C++ (not available for sharing on Github) for parallel processing in a Linux High Performance Computing cluster

## Routine
The find_gaps_and_distances.sh script from the forest-edge-gaps can be run first, then the find_CWS_Norway.sh script from the CWS_Calculation can be run after. Or the two can be combined in a single script which calls on them sequentially.

## License
GPLv3 following the licensing of ForestGALES/fgr.
***

