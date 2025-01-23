Date: 23/01/2025

# Description

In this folder are all the required scripts for calculating the critical wind speed (CWS) for breakage and overturning using the ForestGALES model for each forest pixel uof a forest attribute raster at a national level - here the SR16 data for Norway - using a Linux High Performance Computing cluster.
The scripts can easily be modified and adapted to data from another country.

## Variables
Some variables other than the ones commonly found in forest attribute rasters are required. Namely, the following rasters are required:
- a soil type and rooting depth classification according to the ForestGALES model. A simple translation can be made from national soil type maps to the soil and rooting depth classification from ForestGALES. The output data needs to be extractable for each forest pixel.
- the distance to forest edge and size of the upwind gap for each forest pixel - see the forest-edge-gaps folder.
- (optional: crown snow load for a winter scenario - the script is then appropriately modified to include this extra raster).

## Scripts & routine
The scripts work at the tile level. The Critical Wind Speed for each forest pixel is calculated in parallel with the find_CWS.sh which uses find_CWS_for_tile.sh. The find_CWS_for_tile.sh relies on the R script FGrou_CWS_SR16.R


Step 5: TO DO
Calculate the probability of damage - this may require calculating TOPEX prior using DTM.