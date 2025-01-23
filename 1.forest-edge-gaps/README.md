Date: 23/01/2025

# Description
This code calculates the distance to the forest edge and the size of the upwind gap for each individual forest pixel in a forest raster, in the four cardinal directions. A suite of requirements are chosen and applied in the code to discard forest pixels with insufficient number of trees, small tree height, and to limit the distance over which the forest edge is searched for, and to discard forest pixels. Note that some parameters should be adjusted to your own settings.

The gap size and distance to gap is calculated in parallel with the find_gaps_and_distances.sh which uses find_gaps_and_distances_for_tile.sh relying on the c++ program that currently is named find_gaps_and_distances (compiled with the find_gaps_and_distances.cpp, vertex_grid.cpp, vertex_grid.h, vertex.h, vertex.cpp files).
The files should be copied to project_dir.

## compile cpp program

```shell
cd cpp
project_dir="[your_project_dir]"
scp find_gaps_and_distances.cpp \
    vertex_grid.cpp  vertex.h \
    vertex.cpp \
    vertex_grid.h \
    ${project_dir}/cpp/
```

compile cpp (cd cpp, on server)

```shell
g++ -O3 find_gaps_and_distances.cpp vertex.cpp vertex_grid.cpp \
    -g -Wall -I/usr/include/gdal -L/usr/lib -lgdal -lm -o \
    ../find_gaps_and_distances
```

## run the function

```shell
project_dir="[your_project_dir]"
cd ${project_dir}
mkdir log
{ time ./find_gaps_and_distances.sh ; } > log/YYYMMDD_find_gaps_and_distances.log 2>&1
```
