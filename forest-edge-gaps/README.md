This code calculates the distance to the forest edge and the size of the upwind gap for each individual forest pixel in a forest raster, in the four cardinal directions. A suite of requirements are chosen and applied in the code to discard forest pixels with insufficient number of trees, small tree height, and to limit the distance over which the forest edge is searched for, and to discard forest pixels 

The gap size and distance to gap is calculated in parallel.

find_gaps_and_distances.sh uses find_gaps_and_distances_for_rute.sh,
which uses the c++ program that currently is named
find_gaps_and_distances.

The files should be copied to project_dir.

Note that a lot of stuff in find_gaps_and_distances.sh should be
changed to your own settings.

## compile cpp program

```shell
cd cpp
project_dir="PROJECT"
scp find_gaps_and_distances.cpp \
    vertex_grid.cpp  vertex.h \
    vertex.cpp \
    vertex_grid.h \
    nis@vroom5.int.nibio.no:${project_dir}/cpp/
```

compile cpp (cd cpp, on server)

```shell
g++ -O3 find_gaps_and_distances.cpp vertex.cpp vertex_grid.cpp \
    -g -Wall -I/usr/include/gdal -L/usr/lib -lgdal -lm -o \
    ../find_gaps_and_distances
```

## copy bash scripts

```shell
cd bash
project_dir="PROJECT"
scp find_gaps_and_distances.sh find_gaps_and_distances_for_rute.sh \
    nis@vroom5.int.nibio.no:${project_dir}/
```

## run the function

```shell
project_dir="PROJECT"
cd ${project_dir}
mkdir log
{ time ./find_gaps_and_distances.sh ; } > log/YYYMMDD_find_gaps_and_distances.log 2>&1
```
