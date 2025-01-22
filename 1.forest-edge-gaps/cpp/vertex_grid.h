/* vertex_grid.h
 */

#ifndef VERTEX_GRID_H
  #define VERTEX_GRID_H
  #include <iostream>
  #include <fstream>
  #include <string>
  #include <stdlib.h>
  #include <cstring>
  #include <vector>
  #include <limits>

  #include "vertex.h"

class Vertex;

class Vertex_Grid {
  const unsigned int ns_size;
  const unsigned int ew_size;
  const int projection_epsg;
  const double ew_start;
  const double ew_res;
  const double ew_skew;
  const double ns_start;
  const double ns_skew;
  const double ns_res;
  std::vector<std::vector<Vertex> > vertex_array;


public:
  Vertex_Grid( const unsigned int, const unsigned int, const int,
           double, double, double, double, double, double );
  void CalcEW(float, float, int);
  void CalcNS(float, float, int);
  Vertex* GetVertex( int, int );

};

#endif
