/* vertex_grid.cpp
 */

#include <math.h>
#include <algorithm>

#include "vertex_grid.h"

Vertex_Grid::Vertex_Grid(const unsigned int ns_size_input,
                         const unsigned int ew_size_input,
                         const int vertex_projection_epsg,
                         double vertex_ew_start, double vertex_ew_res,
                         double vertex_ew_skew, double vertex_ns_start,
                         double vertex_ns_skew, double vertex_ns_res)
  : ns_size(ns_size_input), ew_size(ew_size_input),
    projection_epsg(vertex_projection_epsg),
    ew_start(vertex_ew_start), ew_res(vertex_ew_res), ew_skew(vertex_ew_skew),
    ns_start(vertex_ns_start), ns_skew(vertex_ns_skew), ns_res(vertex_ns_res) {

  // epsilon = std::numeric_limits<double>::epsilon();
  // // printf("epsilon: %e", epsilon);

  vertex_array.resize(ns_size);
  for ( unsigned int i = 0; i < ns_size; i++ ) {
    vertex_array[i].resize(ew_size);
  }
}

void Vertex_Grid::CalcEW(float height_nodata_value, float gap_height,
                         int max_distance) {

  int resolution = 16;

  for ( unsigned int i = 0; i < ns_size; i++ ) {

    // first we go from west to east...

    int gap_size = std::numeric_limits<int>::max();
    int distance_to_gap = std::numeric_limits<int>::max();

    for ( unsigned int j = 0; j < ew_size-1; j++ ) {

      Vertex* current_vertex = GetVertex(i, j);
      Vertex* next_vertex = GetVertex(i, j+1);

      float current_height = current_vertex->GetVertexHeight();
      float next_height = next_vertex->GetVertexHeight();

      if (j == 0) {
        if (current_height != height_nodata_value) { // forest
          gap_size = 0;
          distance_to_gap = max_distance;
        } else { // not forest
          gap_size = max_distance;
          distance_to_gap = 0;
        }
      }

      // NOTE: The gap_size and distance_to_gap is set here, from the
      //   previous calculation.

      if (gap_size > max_distance) { gap_size = max_distance; }
      current_vertex->SetGapSizeLeft(gap_size);

      if (distance_to_gap > max_distance) { distance_to_gap = max_distance; }
      current_vertex->SetDistanceToGapLeft(distance_to_gap);
      
      if ((current_height == height_nodata_value) // not forest
          || (current_height < gap_height) // clearcut
          )  {
        if ((next_height == height_nodata_value) // notforest
            || (next_height < gap_height) // clearcut
            ) {
          gap_size += resolution;
          // the distance to gap is not changed (zero)
        } else { // from gap to non clearcut forest
          distance_to_gap = 0;
        }
      } else { // from non clearcut forest
        if ((next_height == height_nodata_value) // notforest
            || (next_height < gap_height) // clearcut
            ) {
          gap_size = resolution;
          distance_to_gap = 0;
        } else { // to non clearcut forest
          distance_to_gap += resolution;
        }
      }
    }


    // now we go from east to west...

    gap_size = std::numeric_limits<int>::max();
    distance_to_gap = std::numeric_limits<int>::max();

    for ( unsigned int j = ew_size-1; j > 0; j-- ) {

      Vertex* current_vertex = GetVertex(i, j);
      Vertex* next_vertex = GetVertex(i, j-1);

      float current_height = current_vertex->GetVertexHeight();
      float next_height = next_vertex->GetVertexHeight();

      if (j == (ew_size-1)) {
        if (current_height != height_nodata_value) { // forest
          gap_size = 0;
          distance_to_gap = max_distance;
        } else { // not forest
          gap_size = max_distance;
          distance_to_gap = 0;
        }
      }

      // NOTE: The gap_size and distance_to_gap is set here, from the
      //   previous calculation.

      if (gap_size > max_distance) { gap_size = max_distance; }
      current_vertex->SetGapSizeRight(gap_size);

      if (distance_to_gap > max_distance) { distance_to_gap = max_distance; }
      current_vertex->SetDistanceToGapRight(distance_to_gap);
      
      if ((current_height == height_nodata_value) // not forest
          || (current_height < gap_height) // clearcut
          )  {
        if ((next_height == height_nodata_value) // notforest
            || (next_height < gap_height) // clearcut
            ) {
          gap_size += resolution;
          // the distance to gap is not changed (zero)
        } else { // from gap to non clearcut forest
          distance_to_gap = 0;
        }
      } else { // from non clearcut forest
        if ((next_height == height_nodata_value) // notforest
            || (next_height < gap_height) // clearcut
            ) {
          gap_size = resolution;
          distance_to_gap = 0;
        } else { // to non clearcut forest
          distance_to_gap += resolution;
        }
      }
    }
  }
}



void Vertex_Grid::CalcNS(float height_nodata_value, float gap_height,
                         int max_distance) {

  int resolution = 16;

  for ( unsigned int j = 0; j < ew_size; j++ ) {

    // first we go from north to south...

    int gap_size = std::numeric_limits<int>::max();
    int distance_to_gap = std::numeric_limits<int>::max();

    for ( unsigned int i = 0; i < ns_size-1; i++ ) {

      Vertex* current_vertex = GetVertex(i, j);
      Vertex* next_vertex = GetVertex(i+1, j);

      float current_height = current_vertex->GetVertexHeight();
      float next_height = next_vertex->GetVertexHeight();

      if (i == 0) {
        if (current_height != height_nodata_value) { // forest
          gap_size = 0;
          distance_to_gap = max_distance;
        } else { // not forest
          gap_size = max_distance;
          distance_to_gap = 0;
        }
      }

      // NOTE: The gap_size and distance_to_gap is set here, from the
      //   previous calculation.

      if (gap_size > max_distance) { gap_size = max_distance; }
      current_vertex->SetGapSizeLeft(gap_size);

      if (distance_to_gap > max_distance) { distance_to_gap = max_distance; }
      current_vertex->SetDistanceToGapLeft(distance_to_gap);
      
      if ((current_height == height_nodata_value) // not forest
          || (current_height < gap_height) // clearcut
          )  {
        if ((next_height == height_nodata_value) // notforest
            || (next_height < gap_height) // clearcut
            ) {
          gap_size += resolution;
          // the distance to gap is not changed (zero)
        } else { // from gap to non clearcut forest
          distance_to_gap = 0;
        }
      } else { // from non clearcut forest
        if ((next_height == height_nodata_value) // notforest
            || (next_height < gap_height) // clearcut
            ) {
          gap_size = resolution;
          distance_to_gap = 0;
        } else { // to non clearcut forest
          distance_to_gap += resolution;
        }
      }
    }


    // now we go from south to north...

    gap_size = std::numeric_limits<int>::max();
    distance_to_gap = std::numeric_limits<int>::max();

    for ( unsigned int i = ns_size-1; i > 0; i-- ) {

      Vertex* current_vertex = GetVertex(i, j);
      Vertex* next_vertex = GetVertex(i-1, j);

      float current_height = current_vertex->GetVertexHeight();
      float next_height = next_vertex->GetVertexHeight();

      if (i == (ns_size-1)) {
        if (current_height != height_nodata_value) { // forest
          gap_size = 0;
          distance_to_gap = max_distance;
        } else { // not forest
          gap_size = max_distance;
          distance_to_gap = 0;
        }
      }

      // NOTE: The gap_size and distance_to_gap is set here, from the
      //   previous calculation.

      if (gap_size > max_distance) { gap_size = max_distance; }
      current_vertex->SetGapSizeRight(gap_size);

      if (distance_to_gap > max_distance) { distance_to_gap = max_distance; }
      current_vertex->SetDistanceToGapRight(distance_to_gap);
      
      if ((current_height == height_nodata_value) // not forest
          || (current_height < gap_height) // clearcut
          )  {
        if ((next_height == height_nodata_value) // notforest
            || (next_height < gap_height) // clearcut
            ) {
          gap_size += resolution;
          // the distance to gap is not changed (zero)
        } else { // from gap to non clearcut forest
          distance_to_gap = 0;
        }
      } else { // from non clearcut forest
        if ((next_height == height_nodata_value) // notforest
            || (next_height < gap_height) // clearcut
            ) {
          gap_size = resolution;
          distance_to_gap = 0;
        } else { // to non clearcut forest
          distance_to_gap += resolution;
        }
      }
    }
  }
}


Vertex* Vertex_Grid::GetVertex( int i, int j ) {
  Vertex* vert;
  try {
    vert = &vertex_array.at(i).at(j);
  } catch (...) { throw; }
  
  return vert;
}
