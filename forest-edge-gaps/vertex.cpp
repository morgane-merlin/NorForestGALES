/* vertex.cpp
 */

#include "vertex.h"

Vertex::Vertex() {
  vertex_height = NULL;
  gap_size_left = std::numeric_limits<int>::max();
  gap_size_right = std::numeric_limits<int>::max();
  distance_to_gap_left = std::numeric_limits<int>::max();
  distance_to_gap_right = std::numeric_limits<int>::max();
}

Vertex::Vertex( float* height ) {
  vertex_height = height;
  gap_size_left = std::numeric_limits<int>::max();
  gap_size_right = std::numeric_limits<int>::max();
  distance_to_gap_left = std::numeric_limits<int>::max();
  distance_to_gap_right = std::numeric_limits<int>::max();
}
