/* vertex.h
 */

#ifndef VERTEX_H
  #define VERTEX_H
  #include <iostream>
  #include <cstdio>
  #include <list>
  #include <limits>
  #include <vector>

// class vector;

class Vertex {
  float* vertex_height;
  int gap_size_left;
  int gap_size_right;
  int distance_to_gap_left;
  int distance_to_gap_right;

 public:
  Vertex();
  Vertex( float* );

  void SetVertexHeightPtr( float* height ) { vertex_height = height; }
  float GetVertexHeight() { return *vertex_height; }

  void SetGapSizeLeft(int gap_size) { gap_size_left = gap_size; };
  int GetGapSizeLeft() { return gap_size_left; }

  void SetGapSizeRight(int gap_size) { gap_size_right = gap_size; };
  int GetGapSizeRight() { return gap_size_right; }

  void SetDistanceToGapLeft(int dist) { distance_to_gap_left = dist; }
  int GetDistanceToGapLeft() { return distance_to_gap_left; }

  void SetDistanceToGapRight(int dist) { distance_to_gap_right = dist; }
  int GetDistanceToGapRight() { return distance_to_gap_right; }

};

#endif
