#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <vector>
#include "../Log/Log.h"

// Geometry class
class Geometry {
private:
  unsigned int counter;
public:
  unsigned int numberOfNodes;
  unsigned int numberOfElementsG;
  std::vector<unsigned int> meshTemp;
  std::vector<float> xDim;
  std::vector<float> yDim;
  float* x;
  float* y;
  unsigned int* mesh;
public:
  Geometry();
  void node(float, float);
  void modelBuild();
  void meshQuadrilateral(int,int,int,int);
  ~Geometry();
};

#endif
