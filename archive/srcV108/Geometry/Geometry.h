#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <vector>
#include <algorithm>
#include "../Log/Log.h"

// Geometry class
class Geometry {
private:
  unsigned int counter;
  std::vector<unsigned int> dofFixTemp;
  std::vector<unsigned int> meshTemp;
  std::vector<float> xDim;
  std::vector<float> yDim;
public:
  unsigned int numberOfNodes;
  unsigned int numberOfElementsG;
  unsigned int dofFixSize;
  unsigned int* dofFix;
  unsigned int dofFreeSize;
  unsigned int* dofFree;
  float* x;
  float* y;
  unsigned int* mesh;
public:
  Geometry();
  ~Geometry();
  void node(float, float);
  void modelBuild();
  void meshQuadrilateral(int,int,int,int);
  void fix(int, int, int);
};

#endif
