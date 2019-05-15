#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <vector>
#include <algorithm>
#include "../Sparse/Sparse.h"
#include "../Log/Log.h"

// Load Struct
struct Load {
  // variables
  std::vector<float> LoadVector_value;
  std::vector<unsigned int> LoadVector_dof_i;
  Sparse* loadSparse;
  // methodes
  Load();
  ~Load();
  void point(int, float, float);
  void build(unsigned int, int);
};
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
  Load* load;
public:
  Geometry();
  ~Geometry();
  void node(float, float);
  void modelBuild();
  void meshQuadrilateral(int,int,int,int);
  void fix(int, int, int);
  Sparse& get_load();
};

#endif
