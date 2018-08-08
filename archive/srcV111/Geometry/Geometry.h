#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cuda.h>
#include <cuda_runtime.h>
#include "../Sparse/Sparse.h"
#include "../Log/Log.h"

// Load Struct
struct Load {
  // variables
  std::vector<float> LoadVector_value;
  std::vector<unsigned int> LoadVector_dof_i;
  float* loadVec; float* loadVec_d;
  unsigned int* loadDof; unsigned int* loadDof_d;
  Sparse* loadSparse_d;
  // methodes
  Load();
  ~Load();
  void point(int, float, float);
  void build(unsigned int);
};

// dof Struct
struct Dof {
  // variables
  std::vector<unsigned int> dofFixTemp;
  unsigned int* fixDofs;
  unsigned int* fixDofs_d;
  unsigned int* freeDofs_d;
  unsigned int fixDofs_d_size;
  unsigned int freeDofs_d_size;
  // methods
  Dof();
  ~Dof();
  void fix(int, int, int);
  void build(unsigned int);
};

// Geometry class
class Geometry {
private:
  unsigned int counter;
  std::vector<unsigned int> meshTemp;
  std::vector<float> xDim;
  std::vector<float> yDim;
  unsigned int numberOfNodes;
  unsigned int numberOfElementsG;
  float* x; float* x_d;
  float* y; float* y_d;
  unsigned int* mesh;unsigned int* mesh_d;
public:
  Load* load;
  Dof* dof;
public:
  Geometry();
  ~Geometry();
  void node(float, float);
  void modelBuild();
  void meshQuadrilateral(int,int,int,int);
  Sparse& get_load();
  float* get_x();
  float* get_y();
  unsigned int get_x_y_size();
  unsigned int get_numberOfElementsG();
  unsigned int get_mesh_Size();
  unsigned int* get_mesh();
  unsigned int* get_freeDofs();
  unsigned int get_freeDofs_size();
  unsigned int* get_fixDofs();
  unsigned int get_fixDofs_size();
};

#endif
