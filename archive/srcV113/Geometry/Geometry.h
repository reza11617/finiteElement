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
  float* loadVec; 
  // methodes
  Load();
  ~Load();
  void point(int, float, float);
  void build(unsigned int*, unsigned int);
  float* get_vector() const;
};

// dof Struct
struct Dof {
  // variables
  unsigned int* free; // if zero means it is fixed if numbered is the new dof value
  unsigned int freeSize;
  unsigned int fixSize;
  std::vector<unsigned int> dofFixTemp;
  // methods
  Dof();
  ~Dof();
  void fix(int, int, int);
  void build(unsigned int);
  unsigned int* get_free() const;
  unsigned int get_freeSize() const;
  unsigned int get_fixSize () const;
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
  float* get_x();
  float* get_y();
  unsigned int get_x_y_size();
  unsigned int get_numberOfElementsG();
  unsigned int get_mesh_Size();
  unsigned int* get_mesh();
  const Dof& get_Dof() const;
  const Load& get_Load() const;
};

#endif
