#include "Geometry.h"

Geometry::Geometry(float dim_x, float dim_y, unsigned int numberOfElementX, unsigned int numberOfElementy)
  : dimentionX_(dim_x), dimentionY_(dim_y), numberOfElementX_(numberOfElementX), numberOfElementY_(numberOfElementy)
{
  nodeDimentionVector();
  meshCantilever();
};

Geometry::Geometry(const Geometry& geometry)
  : dimentionX_(geometry.dimentionX_), dimentionY_(geometry.dimentionY_), numberOfElementX_(geometry.numberOfElementX_), numberOfElementY_(geometry.numberOfElementY_)
{
  Log::Logger().Warning("Copied");
};

Geometry::~Geometry()
{
  Log::Logger().Info("Geometry Deleted");
  delete[] mesh;
  delete[] x;
  delete[] y;
};


void Geometry::nodeDimentionVector() {
  // this function numbers the nodes the use method for numbering is
  /*
    6 --- 7 --- 8
    3 --- 4 --- 5
    0 --- 1 --- 2
  */
  numberOfNodes = (numberOfElementX_ + 1)*(numberOfElementY_ + 1);
  x = new float[numberOfNodes];
  y = new float[numberOfNodes];
  float incr_x = dimentionX_/numberOfElementX_;
  float incr_y = dimentionY_/numberOfElementY_;
  for (unsigned int j = 0; j <= numberOfElementY_; j++) {
    for (unsigned int i = 0; i <= numberOfElementX_; i++) {
      x[i + j*(numberOfElementX_+1)] = i * incr_x;
      y[i + j*(numberOfElementX_+1)] = j * incr_y;
    }
  }
};

void Geometry::meshCantilever()
{
  // this function meshes a contilever
  /*
    6 --- 7 --- 8
    |  3  |  4  |
    3 --- 4 --- 5
    |  1  |  2  |
    0 --- 1 --- 2
  */
  numberOfElementsG = numberOfElementX_*numberOfElementY_;
  unsigned int elementNodes = 4; // how many nodes per element
  unsigned int eleNumber = 0;
  mesh = new unsigned int[numberOfElementsG*elementNodes];
  unsigned int xElement = numberOfElementX_+1;
  for (unsigned int i = 0; i < numberOfElementY_; i++){
    for (unsigned int j = 0; j < numberOfElementX_; j++){
      unsigned int eleNumber = (j + i*numberOfElementY_)*elementNodes;
      mesh[eleNumber++] = i*xElement+j;
      mesh[eleNumber++] = i*xElement+1+j;
      mesh[eleNumber++] = (i+1)*xElement+1+j;
      mesh[eleNumber++] = (i+1)*xElement+j;
    }
  }
};

