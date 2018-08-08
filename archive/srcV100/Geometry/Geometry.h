#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include "../Log/Log.h"

// Geometry class
class Geometry {
private:
  float dimentionX_;
  float dimentionY_;
  unsigned int numberOfElementX_;
  unsigned int numberOfElementY_;
public:
  unsigned int numberOfElementsG;
  unsigned int numberOfNodes;
  unsigned int* mesh;
  float* x;
  float* y;
public:
  Geometry(float, float, unsigned int, unsigned int);
  Geometry(const Geometry& geometry);
  void nodeDimentionVector();
  void meshCantilever();
  ~Geometry();
};

#endif
