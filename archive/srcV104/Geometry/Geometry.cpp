#include "Geometry.h"

Geometry::Geometry() {
  numberOfNodes = 0;
  numberOfElementsG = 0;
};

Geometry::~Geometry() {
  Log::Logger().Info("Geometry Deleted");
  delete[] x;
  delete[] y;
  delete[] mesh;
};

void Geometry::node(float tempX, float tempY) {
  xDim.push_back(tempX);
  yDim.push_back(tempY);
};

void Geometry::modelBuild() {
  numberOfNodes = xDim.size();
  mesh = new unsigned int[meshTemp.size()];
  mesh = &meshTemp[0];
  x = new float[numberOfNodes];
  y = new float[numberOfNodes];
  x = &xDim[0]; // return them as array
  y = &yDim[0];
}

void Geometry::meshQuadrilateral(int node1, int node2, int node3, int node4) {
  numberOfElementsG++;
  meshTemp.push_back(node1);
  meshTemp.push_back(node2);
  meshTemp.push_back(node3);
  meshTemp.push_back(node4);
}

