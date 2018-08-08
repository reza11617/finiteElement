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
  delete[] dofFree;
  delete[] dofFix;
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
  //dofFix
  dofFixSize = dofFixTemp.size();
  // -- sort dofFix
  std::sort(dofFixTemp.begin(), dofFixTemp.end());
  dofFix = new unsigned int[dofFixTemp.size()];
  dofFix = &dofFixTemp[0];
  // dofFree
  dofFreeSize = numberOfNodes*2 - dofFixSize;
  dofFree = new unsigned int[dofFreeSize];
  int j = 0;
  int c = 0;
  for (unsigned int i = 1; i <= numberOfNodes*2; i++) {
    if (i != dofFix[j]) {
    dofFree[c++] = i;
    } else {
      j++;
    }
  }
}

void Geometry::meshQuadrilateral(int node1, int node2, int node3, int node4) {
  numberOfElementsG++;
  meshTemp.push_back(node1);
  meshTemp.push_back(node2);
  meshTemp.push_back(node3);
  meshTemp.push_back(node4);
}

void Geometry::fix(int nodeNumber, int DOFx, int DOFy) {
  if (DOFx == 1) // 1 means fixed
    dofFixTemp.push_back((nodeNumber+1)*2-1);
  if (DOFy == 1) // 1 means fixed
    dofFixTemp.push_back((nodeNumber+1)*2);  
}
