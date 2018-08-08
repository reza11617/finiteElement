#include "Geometry.h"

Geometry::Geometry() {
  numberOfNodes = 0;
  numberOfElementsG = 0;
  load = new Load();
};

Geometry::~Geometry() {
  Log::Logger().Info("Geometry Deleted");
  delete[] x;
  delete[] y;
  delete[] mesh;
  delete[] dofFree;
  delete[] dofFix;
  delete load;
};

void Geometry::node(float tempX, float tempY) {
  xDim.push_back(tempX);
  yDim.push_back(tempY);
};

void Geometry::modelBuild() {
  numberOfNodes = xDim.size();
  //mesh = new unsigned int[meshTemp.size()];
  mesh = &meshTemp[0];
  //x = new float[numberOfNodes];
  //y = new float[numberOfNodes];
  x = &xDim[0]; // return them as array
  y = &yDim[0];
  //dofFix
  dofFixSize = dofFixTemp.size();
  // -- sort dofFix
  std::sort(dofFixTemp.begin(), dofFixTemp.end());
  //dofFix = new unsigned int[dofFixTemp.size()];
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
  // building the load vector
  load->build(numberOfNodes,1);
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

Sparse& Geometry::get_load() {
  return *load->loadSparse;
}

// Struct Laod
Load::Load() {
  loadSparse = new Sparse();
}
Load::~Load() {
  delete loadSparse;
}
void Load::point(int NodeNumber, float loadValue_x, float loadValue_y) {
  LoadVector_dof_i.push_back(NodeNumber*2-1);
  LoadVector_value.push_back(loadValue_x);
  LoadVector_dof_i.push_back(NodeNumber*2);
  LoadVector_value.push_back(loadValue_y); 
}
void Load::build(unsigned int nOfNodes, int deviceID) {
  loadSparse->set_valueSize(LoadVector_dof_i.size());
  loadSparse->i = &LoadVector_dof_i[0];
  loadSparse->value = &LoadVector_value[0];
  loadSparse->set_numberOfRows(2*nOfNodes);
  loadSparse->set_numberOfColumns(1);
  loadSparse->set_device(deviceID);
}

