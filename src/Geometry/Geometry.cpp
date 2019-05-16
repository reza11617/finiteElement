#include "Geometry.h"

Geometry::Geometry() {
  Log::Logger().Info("Geometry Created");
  numberOfNodes = 0;
  numberOfElementsG = 0;
  load = new Load();
  dof = new Dof();
};

Geometry::~Geometry() {
  Log::Logger().Info("Geometry Deleted");
  cudaFree(x_d);
  cudaFree(y_d);
  cudaFree(mesh_d);
  cudaFree(thickness_array_d);
  
  delete load;
  delete dof;
};

void Geometry::node(double tempX, double tempY) {
  xDim.push_back(tempX);
  yDim.push_back(tempY);
};

void Geometry::meshQuadrilateral(int node1, int node2, int node3, int node4) {
  numberOfElementsG++;
  meshTemp.push_back(node1);
  meshTemp.push_back(node2);
  meshTemp.push_back(node3);
  meshTemp.push_back(node4);
}

void Geometry::modelBuild() {
  // copy mesh to unified memory
  mesh = &meshTemp[0];
  cudaMallocManaged(&mesh_d, meshTemp.size()*sizeof(unsigned int));
  cudaMemcpy(mesh_d, mesh, meshTemp.size()*sizeof(unsigned int), cudaMemcpyHostToDevice);
  // copy x and y to unified memory
  numberOfNodes = xDim.size();
  x = &xDim[0]; // return them as array
  y = &yDim[0];
  cudaMallocManaged(&x_d, numberOfNodes*sizeof(double));
  cudaMemcpy(x_d, x, numberOfNodes*sizeof(double), cudaMemcpyHostToDevice);
  cudaMallocManaged(&y_d, numberOfNodes*sizeof(double));
  cudaMemcpy(y_d, y, numberOfNodes*sizeof(double), cudaMemcpyHostToDevice);
  // Dof
  dof->build(numberOfNodes);
  // building the load vector
  load->build(dof->get_free(), dof->get_freeSize());

  // creating the thickness array
  cudaMallocManaged(&thickness_array_d, meshTemp.size()*sizeof(double));
  for(int iter = 0; iter < meshTemp.size(); iter++) {thickness_array_d[iter] = thickness;}

  // cuda device sync
  cudaDeviceSynchronize();
}

void Geometry::set_thickness(double t) {  thickness = t;}
const double* Geometry::get_thickness() const {  return thickness_array_d;}
double* Geometry::get_x() { return x_d;}
double* Geometry::get_y() { return y_d;}
unsigned int* Geometry::get_mesh() { return mesh_d;}
unsigned int Geometry::get_numberOfElementsG() {return numberOfElementsG;}
unsigned int Geometry::get_x_y_size() {return numberOfNodes;}
unsigned int Geometry::get_mesh_Size() {return meshTemp.size();}

const Dof&  Geometry::get_Dof()  const { return *dof; };
const Load& Geometry::get_Load() const { return *load;};
// Struct Laod
Load::Load() {
}
Load::~Load() {
  cudaFree(loadVec);
}
void Load::point(int NodeNumber, double loadValue_x, double loadValue_y) {
  LoadVector_dof_i.push_back(NodeNumber*2-1);
  LoadVector_value.push_back(loadValue_x);
  LoadVector_dof_i.push_back(NodeNumber*2);
  LoadVector_value.push_back(loadValue_y); 
}
void Load::build(unsigned int* dofFree, unsigned int freeSize) {
  cudaMallocManaged(&loadVec, freeSize*sizeof(double));
  cudaMemset(loadVec, 0, freeSize*sizeof(double));
  for (unsigned int i = 0; i < LoadVector_dof_i.size(); i++) {
    loadVec[dofFree[LoadVector_dof_i[i]-1]-1] = LoadVector_value[i];
  }
}

double* Load::get_vector() const { return loadVec;}

// Struct Dof
Dof::Dof() { };

Dof::~Dof() {
  cudaFree(free);
}

void Dof::fix(int nodeNumber, int DOFx, int DOFy) {
  if (DOFx == 1) // 1 means fixed
    dofFixTemp.push_back((nodeNumber+1)*2-1);
  if (DOFy == 1) // 1 means fixed
    dofFixTemp.push_back((nodeNumber+1)*2);  
}

void Dof::build (unsigned int nNode) {
  cudaMallocManaged(&free, 2*nNode*sizeof(unsigned int));
  cudaMemset(free, 1,  2*nNode*sizeof(unsigned int));
  fixSize = dofFixTemp.size();
  freeSize = nNode*2-fixSize;
  for (unsigned int i = 0; i < fixSize; i++) {
    free[dofFixTemp[i]-1] = 0;
  }
  unsigned int counter = 1;
  for (unsigned int i = 0; i < 2*nNode; i++) {
    if (free[i]) {free[i] = counter++;} 
  }
}

unsigned int* Dof::get_free() const { return free;} 
unsigned int Dof::get_freeSize() const {return freeSize; }
unsigned int Dof::get_fixSize () const {return fixSize;  }



