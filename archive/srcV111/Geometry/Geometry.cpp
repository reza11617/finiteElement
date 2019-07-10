#include "Geometry.h"

Geometry::Geometry() {
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
  delete load;
  delete dof;
};

void Geometry::node(float tempX, float tempY) {
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
  cudaMallocManaged(&x_d, numberOfNodes*sizeof(float));
  cudaMemcpy(x_d, x, numberOfNodes*sizeof(float), cudaMemcpyHostToDevice);
  cudaMallocManaged(&y_d, numberOfNodes*sizeof(float));
  cudaMemcpy(y_d, y, numberOfNodes*sizeof(float), cudaMemcpyHostToDevice);
  // Dof
  dof->build(numberOfNodes);
  // building the load vector
  load->build(numberOfNodes);
  cudaDeviceSynchronize();
}

Sparse& Geometry::get_load() { return *(load->loadSparse_d);}
float* Geometry::get_x() { return x_d;}
float* Geometry::get_y() { return y_d;}
unsigned int* Geometry::get_mesh() { return mesh_d;}
unsigned int* Geometry::get_fixDofs() { return dof->fixDofs_d;};
unsigned int Geometry::get_fixDofs_size() { return dof->fixDofs_d_size;};
unsigned int* Geometry::get_freeDofs() { return dof->freeDofs_d;};
unsigned int Geometry::get_freeDofs_size() { return dof->freeDofs_d_size;};
unsigned int Geometry::get_numberOfElementsG() {return numberOfElementsG;}
unsigned int Geometry::get_x_y_size() {return numberOfNodes;}
unsigned int Geometry::get_mesh_Size() {return meshTemp.size();}

// Struct Laod
Load::Load() {
  loadSparse_d = new Sparse();
}
Load::~Load() {
  //delete loadSparse_d;
  // loadVec_d -> will be free in the sparse class
  // loadDof_d
}
void Load::point(int NodeNumber, float loadValue_x, float loadValue_y) {
  LoadVector_dof_i.push_back(NodeNumber*2-1);
  LoadVector_value.push_back(loadValue_x);
  LoadVector_dof_i.push_back(NodeNumber*2);
  LoadVector_value.push_back(loadValue_y); 
}
void Load::build(unsigned int nOfNodes) {
  loadVec = &LoadVector_value[0];
  loadDof = &LoadVector_dof_i[0];
  cudaMallocManaged(&loadVec_d, nOfNodes*sizeof(float));
  cudaMemcpy(loadVec_d, loadVec, nOfNodes*sizeof(float), cudaMemcpyHostToDevice);
  cudaMallocManaged(&loadDof_d, nOfNodes*sizeof(unsigned int));
  cudaMemcpy(loadDof_d, loadDof, nOfNodes*sizeof(unsigned int), cudaMemcpyHostToDevice);
  loadSparse_d->set_valueSize(LoadVector_dof_i.size());
  loadSparse_d->set_i(loadDof_d);
  loadSparse_d->set_x(loadVec_d);
  loadSparse_d->set_numberOfRows(2*nOfNodes);
  loadSparse_d->set_numberOfColumns(1);

}

// Struct Dof
void Dof::fix(int nodeNumber, int DOFx, int DOFy) {
  if (DOFx == 1) // 1 means fixed
    dofFixTemp.push_back((nodeNumber+1)*2-1);
  if (DOFy == 1) // 1 means fixed
    dofFixTemp.push_back((nodeNumber+1)*2);  
}

void Dof::build(unsigned int numberOfNodes) {
  fixDofs_d_size = dofFixTemp.size();
  // -- sort fixDofs_d
  std::sort(dofFixTemp.begin(), dofFixTemp.end());
  fixDofs = &dofFixTemp[0];
  // -- copy dof to uniffied memory
  cudaMallocManaged(&fixDofs_d, fixDofs_d_size*sizeof(unsigned int));
  cudaMemcpy(fixDofs_d, fixDofs, fixDofs_d_size*sizeof(float), cudaMemcpyHostToDevice);
  // -- dof free
  freeDofs_d_size = numberOfNodes*2 - fixDofs_d_size;
  cudaMallocManaged(&freeDofs_d, freeDofs_d_size*sizeof(unsigned int));
  int j = 0;
  int c = 0;
  for (unsigned int i = 1; i <= numberOfNodes*2; i++) {
    if (i != fixDofs_d[j]) {
    freeDofs_d[c++] = i;
    } else {
      j++;
    }
  }
}

Dof::Dof() {
};

Dof::~Dof() {
  cudaFree(fixDofs_d);
  cudaFree(freeDofs_d);
}
