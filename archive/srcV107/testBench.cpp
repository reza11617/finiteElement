#include "testBench.h"

int main() {
  // start logging
  Log::Logger().setLevel(Log::LevelInfo);
  // GEOMETRY (building the contiliver)
  Geometry* p_cantilever;
  p_cantilever = geometry();
  // Build material vector (E, nu)
  Material mat(1000.0f, 0.3f);
  // Building stiffness matrix on CPU single core, CPU multiple core & GPU.
  StiffnessMatrixFirstOrder* stiffMat       = new StiffnessMatrixSingleCPU(mat,*p_cantilever,4);
  StiffnessMatrixFirstOrder* stiffMatParCPU = new StiffnessMatrixParallelCPU(mat,*p_cantilever,4,5);
  StiffnessMatrixFirstOrder* stiffMatGPU    = new StiffnessMatrixGPU(mat,*p_cantilever,4);
  Sparse* k = stiffMat->GetStiffnessMatrix();
  Sparse* kCPU = stiffMatParCPU->GetStiffnessMatrix();
  Sparse* kGPU = stiffMatGPU->GetStiffnessMatrix();
  Sparse::Printer(*k);
  std::cout << Sparse::Compare(*k,*kGPU) << " is the Max Error"<<std::endl;  
  Log::Logger().Info(p_cantilever->dofFix ,p_cantilever->dofFixSize);
  Log::Logger().Info(p_cantilever->dofFree,p_cantilever->dofFreeSize);
  k->Assemble(p_cantilever->dofFree,p_cantilever->dofFreeSize);
  Sparse::Printer(*k);
  delete stiffMat;
  delete stiffMatParCPU;
  delete stiffMatGPU;
  delete p_cantilever;
  return 0;
}

Geometry* geometry() {
  // GEOMETRY (building the contiliver)
  float dimentionX = 200.0; int numberOfElementX = 1; // dimention of x and number of element in x
  float dimentionY = 100.0; int numberOfElementY = 1; 
  float incr_x = dimentionX/numberOfElementX; // increment between each node 
  float incr_y = dimentionY/numberOfElementY;
  Geometry* cantilever = new Geometry();
  /*
  6 --- 7 --- 8
  3 --- 4 --- 5
  0 --- 1 --- 2
  */
  for (unsigned int j = 0; j <= numberOfElementY; j++) {
    for (unsigned int i = 0; i <= numberOfElementX; i++)
      cantilever->node(i * incr_x, j * incr_y);
      // Node(x,y): build nodes for the geometry
  }
  for (unsigned int j = 0; j <= numberOfElementY; j++)
    cantilever->fix(j*(numberOfElementX+1),1,1);
  // Mesh the Cantilever beam
    /*
    6 --- 7 --- 8
    |  3  |  4  |
    3 --- 4 --- 5
    |  1  |  2  |
    0 --- 1 --- 2
  */
  int xElement = numberOfElementX+1;
  for (unsigned int i = 0; i < numberOfElementY; i++){
    for (unsigned int j = 0; j < numberOfElementX; j++){
      cantilever->meshQuadrilateral(i*xElement+j,i*xElement+1+j,(i+1)*xElement+1+j,(i+1)*xElement+j); 
      // MeshQuadrilateral({Node1,Node2,Node3,Node4})
    }
  }
  cantilever->modelBuild();
  return cantilever;
}



