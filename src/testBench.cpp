#include "testBench.h"

int main() {
  // -- start logging
  Log::Logger().setLevel(Log::LevelInfo);
  // -- GEOMETRY (building the contiliver)
  Geometry p_cantilever = geometry();
  // -- Build material vector (E, nu)
  Material mat(1000.0d, 0.3d);
  // -- Building stiffness matrix on CPU single core, CPU multiple core & GPU.
  StiffnessMatrixFirstOrder* stiffMat       = new StiffnessMatrixSingleCPU(mat,p_cantilever,2);
  StiffnessMatrixFirstOrder* stiffMatParCPU = new StiffnessMatrixParallelCPU(mat,p_cantilever,2,10);
  StiffnessMatrixFirstOrder* stiffMatGPU    = new StiffnessMatrixGPU(mat,p_cantilever,2);
  StiffnessMatrixFirstOrder* stiffMatGPU2    = new StiffnessMatrixGPU(mat,p_cantilever,2);
  // -- Calculate the stiffness matrix and do the assembly
  Sparse &k1 = stiffMat->GetStiffnessMatrix();
  Sparse &k2 = stiffMatParCPU->GetStiffnessMatrix();
  Sparse &k = stiffMatGPU->GetStiffnessMatrix();
  Sparse &k3 = stiffMatGPU2->GetStiffnessMatrix();
  Recorder::File().SparseMatrix("kSingleCPU.out", k);
  k.STLAssemble(0.001);
  k1.ThrustAssemble(0.001);
  k2.STLAssemble2(0.001);
  k3.ThrustAssemble2(0.001);
  // -- Solver
  //double* displacement;
  //cudaMallocManaged(&displacement, k.get_numberOfRows()*sizeof(double));
  //SolverSp(k.get_value(), k.get_i(), k.get_j(), k.get_valueSize(), k.get_numberOfRows(),
  //p_cantilever.get_Load().get_vector(), displacement);
  // -- Recorder
  //Log::Logger().Info(displacement[k.get_numberOfRows()-1]);
  //Recorder::File().matrix("displacement.out", displacement, k.get_numberOfRows());
  // -- Delete all variables in Heap memory
  //cudaFree(displacement);
  delete stiffMat;
  delete stiffMatParCPU;
  delete stiffMatGPU;
  

  //Sparse kCPU = stiffMatParCPU->GetStiffnessMatrix();
  //
  
  //std::cout << Sparse::Compare(k,kGPU) << " is the Max Error"<<std::endl;
  
  // Log::Logger().Info(f.value[f.valueSize-1]);

  return 0;
}

Geometry& geometry() {
  // GEOMETRY (building the contiliver)
  double dimentionX = 10.0; int numberOfElementX = 2000; // dimention of x and number of element in x
  double dimentionY =  1.0; int numberOfElementY = 200; 
  double incr_x = dimentionX/numberOfElementX; // increment between each node 
  double incr_y = dimentionY/numberOfElementY;
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
    cantilever->dof->fix(j*(numberOfElementX+1),1,1);
  cantilever->load->point((numberOfElementX+1)*(numberOfElementY+1),0.0d,-0.1d); // pointLoad(NodeNumber,DOF(1 for x 2 for y), value(N))
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
  return *cantilever;
}

