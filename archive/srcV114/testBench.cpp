#include "testBench.h"

int main() {
  // -- start logging
  Log::Logger().setLevel(Log::LevelInfo);
  // -- GEOMETRY (building the contiliver)
  Geometry p_cantilever = geometry();
  // -- Build material vector (E, nu)
  Material mat(1000.0f, 0.3f);
  // -- Building stiffness matrix on CPU single core, CPU multiple core & GPU.
  StiffnessMatrixFirstOrder* stiffMat       = new StiffnessMatrixSingleCPU(mat,p_cantilever,4);
  // -- Calculate the stiffness matrix and do the assembly
  Sparse &k = stiffMat->GetStiffnessMatrix();
  Recorder::File().SparseMatrix("kSingleCPUBefor.out", k);
  k.Assemble(0.001);
  Recorder::File().SparseMatrix("kSingleCPU.out", k);
  // -- Solver
  float* displacement;
  cudaMallocManaged(&displacement, k.get_numberOfRows()*sizeof(float));
  SolverSp(k.get_value(), k.get_i(), k.get_j(), k.get_valueSize(), k.get_numberOfRows(),
	   p_cantilever.get_Load().get_vector(), displacement);
  Log::Logger().Info(displacement[k.get_numberOfRows()-1]);
  Recorder::File().matrix("displacement.out", displacement, k.get_numberOfRows());
  
  //StiffnessMatrixFirstOrder* stiffMatParCPU = new StiffnessMatrixParallelCPU(mat,p_cantilever,4,5);
  //StiffnessMatrixFirstOrder* stiffMatGPU    = new StiffnessMatrixGPU(mat,p_cantilever,2);

  //Sparse kCPU = stiffMatParCPU->GetStiffnessMatrix();
  //Sparse kGPU = stiffMatGPU->GetStiffnessMatrix();
  
  
  //std::cout << Sparse::Compare(k,kGPU) << " is the Max Error"<<std::endl;
  
  // Log::Logger().Info(f.value[f.valueSize-1]);

  // -- Delete all variables in Heap memory
  cudaFree(displacement);
  delete stiffMat;
  //delete stiffMatParCPU;
  //delete stiffMatGPU;
  return 0;
}

Geometry& geometry() {
  // GEOMETRY (building the contiliver)
  float dimentionX = 25.0; int numberOfElementX = 5; // dimention of x and number of element in x
  float dimentionY =  1.5; int numberOfElementY = 3; 
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
    cantilever->dof->fix(j*(numberOfElementX+1),1,1);
  cantilever->load->point((numberOfElementX+1)*(numberOfElementY+1),0.0f,-0.1f); // pointLoad(NodeNumber,DOF(1 for x 2 for y), value(N))
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

