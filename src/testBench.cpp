#include "testBench.h"
int main() {
  // -- start logging
  Log::Logger().setLevel(Log::LevelInfo);
  // -- GEOMETRY (building the contiliver)
  Geometry p_cantilever = geometry();
  // -- Build material vector (E, nu)
  double* d_matrix;
  cudaMallocManaged(&d_matrix, p_cantilever.get_numberOfElementsG()*6*sizeof(double));
  std::vector<Material> mat;
  mat.reserve(p_cantilever.get_numberOfElementsG());
  for(int iter = 0; iter < p_cantilever.get_numberOfElementsG(); iter++)
  {
    mat.emplace_back(1000.0, 0.3);
    for (int j = 0; j<6;j++)
      d_matrix[iter*6+j] = mat[iter].materialMatrix[j];
  }
  // -- Building stiffness matrix on CPU single core, CPU multiple core & GPU.
  //StiffnessMatrixFirstOrder* stiffMat       = new StiffnessMatrixSingleCPU(mat,p_cantilever,4); 
  // StiffnessMatrixFirstOrder* stiffMatParCPU = new StiffnessMatrixParallelCPU(mat,p_cantilever,4,20);
  StiffnessMatrixFirstOrder* stiffMatGPU    = new StiffnessMatrixGPU(d_matrix,p_cantilever,4);
  //StiffnessMatrixFirstOrder* stiffMatGPU2    = new StiffnessMatrixGPU(mat,p_cantilever,4);
  // -- Calculate the stiffness matrix and do the assembly
  //Sparse &k = stiffMat->GetStiffnessMatrix();
  Sparse &k = stiffMatGPU->GetStiffnessMatrix();
  //Sparse &k2 = stiffMatGPU2->GetStiffnessMatrix();
  //std::cout << Sparse::Compare(k,k2) << " is the Max Error"<<std::endl;
  //k2.STLAssemble2(0.001);
  //Recorder::File().SparseMatrix("assemblyBefore.out", k);
  //k.ThrustAssemble2(0.001);
  //Recorder::File().SparseMatrix("assembly1.out", k);
  
  /* method for doing the assembly in single and multiple cpu */
  //Assembly *a = new AssemblySingleCpu(k);
  //Assembly *a = new AssemblyParCpu(k,20);
  //a->calculateAssembly();
  //Recorder::File().SparseMatrix("assembly2.out", a->stiffMat);
  //delete(a);

  k.ThrustAssemble2(0.001);
  Recorder::File().SparseMatrix("k.out", k);
  // -- Solver
  
  double* displacement;
  cudaMallocManaged(&displacement, k.get_numberOfRows()*sizeof(double));
  cudaMallocManaged(&displacement, k.get_numberOfRows()*sizeof(double));
  //Solver(N, nz, val, rowPtr, colIndex, x, rhs); 
  Solver(k.get_numberOfRows(), k.get_valueSize(), k.get_value(), (int*) k.get_i(), (int*) k.get_j(),
  	   displacement, p_cantilever.get_Load().get_vector());
  // -- Recorder
  Log::Logger().Info(displacement[k.get_numberOfRows()-1]);
  Recorder::File().matrix("displacement.out", displacement, k.get_numberOfRows());
  



  //Recorder::File().SparseMatrix("k2.out", k2);
  //Recorder::File().SparseMatrix("k3.out", k3);
  // -- Delete all variables in Heap memory
  //cudaFree(displacement);
  //delete stiffMat;
  //delete stiffMatParCPU;
  delete stiffMatGPU;
  

  //Sparse kCPU = stiffMatParCPU->GetStiffnessMatrix();
  //
  
  //std::cout << Sparse::Compare(k,kGPU) << " is the Max Error"<<std::endl;
  
  // Log::Logger().Info(f.value[f.valueSize-1]);

  return 0;
}

Geometry& geometry() {
  // GEOMETRY (building the cantilever)
  double dimensionX = 100.0; int numberOfElementX = 80; // dimension of x and number of element in x
  double dimensionY =  10.0; int numberOfElementY = 8; 
  double incr_x = dimensionX/numberOfElementX; // increment between each node 
  double incr_y = dimensionY/numberOfElementY;
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
  cantilever->set_thickness(1.00);

  for (unsigned int j = 0; j <= numberOfElementY; j++)
    cantilever->dof->fix(j*(numberOfElementX+1),1,1);
  cantilever->load->point((numberOfElementX+1)*(numberOfElementY+1),0.0,-0.1); // pointLoad(NodeNumber,DOF(1 for x 2 for y), value(N))
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

