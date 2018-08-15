#include "StiffnessMatrixFirstOrder.h"

void StiffnessMatrixFirstOrder::constantCreator(unsigned int numberElement, double* c, double* x, double* y, unsigned int* mesh)
{
  unsigned int i = numberElement*6;
  c[i++] = (x[mesh[numberElement*4+0]] - x[mesh[numberElement*4+1]] + x[mesh[numberElement*4+2]] - x[mesh[numberElement*4+3]])/4; 
  c[i++] = (x[mesh[numberElement*4+0]] - x[mesh[numberElement*4+1]] - x[mesh[numberElement*4+2]] + x[mesh[numberElement*4+3]])/4;
  c[i++] = (x[mesh[numberElement*4+0]] - x[mesh[numberElement*4+3]] - x[mesh[numberElement*4+2]] + x[mesh[numberElement*4+1]])/4;
  c[i++] = (y[mesh[numberElement*4+0]] - y[mesh[numberElement*4+1]] + y[mesh[numberElement*4+2]] - y[mesh[numberElement*4+3]])/4;
  c[i++] = (y[mesh[numberElement*4+0]] - y[mesh[numberElement*4+1]] - y[mesh[numberElement*4+2]] + y[mesh[numberElement*4+3]])/4;
  c[i++] = (y[mesh[numberElement*4+0]] - y[mesh[numberElement*4+3]] - y[mesh[numberElement*4+2]] + y[mesh[numberElement*4+1]])/4;
  // defined the constants c1x to c3y
};

void StiffnessMatrixFirstOrder::stiffnessMatrixCalculation(unsigned int numberElement, unsigned int nip ,double* in, unsigned int* ip, double* iw, double* c, double* D, unsigned int* mesh, double* k, unsigned int* i_index, unsigned int* j_index, unsigned int* dofFree)
// numberElement -> the element number needed to be calculated
// nip is the number of integration point squared.
// in is the integrationNode
// ip -> integrationPos
// iw -> integrationWeight
// c -> constants
// mesh -> a pointer to mesh array 
// D -> material matrix
// k -> pointer to stiffness matrix
// i_index -> pointer to DOFi
// j_index -> pointer to DOFj
// dofFree -> lists the free dofs and value for new Dofs
{
  // define local stiffness Matrix
  double kLocal[36] = {}; 
  for (unsigned int noIP = 0; noIP < nip; noIP++) { 
    // noIP -> the integration point number needed to be calculated
    double XI = in[ip[2*noIP]]; double YI = in[ip[2*noIP+1]];
    // Jacobian
    double J11 = c[numberElement*6+0]*YI-c[numberElement*6+1]; double J12 = c[numberElement*6+3]*YI-c[numberElement*6+4];
    double J21 = c[numberElement*6+0]*XI-c[numberElement*6+2]; double J22 = c[numberElement*6+3]*XI-c[numberElement*6+5];
    double detJ = J11*J22-J12*J21;
    double WeightPerDetJ = (iw[ip[2*noIP]]*iw[ip[2*noIP+1]])/detJ;
    // derveativs of the shape function N1x N2x ... N1y N2y ...
    double Ni[8] = {J22*( YI-1)/4 -  J12*( XI-1)/4, J22*(-YI+1)/4 -  J12*(-XI-1)/4, \
		    J22*( YI+1)/4 -  J12*( XI+1)/4, J22*(-YI-1)/4 -  J12*(-XI+1)/4, \
		    J11*( XI-1)/4 -  J21*( YI-1)/4, J11*(-XI-1)/4 -  J21*(-YI+1)/4, \
		    J11*( XI+1)/4 -  J21*( YI+1)/4, J11*(-XI+1)/4 -  J21*(-YI-1)/4};
    // multiplication of shape functions N1x^2 N1x*N2x ....
    double N[36];
    unsigned int counterN = 0;
    for (unsigned int i = 0; i < 8; i++) {
      for (unsigned int j = i; j < 8 ; j++)
	N[counterN++] = Ni[i]*Ni[j];
    };
    
    // find the position to start filling the stiffness matrix
    // writes all 36 components of the 8 by 8 symmteric stiffness Matrix (storing upper triangle)
    kLocal[0]  = kLocal[0]  + WeightPerDetJ*(D[0]*N[0] + 2*D[4]*N[4] + D[2]*N[26]); //1,1

    kLocal[1]  = kLocal[1]  + WeightPerDetJ*(D[4]*N[0] + D[5]*N[26] + D[3]*N[4] + D[2]*N[4]); //2,1
    kLocal[2]  = kLocal[2]  + WeightPerDetJ*(D[2]*N[0] + 2*D[5]*N[4] + D[1]*N[26]);           //2,2

    kLocal[3]  = kLocal[3]  + WeightPerDetJ*(D[0]*N[1] + D[4]*N[5] + D[4]*N[11] + D[2]*N[27]);//3,1
    kLocal[4]  = kLocal[4]  + WeightPerDetJ*(D[4]*N[1] + D[3]*N[11] + D[2]*N[5] + D[5]*N[27]);//3,2
    kLocal[5]  = kLocal[5]  + WeightPerDetJ*(D[0]*N[8] + 2*D[4]*N[12] + D[2]*N[30]);          //3,3

    kLocal[6]  = kLocal[6]  + WeightPerDetJ*(D[4]*N[1] + D[3]*N[5] + D[2]*N[11] + D[5]*N[27]);//4,1
    kLocal[7]  = kLocal[7]  + WeightPerDetJ*(D[2]*N[1] + D[5]*N[5] + D[5]*N[11] + D[1]*N[27]);//4,2
    kLocal[8]  = kLocal[8]  + WeightPerDetJ*(D[4]*N[8] + D[5]*N[30] + D[3]*N[12] + D[2]*N[12]);//4,3
    kLocal[9]  = kLocal[9]  + WeightPerDetJ*(D[2]*N[8] + 2*D[5]*N[12] + D[1]*N[30]);           //4,4

    kLocal[10] = kLocal[10] + WeightPerDetJ*(D[0]*N[2] + D[4]*N[6] + D[4]*N[17] + D[2]*N[28]);//5,1
    kLocal[11] = kLocal[11] + WeightPerDetJ*(D[4]*N[2] + D[3]*N[17] + D[2]*N[6] + D[5]*N[28]);//5,2
    kLocal[12] = kLocal[12] + WeightPerDetJ*(D[0]*N[9] + D[4]*N[13] + D[4]*N[18] + D[2]*N[31]);//5,3
    kLocal[13] = kLocal[13] + WeightPerDetJ*(D[4]*N[9] + D[3]*N[18] + D[2]*N[13] + D[5]*N[31]);//5,4
    kLocal[14] = kLocal[14] + WeightPerDetJ*(D[0]*N[15] + 2*D[4]*N[19] + D[2]*N[33]);//5,5

    kLocal[15] = kLocal[15] + WeightPerDetJ*(D[4]*N[2] + D[3]*N[6] + D[2]*N[17] + D[5]*N[28]);//6,1
    kLocal[16] = kLocal[16] + WeightPerDetJ*(D[2]*N[2] + D[5]*N[6] + D[5]*N[17] + D[1]*N[28]);//6,2
    kLocal[17] = kLocal[17] + WeightPerDetJ*(D[4]*N[9] + D[3]*N[13] + D[2]*N[18] + D[5]*N[31]);//6,3
    kLocal[18] = kLocal[18] + WeightPerDetJ*(D[2]*N[9] + D[5]*N[13] + D[5]*N[18] + D[1]*N[31]);//6,4
    kLocal[19] = kLocal[19] + WeightPerDetJ*(D[4]*N[15] + D[5]*N[33] + D[3]*N[19] + D[2]*N[19]);//6,5
    kLocal[20] = kLocal[20] + WeightPerDetJ*(D[2]*N[15] + 2*D[5]*N[19] + D[1]*N[33]);//6,6

    kLocal[21] = kLocal[21]  + WeightPerDetJ*(D[0]*N[3] + D[4]*N[7] + D[4]*N[22] + D[2]*N[29]); //7,1
    kLocal[22] = kLocal[22] + WeightPerDetJ*(D[4]*N[3] + D[3]*N[22] + D[2]*N[7] + D[5]*N[29]); //7,2
    kLocal[23] = kLocal[23] + WeightPerDetJ*(D[0]*N[10] + D[4]*N[14] + D[4]*N[23] + D[2]*N[32]);//7,3
    kLocal[24] = kLocal[24] + WeightPerDetJ*(D[4]*N[10] + D[3]*N[23] + D[2]*N[14] + D[5]*N[32]);//7,4
    kLocal[25] = kLocal[25] + WeightPerDetJ*(D[0]*N[16] + D[4]*N[20] + D[4]*N[24] + D[2]*N[34]);//7,5
    kLocal[26] = kLocal[26] + WeightPerDetJ*(D[4]*N[16] + D[3]*N[24] + D[2]*N[20] + D[5]*N[34]);//7,6
    kLocal[27] = kLocal[27] + WeightPerDetJ*(D[0]*N[21] + 2*D[4]*N[25] + D[2]*N[35]);//7,7

    kLocal[28] = kLocal[28] + WeightPerDetJ*(D[4]*N[3] + D[3]*N[7] + D[2]*N[22] + D[5]*N[29]);//8,1
    kLocal[29] = kLocal[29] + WeightPerDetJ*(D[2]*N[3] + D[5]*N[7] + D[5]*N[22] + D[1]*N[29]);//8,2
    kLocal[30] = kLocal[30] + WeightPerDetJ*(D[4]*N[10] + D[3]*N[14] + D[2]*N[23] + D[5]*N[32]);//8,3
    kLocal[31] = kLocal[31] + WeightPerDetJ*(D[2]*N[10] + D[5]*N[14] + D[5]*N[23] + D[1]*N[32]);//8,4
    kLocal[32] = kLocal[32] + WeightPerDetJ*(D[4]*N[16] + D[3]*N[20] + D[2]*N[24] + D[5]*N[34]);//8,5
    kLocal[33] = kLocal[33] + WeightPerDetJ*(D[2]*N[16] + D[5]*N[20] + D[5]*N[24] + D[1]*N[34]);//8,6
    kLocal[34] = kLocal[34] + WeightPerDetJ*(D[4]*N[21] + D[5]*N[35] + D[3]*N[25] + D[2]*N[25]);//8,7
    kLocal[35] = kLocal[35] + WeightPerDetJ*(D[2]*N[21] + 2*D[5]*N[25] + D[1]*N[35]);//8,8
  }
  unsigned int counter = 36*(numberElement);
  unsigned int kcounter = 0, dof_i,dof_j;
  for (unsigned int i = 0; i < 8; i++) {
    for (unsigned int j = 0; j <= i; j++) {
      dof_i = (mesh[numberElement*4+i/2])*2+i%2+1; // this (mesh[numberElement*4+i/2])*2+i%2+1 the required dof however we need to look at dofFree which the index starts at zero
      dof_j = (mesh[numberElement*4+j/2])*2+j%2+1;
      if (dofFree[dof_i-1] && dofFree[dof_j-1]) {
	i_index[counter] = dofFree[max(dof_i-1,dof_j-1)];
	j_index[counter] = dofFree[min(dof_i-1,dof_j-1)];
	k[counter] = kLocal[kcounter];
      }
      //rowPtr[i_index[counter]] = rowPtr[i_index[counter]] + 1;
      counter++;
      kcounter++;
    }
  }
}

StiffnessMatrixFirstOrder::StiffnessMatrixFirstOrder(Material& mat, Geometry& geo, unsigned int n)
  :StiffnessMatrix(mat,geo,n)
{
  Log::Logger().Info("StiffnessMatrixFirstOrder Created by CPU");
  sizeStiffMatPerEle = 36;
  stiffMatSize = numberOfElements*sizeStiffMatPerEle;
  globalStiffMatSize = globalStiffMatSizeCalculator(geometry->get_mesh(), geometry->get_mesh_Size(), geometry->get_x_y_size());
  simulationSize = numberOfElements;
  stiffMat = new Sparse(stiffMatSize,geometry->get_Dof().get_freeSize());
  //  stiffMat = new Sparse(stiffMatSize,geometry->get_x_y_size()*2);
  cudaMallocManaged(&c,numberOfElements*6*sizeof(double));
};

StiffnessMatrixFirstOrder::~StiffnessMatrixFirstOrder() {
  Log::Logger().Info("StiffnessMatrixFirstOrder Deleted by CPU");
  cudaFree(c);
  delete stiffMat;
}

int StiffnessMatrixFirstOrder::GetStiffnessMatrixSize()
{
  return stiffMatSize;
}















