#include "StiffnessMatrixFirstOrder.h"

void StiffnessMatrixFirstOrder::constantCreator(unsigned int numberElement, float* c, float* x, float* y, unsigned int* mesh)
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

void StiffnessMatrixFirstOrder::stiffnessMatrixCalculation(unsigned int numberElement, unsigned int nip ,float* in, unsigned int* ip, float* iw, float* c, float* D, unsigned int* mesh, float* k, unsigned int* i_index, unsigned int* j_index)
// numberElement -> the element number needed to be calculated
// nip is the number of integration point squared.
// in is the integrationNode
// ip -> integrationPos
// iw -> integrationWeight
// c -> constants
// D -> material matrix
// k -> stiffness matrix
{
  unsigned int counter = 36*(numberElement);
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
    for (unsigned int i = 0; i < 8; i++)
      {
	for (unsigned int j = i; j < 8 ; j++)
	  N[counterN++] = Ni[i]*Ni[j];
      };
    // find the position to start filling the stiffness matrix
    // writes all 36 components of the 8 by 8 stiffness Matrix considering symmetry
    k[counter+0]  = k[counter+0]  + WeightPerDetJ*(D[0]*N[0] + 2*D[4]*N[4] + D[2]*N[26]);
    k[counter+1]  = k[counter+1]  + WeightPerDetJ*(D[4]*N[0] + D[5]*N[26] + D[3]*N[4] + D[2]*N[4]);
    k[counter+2]  = k[counter+2]  + WeightPerDetJ*(D[2]*N[0] + 2*D[5]*N[4] + D[1]*N[26]);
    k[counter+3]  = k[counter+3]  + WeightPerDetJ*(D[0]*N[1] + D[4]*N[5] + D[4]*N[11] + D[2]*N[27]);
    k[counter+4]  = k[counter+4]  + WeightPerDetJ*(D[4]*N[1] + D[3]*N[11] + D[2]*N[5] + D[5]*N[27]);
    k[counter+5]  = k[counter+5]  + WeightPerDetJ*(D[0]*N[8] + 2*D[4]*N[12] + D[2]*N[30]);
    k[counter+6]  = k[counter+6]  + WeightPerDetJ*(D[4]*N[1] + D[3]*N[5] + D[2]*N[11] + D[5]*N[27]);
    k[counter+7]  = k[counter+7]  + WeightPerDetJ*(D[2]*N[1] + D[5]*N[5] + D[5]*N[11] + D[1]*N[27]);
    k[counter+8]  = k[counter+8]  + WeightPerDetJ*(D[4]*N[8] + D[5]*N[30] + D[3]*N[12] + D[2]*N[12]);
    k[counter+9]  = k[counter+9]  + WeightPerDetJ*(D[2]*N[8] + 2*D[5]*N[12] + D[1]*N[30]);
    k[counter+10] = k[counter+10] + WeightPerDetJ*(D[0]*N[2] + D[4]*N[6] + D[4]*N[17] + D[2]*N[28]);
    k[counter+11] = k[counter+11] + WeightPerDetJ*(D[4]*N[2] + D[3]*N[17] + D[2]*N[6] + D[5]*N[28]);
    k[counter+12] = k[counter+12] + WeightPerDetJ*(D[0]*N[9] + D[4]*N[13] + D[4]*N[18] + D[2]*N[31]);
    k[counter+13] = k[counter+13] + WeightPerDetJ*(D[4]*N[9] + D[3]*N[18] + D[2]*N[13] + D[5]*N[31]);
    k[counter+14] = k[counter+14] + WeightPerDetJ*(D[0]*N[15] + 2*D[4]*N[19] + D[2]*N[33]);
    k[counter+15] = k[counter+15] + WeightPerDetJ*(D[4]*N[2] + D[3]*N[6] + D[2]*N[17] + D[5]*N[28]);
    k[counter+16] = k[counter+16] + WeightPerDetJ*(D[2]*N[2] + D[5]*N[6] + D[5]*N[17] + D[1]*N[28]);
    k[counter+17] = k[counter+17] + WeightPerDetJ*(D[4]*N[9] + D[3]*N[13] + D[2]*N[18] + D[5]*N[31]);
    k[counter+18] = k[counter+18] + WeightPerDetJ*(D[2]*N[9] + D[5]*N[13] + D[5]*N[18] + D[1]*N[31]);
    k[counter+19] = k[counter+19] + WeightPerDetJ*(D[4]*N[15] + D[5]*N[33] + D[3]*N[19] + D[2]*N[19]);
    k[counter+20] = k[counter+20] + WeightPerDetJ*(D[2]*N[15] + 2*D[5]*N[19] + D[1]*N[33]);
    k[counter+21] = k[counter+21] + WeightPerDetJ*(D[0]*N[3] + D[4]*N[7] + D[4]*N[22] + D[2]*N[29]);
    k[counter+22] = k[counter+22] + WeightPerDetJ*(D[4]*N[3] + D[3]*N[22] + D[2]*N[7] + D[5]*N[29]);
    k[counter+23] = k[counter+23] + WeightPerDetJ*(D[0]*N[10] + D[4]*N[14] + D[4]*N[23] + D[2]*N[32]);
    k[counter+24] = k[counter+24] + WeightPerDetJ*(D[4]*N[10] + D[3]*N[23] + D[2]*N[14] + D[5]*N[32]);
    k[counter+25] = k[counter+25] + WeightPerDetJ*(D[0]*N[16] + D[4]*N[20] + D[4]*N[24] + D[2]*N[34]);
    k[counter+26] = k[counter+26] + WeightPerDetJ*(D[4]*N[16] + D[3]*N[24] + D[2]*N[20] + D[5]*N[34]);
    k[counter+27] = k[counter+27] + WeightPerDetJ*(D[0]*N[21] + 2*D[4]*N[25] + D[2]*N[35]);
    k[counter+28] = k[counter+28] + WeightPerDetJ*(D[4]*N[3] + D[3]*N[7] + D[2]*N[22] + D[5]*N[29]);
    k[counter+29] = k[counter+29] + WeightPerDetJ*(D[2]*N[3] + D[5]*N[7] + D[5]*N[22] + D[1]*N[29]);
    k[counter+30] = k[counter+30] + WeightPerDetJ*(D[4]*N[10] + D[3]*N[14] + D[2]*N[23] + D[5]*N[32]);
    k[counter+31] = k[counter+31] + WeightPerDetJ*(D[2]*N[10] + D[5]*N[14] + D[5]*N[23] + D[1]*N[32]);
    k[counter+32] = k[counter+32] + WeightPerDetJ*(D[4]*N[16] + D[3]*N[20] + D[2]*N[24] + D[5]*N[34]);
    k[counter+33] = k[counter+33] + WeightPerDetJ*(D[2]*N[16] + D[5]*N[20] + D[5]*N[24] + D[1]*N[34]);
    k[counter+34] = k[counter+34] + WeightPerDetJ*(D[4]*N[21] + D[5]*N[35] + D[3]*N[25] + D[2]*N[25]);
    k[counter+35] = k[counter+35] + WeightPerDetJ*(D[2]*N[21] + 2*D[5]*N[25] + D[1]*N[35]);
}
  unsigned int count = counter;
  unsigned int xi, xj, yi, yj;
  unsigned int dof_i, dof_j;
  for (unsigned int i = 0; i<8; i++)
    for (unsigned int j = 0; j<i+1; j++) {
      xi = i/2; yi = i%2-1; xj = j/2; yj = j%2-1;
      dof_i = (mesh[numberElement*4+xi]+1)*2+yi;
      dof_j = (mesh[numberElement*4+xj]+1)*2+yj;
      i_index[count]   = max(dof_i,dof_j);
      j_index[count++] = min(dof_i,dof_j);

    }
  /*
  unsigned int count = counter;
  for (unsigned int i = 0; i<4; i++) {
    for (unsigned int j = 0; j<i+1; j++) {
      for (unsigned int k = 1; k < 3; k++) {
	dof_i = (mesh[numberElement*4+i])*2+k;
	dof_j = (mesh[numberElement*4+j])*2+k;
      }
    }
  } 
  */
}

Sparse* StiffnessMatrixFirstOrder::assembler(Sparse* k, unsigned int* freeDofs, unsigned int freeDofSize)
{
  unsigned int kg_size = k->get_numberOfRows();
  Sparse* K = new Sparse(kg_size*kg_size,kg_size);
  unsigned int posGlobal;
  // assemble value
  int totalElement = k->get_valueSize()/36;
  unsigned int totalDOF = k->get_numberOfRows();
  int counter = 0;
  for (unsigned int nOfE = 0; nOfE< totalElement; nOfE++) {
    counter = nOfE*36;
    for (unsigned int c = counter; c < counter+36; c++) {
      posGlobal = totalDOF*(k->i[c]-1)+k->j[c]-1;
      K->value[posGlobal] = K->value[posGlobal] + k->value[c];
      K->i[posGlobal] = k->i[c];
      K->j[posGlobal] = k->j[c];
      if (k->j[c] != k->i[c]) {
	posGlobal = totalDOF*(k->j[c]-1)+k->i[c]-1;
	K->value[posGlobal] = K->value[posGlobal] + k->value[c];
	K->i[posGlobal] = k->i[c];
	K->j[posGlobal] = k->j[c];
      }
    }
  }
  /*
  // assemble Dofs
  counter = 0;
  for (unsigned int i = 0; i < freeDofSize; i++) {
    for (unsigned int j = 0; j < freeDofSize; j++) {
      K->i[counter]       = freeDofs[i];
      K->j[counter]       = freeDofs[j];
      K->value[counter] = Ktemp[totalNode*2*(K->i[counter]-1)+K->j[counter]-1];
      //Log::Logger().Info(totalNode*2*(K->i[counter]-1)+K->j[counter]-1);
      counter++;
    }
  }
  */
  delete k;
  //delete[] Ktemp;
  return K;
}
  
StiffnessMatrixFirstOrder::StiffnessMatrixFirstOrder(Material& mat, Geometry& geo, unsigned int n)
  :StiffnessMatrix(mat,geo,n)
{
  Log::Logger().Info("StiffnessMatrixFirstOrder Created by CPU");
  sizeStiffMatPerEle = 36;
  stiffMatSize = numberOfElements*sizeStiffMatPerEle;
  globalStiffMatSize = globalStiffMatSizeCalculator(geometry->get_mesh(), geometry->get_mesh_Size(), geometry->get_x_y_size());
  simulationSize = numberOfElements;
  stiffMat = new Sparse(stiffMatSize,geometry->get_x_y_size()*2);
  cudaMallocManaged(&c,numberOfElements*6*sizeof(float));
};

StiffnessMatrixFirstOrder::~StiffnessMatrixFirstOrder()
{
  Log::Logger().Info("StiffnessMatrixFirstOrder Deleted by CPU");
  cudaFree(c);
}

int StiffnessMatrixFirstOrder::GetStiffnessMatrixSize()
{
  return stiffMatSize;
}















