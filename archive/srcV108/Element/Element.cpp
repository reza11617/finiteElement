#include "Element.h"
#define D material->materialMatrix

Element::Element(Material& mat, Geometry& geo, unsigned int n)
  : material(&mat), geometry(&geo), numberOfIntegrationPoint(n)
{

  // single core processing
};

Element::~Element()
{
  Log::Logger().Info("Element Deleted");
  delete[] integrationNode;
  delete[] integrationPos;
  delete[] integrationWeight;
  delete[] DOF_i;
  delete[] DOF_j;
  delete[] stiffMat;
  delete[] c;
}

void Element::Initialize() {
  // Initialize
  numberOfElements = geometry->numberOfElementsG;
  nipSquared = numberOfIntegrationPoint*numberOfIntegrationPoint;
  stiffMatSize = numberOfElements*sizeStiffMatPerEle*nipSquared;
  DOF_i = new unsigned int[stiffMatSize];
  DOF_j = new unsigned int[stiffMatSize];
  stiffMat = new float[stiffMatSize];
  c = new float[numberOfElements*6];
  integrationNode   = new float[numberOfIntegrationPoint];
  integrationPos    = new unsigned int[numberOfIntegrationPoint*dimention*numberOfIntegrationPoint];
  integrationWeight = new float[numberOfIntegrationPoint];

  // fuctions
  integrationPoint();  
}

void Element::SingleCoreSimulation()
{
  Initialize();
  Log::Logger().Info("time for SingleCoreSimulation function");
  Timer timer;
  for (unsigned int i = 0; i<numberOfElements; i++)
    {
      for (unsigned int j = 0; j<nipSquared; j++)
	{
	  stiffnessMatrixFirstOrder(i,j);
	}
    }
}

void Element::stiffnessMatrixFirstOrder(unsigned int numberElement, unsigned int noIP)
{
  if (noIP == 0)
    constantCreator(numberElement);
  DOFCreator(numberElement,noIP);
  double XI = integrationNode[integrationPos[2*noIP]]; double YI = integrationNode[integrationPos[2*noIP+1]];
  // Jacobian
  double J11 = c[numberElement*6+0]*YI-c[numberElement*6+1]; double J12 = c[numberElement*6+3]*YI-c[numberElement*6+4];
  double J21 = c[numberElement*6+0]*XI-c[numberElement*6+2]; double J22 = c[numberElement*6+3]*XI-c[numberElement*6+5];
  double detJ = J11*J22-J12*J21;
  double WeightPerDetJ = (integrationWeight[integrationPos[2*noIP]]*integrationWeight[integrationPos[2*noIP+1]])/detJ;
  // derveativs of the shape function N1x N2x ... N1y N2y ...
  double Ni[8] = {J22*( YI-1)/4 -  J12*( XI-1)/4, J22*(-YI+1)/4 -  J12*(-XI-1)/4,\
		  J22*( YI+1)/4 -  J12*( XI+1)/4, J22*(-YI-1)/4 -  J12*(-XI+1)/4,\
		  J11*( XI-1)/4 -  J21*( YI-1)/4, J11*(-XI-1)/4 -  J21*(-YI+1)/4, \
		  J11*( XI+1)/4 -  J21*( YI+1)/4, J11*(-XI+1)/4 -  J21*(-YI-1)/4};
  // multiplication of the stiffness matrix N1x^2 N1x*N2x ....
  double N[36];
  unsigned int counterN = 0;
  for (unsigned int i = 0; i < 8; i++)
    {
      for (unsigned int j = i; j < 8 ; j++)
	N[counterN++] = Ni[i]*Ni[j];
    };
  // find the position to start filling the stiffness matrix
  unsigned int counter = sizeStiffMatPerEle*(numberElement*numberOfIntegrationPoint*numberOfIntegrationPoint+noIP);
  // writes all 36 components of the 8 by 8 stiffness Matrix considering symmetry
  stiffMat[counter+0]  = WeightPerDetJ*(D[0]*N[0] + 2*D[4]*N[4] + D[2]*N[26]);
  stiffMat[counter+1]  = WeightPerDetJ*(D[4]*N[0] + D[5]*N[26] + D[3]*N[4] + D[2]*N[4]);
  stiffMat[counter+2]  = WeightPerDetJ*(D[2]*N[0] + 2*D[5]*N[4] + D[1]*N[26]);
  stiffMat[counter+3]  = WeightPerDetJ*(D[0]*N[1] + D[4]*N[5] + D[4]*N[11] + D[2]*N[27]);
  stiffMat[counter+4]  = WeightPerDetJ*(D[4]*N[1] + D[3]*N[11] + D[2]*N[5] + D[5]*N[27]);
  stiffMat[counter+5]  = WeightPerDetJ*(D[0]*N[8] + 2*D[4]*N[12] + D[2]*N[30]);
  stiffMat[counter+6]  = WeightPerDetJ*(D[4]*N[1] + D[3]*N[5] + D[2]*N[11] + D[5]*N[27]);
  stiffMat[counter+7]  = WeightPerDetJ*(D[2]*N[1] + D[5]*N[5] + D[5]*N[11] + D[1]*N[27]);
  stiffMat[counter+8]  = WeightPerDetJ*(D[4]*N[8] + D[5]*N[30] + D[3]*N[12] + D[2]*N[12]);
  stiffMat[counter+9]  = WeightPerDetJ*(D[2]*N[8] + 2*D[5]*N[12] + D[1]*N[30]);
  stiffMat[counter+10] = WeightPerDetJ*(D[0]*N[2] + D[4]*N[6] + D[4]*N[17] + D[2]*N[28]);
  stiffMat[counter+11] = WeightPerDetJ*(D[4]*N[2] + D[3]*N[17] + D[2]*N[6] + D[5]*N[28]);
  stiffMat[counter+12] = WeightPerDetJ*(D[0]*N[9] + D[4]*N[13] + D[4]*N[18] + D[2]*N[31]);
  stiffMat[counter+13] = WeightPerDetJ*(D[4]*N[9] + D[3]*N[18] + D[2]*N[13] + D[5]*N[31]);
  stiffMat[counter+14] = WeightPerDetJ*(D[0]*N[15] + 2*D[4]*N[19] + D[2]*N[33]);
  stiffMat[counter+15] = WeightPerDetJ*(D[4]*N[2] + D[3]*N[6] + D[2]*N[17] + D[5]*N[28]);
  stiffMat[counter+16] = WeightPerDetJ*(D[2]*N[2] + D[5]*N[6] + D[5]*N[17] + D[1]*N[28]);
  stiffMat[counter+17] = WeightPerDetJ*(D[4]*N[9] + D[3]*N[13] + D[2]*N[18] + D[5]*N[31]);
  stiffMat[counter+18] = WeightPerDetJ*(D[2]*N[9] + D[5]*N[13] + D[5]*N[18] + D[1]*N[31]);
  stiffMat[counter+19] = WeightPerDetJ*(D[4]*N[15] + D[5]*N[33] + D[3]*N[19] + D[2]*N[19]);
  stiffMat[counter+20] = WeightPerDetJ*(D[2]*N[15] + 2*D[5]*N[19] + D[1]*N[33]);
  stiffMat[counter+21] = WeightPerDetJ*(D[0]*N[3] + D[4]*N[7] + D[4]*N[22] + D[2]*N[29]);
  stiffMat[counter+22] = WeightPerDetJ*(D[4]*N[3] + D[3]*N[22] + D[2]*N[7] + D[5]*N[29]);
  stiffMat[counter+23] = WeightPerDetJ*(D[0]*N[10] + D[4]*N[14] + D[4]*N[23] + D[2]*N[32]);
  stiffMat[counter+24] = WeightPerDetJ*(D[4]*N[10] + D[3]*N[23] + D[2]*N[14] + D[5]*N[32]);
  stiffMat[counter+25] = WeightPerDetJ*(D[0]*N[16] + D[4]*N[20] + D[4]*N[24] + D[2]*N[34]);
  stiffMat[counter+26] = WeightPerDetJ*(D[4]*N[16] + D[3]*N[24] + D[2]*N[20] + D[5]*N[34]);
  stiffMat[counter+27] = WeightPerDetJ*(D[0]*N[21] + 2*D[4]*N[25] + D[2]*N[35]);
  stiffMat[counter+28] = WeightPerDetJ*(D[4]*N[3] + D[3]*N[7] + D[2]*N[22] + D[5]*N[29]);
  stiffMat[counter+29] = WeightPerDetJ*(D[2]*N[3] + D[5]*N[7] + D[5]*N[22] + D[1]*N[29]);
  stiffMat[counter+30] = WeightPerDetJ*(D[4]*N[10] + D[3]*N[14] + D[2]*N[23] + D[5]*N[32]);
  stiffMat[counter+31] = WeightPerDetJ*(D[2]*N[10] + D[5]*N[14] + D[5]*N[23] + D[1]*N[32]);
  stiffMat[counter+32] = WeightPerDetJ*(D[4]*N[16] + D[3]*N[20] + D[2]*N[24] + D[5]*N[34]);
  stiffMat[counter+33] = WeightPerDetJ*(D[2]*N[16] + D[5]*N[20] + D[5]*N[24] + D[1]*N[34]);
  stiffMat[counter+34] = WeightPerDetJ*(D[4]*N[21] + D[5]*N[35] + D[3]*N[25] + D[2]*N[25]);
  stiffMat[counter+35] = WeightPerDetJ*(D[2]*N[21] + 2*D[5]*N[25] + D[1]*N[35]);
};

void Element::constantCreator(unsigned int numberElement)
{
  unsigned int i = numberElement*6;
  c[i++] = (geometry->x[geometry->mesh[numberElement*4+0]] - geometry->x[geometry->mesh[numberElement*4+1]] + geometry->x[geometry->mesh[numberElement*4+2]] - geometry->x[geometry->mesh[numberElement*4+3]])/4;
  c[i++] = (geometry->x[geometry->mesh[numberElement*4+0]] - geometry->x[geometry->mesh[numberElement*4+1]] - geometry->x[geometry->mesh[numberElement*4+2]] + geometry->x[geometry->mesh[numberElement*4+3]])/4;
  c[i++] = (geometry->x[geometry->mesh[numberElement*4+0]] - geometry->x[geometry->mesh[numberElement*4+3]] - geometry->x[geometry->mesh[numberElement*4+2]] + geometry->x[geometry->mesh[numberElement*4+1]])/4;
  c[i++] = (geometry->y[geometry->mesh[numberElement*4+0]] - geometry->y[geometry->mesh[numberElement*4+1]] + geometry->y[geometry->mesh[numberElement*4+2]] - geometry->y[geometry->mesh[numberElement*4+3]])/4;
  c[i++] = (geometry->y[geometry->mesh[numberElement*4+0]] - geometry->y[geometry->mesh[numberElement*4+1]] - geometry->y[geometry->mesh[numberElement*4+2]] + geometry->y[geometry->mesh[numberElement*4+3]])/4;
  c[i++] = (geometry->y[geometry->mesh[numberElement*4+0]] - geometry->y[geometry->mesh[numberElement*4+3]] - geometry->y[geometry->mesh[numberElement*4+2]] + geometry->y[geometry->mesh[numberElement*4+1]])/4;
  // defined the constants c1x to c3y
};

void Element::DOFCreator(unsigned int numberElement, unsigned int noIP)
{
  // size of DOF is like the size of stiffness Matrix
  unsigned int counter = sizeStiffMatPerEle*(numberElement*numberOfIntegrationPoint*numberOfIntegrationPoint+noIP);
  unsigned int xi, xj, yi, yj;
  for (unsigned int i = 0; i<8; i++)
    for (unsigned int j = 0; j<i+1; j++)
      {
	xi = i/2; yi = i%2-1; xj = j/2; yj = j%2-1;
	DOF_i[counter]   = (geometry->mesh[numberElement*4+xi]+1)*2+yi;
	DOF_j[counter++] = (geometry->mesh[numberElement*4+xj]+1)*2+yj;
      }  
};

void Element::integrationPoint()
// Creats the integration points
// XI = integrationNode[integrationPos[i]] YI = integrationNode[integrationPos[i+1]] 
{

  unsigned int counter = 0;
  for (unsigned int i = 0; i < numberOfIntegrationPoint; i++)
    for (unsigned int j = 0; j < numberOfIntegrationPoint; j++)
      {
	integrationPos[counter++] = i;
	integrationPos[counter++] = j;
      };
  
  if (numberOfIntegrationPoint == 1) {
    integrationNode[0] = 0; integrationWeight[0] = 4;
  } else if (numberOfIntegrationPoint == 2) {
    integrationNode[0] = -0.57735; integrationWeight[0] = 1.0;
    integrationNode[1] =  0.57735; integrationWeight[1] = 1.0;
  } else if (numberOfIntegrationPoint == 3) {
    integrationNode[0] = -0.774596; integrationWeight[0] = 0.555556;
    integrationNode[1] =  0.0     ; integrationWeight[1] = 0.888889;
    integrationNode[2] =  0.774596; integrationWeight[2] = 0.555556;
  } else if (numberOfIntegrationPoint == 4) {
    integrationNode[0] = -0.861136; integrationWeight[0] = 0.347855;
    integrationNode[1] = -0.339981; integrationWeight[1] = 0.652145;
    integrationNode[2] =  0.339981; integrationWeight[2] = 0.652145;
    integrationNode[3] =  0.861136; integrationWeight[3] = 0.347855;
  } else if (numberOfIntegrationPoint == 5) {
    integrationNode[0] = -0.90618;  integrationWeight[0] = 0.236927;
    integrationNode[1] = -0.538469; integrationWeight[1] = 0.478629;
    integrationNode[2] =  0.0     ; integrationWeight[2] = 0.568889;
    integrationNode[3] =  0.538469; integrationWeight[3] = 0.478629;
    integrationNode[4] =  0.90618;  integrationWeight[4] = 0.236927;
  } else {
    Log::Logger().Error("Integration points more than five is under construction");
  }
};

unsigned int Element::getStiffMatSize() {
  return stiffMatSize;
};
