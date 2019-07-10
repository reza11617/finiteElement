#include <iostream>
#include "cs.h"
#include "Log/Log.h"
#include "Timer/Timer.h"
#include "Sparse/Sparse.h"
#include "Geometry/Geometry.h"
#include "Material/Material.h"
#include "Element/Element.h"
#include "ElementCuda/ElementCuda.h"

using namespace std;

int main() {
  // start logging
  Log::Logger().setLevel(Log::LevelInfo);
  //Timer timer;
  // build and Mesh the Cantilever beam
  //Geometry (dimention x, dimention y, #ofElement x, #ofElement y) 
  Geometry cantilever(2000000.0f,1000.0f,40000,100); 
  // Build material vector (E, nu)
  Material mat(1000.0f, 0.3f);
  // Element class (Material, Geometry,IP)
  Element ele(mat,cantilever,4);
  // Element class on GPU
  ElementCuda eleCuda(mat,cantilever,4);
  // prints
  Log::Logger().Info(ele.stiffMat,4);
  Log::Logger().Info(eleCuda.stiffMat,4);
  //Sparse::Printer(mat.materialMatrix);
  return 0;
}



