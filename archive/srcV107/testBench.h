#ifndef TESTBENCH_H
#define TESTBENCH_H

#include <iostream>
#include "cs.h"
#include "Log/Log.h"
#include "Timer/Timer.h"
#include "Sparse/Sparse.h"
#include "Geometry/Geometry.h"
#include "Material/Material.h"
#include "StiffnessMatrix/StiffnessMatrix.h"
#include "StiffnessMatrix/StiffnessMatrixFirstOrder/StiffnessMatrixFirstOrder.h"
#include "StiffnessMatrix/StiffnessMatrixParallelCPU/StiffnessMatrixParallelCPU.h"
#include "StiffnessMatrix/StiffnessMatrixSingleCPU/StiffnessMatrixSingleCPU.h"
#include "StiffnessMatrix/StiffnessMatrixGPU/StiffnessMatrixGPU.h"

int main();
Geometry* geometry();
float maxError(StiffnessMatrix& , StiffnessMatrix&);
#endif
