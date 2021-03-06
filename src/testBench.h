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
#include "Recorder/Recorder.h"
#include "SolverSp/SolverSp.h"
#include "Assembly/Assembly.h"
#include "Assembly/AssemblySingleCpu.h"
#include "Assembly/AssemblyParCpu.h"
#include "Solver/Solver.h"


int main();
Geometry& geometry();
double maxError(StiffnessMatrix& , StiffnessMatrix&);
#endif
