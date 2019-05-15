#ifndef TESTBENCH_H
#define TESTBENCH_H

#include <iostream>
#include "cs.h"
#include "Log/Log.h"
#include "Timer/Timer.h"
#include "Sparse/Sparse.h"
#include "Geometry/Geometry.h"
#include "Material/Material.h"
#include "Element/Element.h"
#include "ElementParallel/ElementParallel.h"
#include "ElementCuda/ElementCuda.h"

int main();
Geometry &geometry();

#endif
