#pragma once

#include "petsc.h"
#include "../third_party/sfem_petsc.h"

namespace sfem::la
{
    std::pair<std::vector<Float>, std::vector<Vector *>> SlepcSolve(Matrix *A, Matrix *B);
}
