#pragma once

#include "solver.h"

namespace sfem::solver
{
    struct StaticSolver: Solver
    {
        StaticSolver(assembly::MonoFieldAssembler* assembler);

        ~StaticSolver();

        void _Run();

        void PostRun();
    };
}