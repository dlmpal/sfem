#pragma once

#include "solver.h"

namespace sfem::solver
{
    struct ModalSolver : Solver
    {
        ModalSolver(assembly::MonoFieldAssembler *assembler);

        ~ModalSolver();

        void _Run();
    };
}