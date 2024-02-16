#include "implicit_euler.h"

namespace sfem::time_integrator
{
    ImplicitEuler::ImplicitEuler(assembly::MonoFieldAssembler *assembler, Float dt) : TimeIntegrator(assembler, dt)
    {
    }

    ImplicitEuler::~ImplicitEuler()
    {
    }

    void ImplicitEuler::Step()
    {
        MatScale(assembler->M->_matA, 1 / dt);

        MatAXPY(assembler->K->_matA, 1.0, assembler->M->_matA, SAME_NONZERO_PATTERN);

        MatMultAdd(assembler->M->_matA, assembler->U->_x, assembler->F->_x, assembler->F->_x);

        la::PetscSolve(assembler->K, assembler->F, assembler->U);
    }
}