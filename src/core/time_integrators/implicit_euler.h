#pragma once

#include "time_integrator.h"

namespace sfem::time_integrator
{   
    static int const IMPLICIT_EULER = 0;

    struct ImplicitEuler : TimeIntegrator
    {
        ImplicitEuler(assembly::MonoFieldAssembler *assembler, Float dt);
        
        ~ImplicitEuler();

        void Step();
    };
}