#pragma once

#include "../mono_field_assembler.h"

namespace sfem::time_integrator
{
    struct TimeIntegrator
    {
        int type;

        /* Time step size */
        Float dt;

        /* Matrix assembly */
        assembly::MonoFieldAssembler *assembler;

        /* Constructor */
        TimeIntegrator(assembly::MonoFieldAssembler *assembler, Float dt);

        /* Destructor */
        ~TimeIntegrator();

        /* Advance (step) the solution for one timestep */
        virtual void Step() = 0;
    };

    TimeIntegrator* CreateIntegrator(assembly::MonoFieldAssembler* assembler, Float dt, int type);
}

#include "implicit_euler.h"