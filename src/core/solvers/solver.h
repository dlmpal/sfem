#pragma once

#include "../mono_field_assembler.h"

namespace sfem::solver
{
    struct Solver
    {
        /* Matrix assembly */
        assembly::MonoFieldAssembler *assembler;

        Solver(assembly::MonoFieldAssembler *assembler);

        ~Solver();

        virtual void Run();
        virtual void _Run() = 0;
        virtual void PreRun();
        virtual void PostRun();
    };
}