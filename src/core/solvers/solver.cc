#include "solver.h"

namespace sfem::solver
{
    Solver::Solver(assembly::MonoFieldAssembler *assembler)
    {
        this->assembler = assembler;
    }

    Solver::~Solver()
    {
    }

    void Solver::Run()
    {
        PreRun();

        _Run();

        PostRun();
    }

    void Solver::PreRun()
    {
    }

    void Solver::PostRun()
    {
    }
}