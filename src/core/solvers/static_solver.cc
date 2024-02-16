#include "static_solver.h"
#include "chrono"

namespace sfem::solver
{
    StaticSolver::StaticSolver(assembly::MonoFieldAssembler *assembler) : Solver(assembler)
    {
    }

    StaticSolver::~StaticSolver()
    {
    }

    void StaticSolver::_Run()
    {
        auto t1 = std::chrono::high_resolution_clock::now();
        assembler->AssembleSystem();
        auto t2 = std::chrono::high_resolution_clock::now();
        auto dur1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        t1 = std::chrono::high_resolution_clock::now();
        la::PetscSolve(assembler->GetStiffnessMatrix(),
                       assembler->GetLoadVector(),
                       assembler->GetSolutionVector());
        t2 = std::chrono::high_resolution_clock::now();
        auto dur2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << dur1.count() << " " << dur2.count() << "\n";
        assembler->UpdateFieldValues();
    }

    void StaticSolver::PostRun()
    {
    }
}