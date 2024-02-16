#include "modal_solver.h"

namespace sfem::solver
{
    ModalSolver::ModalSolver(assembly::MonoFieldAssembler *assembler) : Solver(assembler)
    {
    }

    ModalSolver::~ModalSolver()
    {
    }

    void ModalSolver::_Run()
    {

        assembler->operator()();

        auto eigen_pairs = la::SlepcSolve(assembler->K, assembler->M);

        std::string message = "Number of eigen-pairs: " + std::to_string(eigen_pairs.first.size()) + "\n";

        for (auto i = 0; i < eigen_pairs.first.size(); i++)
        {
            message += std::to_string(sqrt(eigen_pairs.first[i]) / 2.0 / M_PI) + " [Hz]\n";

            field::FieldMPI eigen_mode(assembler->field->name + "_mode_" + std::to_string(i), assembler->field->n_vars, assembler->field->mesh);

            eigen_mode.UpdateValues(eigen_pairs.second[i]);

            eigen_mode.WriteValues();

            delete eigen_pairs.second[i];
        }

        if (Logger::GetInstance().proc_rank == SFEM_ROOT)

            Logger::GetInstance().LogMessage(message, Logger::INFO);
    }

}