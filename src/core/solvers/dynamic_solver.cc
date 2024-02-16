#include "dynamic_solver.h"

namespace sfem::solver
{
    DynamicSolver::DynamicSolver(assembly::MonoFieldAssembler *assembler, Float t_start, Float t_final, Float dt, int integrator_type)
        : Solver(assembler)
    {
        this->t_start = t_start;

        this->t_final = t_final;

        this->dt = dt;

        this->integrator = time_integrator::CreateIntegrator(assembler, dt, integrator_type);
    }

    DynamicSolver::~DynamicSolver()
    {
        delete integrator;
    }

    int DynamicSolver::GetNumTimeSteps()
    {
        int n_timesteps = static_cast<int>((t_final - t_start) / dt) + 1;

        return n_timesteps;
    }

    void DynamicSolver::ResetTime()
    {
        time = std::make_pair(0, t_start);
    }

    std::pair<int, Float> DynamicSolver::GetCurrentTime()
    {
        return time;
    }

    void DynamicSolver::IncrementTime()
    {
        time.first += 1;

        time.second += dt;
    }

    void DynamicSolver::Run()
    {
        ResetTime();

        auto n_timesteps = GetNumTimeSteps();

        for (auto i = 0; i < n_timesteps; i++)
        {
            PreRun();

            _Run();

            PostRun();
        }
    }

    void DynamicSolver::_Run()
    {
        assembler->operator()();

        integrator->Step();

        assembler->field->UpdateValues(assembler->U);
    }

    void DynamicSolver::PostRun()
    {
        assembler->field->WriteValues(GetCurrentTime().first);

        IncrementTime();
    }
}