#include "femxx/la/slepc.h"

namespace sfem::la
{
    std::pair<std::vector<Float>, std::vector<Vector *>> SlepcSolve(Matrix *A, Matrix *B)
    {
        EPS solver;
        EPSCreate(SFEM_COMM_WORLD, &solver);
        if (B == nullptr)
        {
            EPSSetOperators(solver, A->Get(), nullptr);
            EPSSetProblemType(solver, EPS_HEP);
        }
        else
        {
            EPSSetOperators(solver, A->Get(), B->Get());
            EPSSetProblemType(solver, EPS_GHEP);
        }
        EPSSetFromOptions(solver);
        EPSSolve(solver);

        int n_pairs;
        EPSGetConverged(solver, &n_pairs);
        std::vector<Float> eig_vals(n_pairs);
        std::vector<Vector *> eig_vecs(n_pairs);
        Float eig_real, eig_imag;
        for (auto i = 0; i < n_pairs; i++)
        {
            eig_vecs[i] = new Vector(A->GetLocalSize(), A->GetGlobalSize(), {});
            EPSGetEigenpair(solver, i, &eig_real, nullptr, eig_vecs[i]->Get(), nullptr);
            eig_vals[i] = eig_real;
        }

        EPSDestroy(&solver);
        return std::make_pair(eig_vals, eig_vecs);
    }
}