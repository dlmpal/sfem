// Solve the equation of linear elasticity in 3d.

#include "sfem.h"

using namespace sfem;

extern void ElasitcitySolver(Float rho, Float E, Float nu);

int main(int argc, char **argv)
{
    Initialize("ELASTICITY_SOLVER", &argc, &argv);
    Float rho = 1e-3;
    Float E = 1e5;
    Float nu = 0.2;
    ElasitcitySolver(rho, E, nu);
    Finalize();
    return 0;
}

void ElasitcitySolver(Float rho, Float E, Float nu)
{
    auto mesh = io::ReadMesh("Mesh3D");

    // Displacement field
    field::Field U("U", 3, &mesh);
    for (auto i = 0; i < 3; i++)
        U.AddFixedDof("Fixed", i, 0);

    kernel::elasticity::Elasticity3D elasticity_kernel(E, nu);
    kernel::basic::Mass mass_kernel(rho);
    kernel::elasticity::Gravity gravity_kernel(rho, 1);

    assembly::MonoFieldAssembler assembler(&U, assembly::AssemblyType::DYNAMIC);
    assembler.AddKernel("Body", &elasticity_kernel);
    assembler.AddKernel("Body", &mass_kernel);
    assembler.AddKernel("Body", &gravity_kernel);

    solver::StaticSolver solver(&assembler);
    solver.Run();

    io::WriteFieldValues("fields/U_0", U, true);
    // io::WriteVTK("fields/beam_" + std::to_string(::Logger::GetInstance().proc_rank) + ".vtk", mesh, {U});
}