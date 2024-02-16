#include "elasticity.h"

namespace sfem::kernel::elasticity
{
    class Elasticity2D : public Elasticity
    {
    public:
        Elasticity2D(Float E, Float nu);

    private:
        std::pair<int, int> ShapeDerivMatrixSize() override;
        void ShapeDerivMatrix(Float B[]) override;
    };

    class PlaneStress : public Elasticity2D
    {
    public:
        PlaneStress(Float E, Float nu);

    private:
        void Constitutive(Float D[]) override;
    };

    class PlaneStrain : public Elasticity2D
    {
    public:
        PlaneStrain(Float E, Float nu);

    private:
        void Constitutive(Float D[]) override;
    };
}