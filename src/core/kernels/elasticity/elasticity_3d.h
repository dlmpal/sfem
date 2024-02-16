#include "elasticity.h"

namespace sfem::kernel::elasticity
{
    class Elasticity3D : public Elasticity
    {
    public:
        Elasticity3D(Float E, Float nu);

    private:
        std::pair<int, int> ShapeDerivMatrixSize() override;
        void ShapeDerivMatrix(Float B[]) override;
        void Constitutive(Float D[]) override;
    };
}