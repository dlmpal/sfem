#pragma once

#include "elasticity.h"

namespace sfem::kernel::elasticity
{

    struct ShellElement
    {
        /* Underlying geometrical shape */
        const geo::Shape *shape;

        /* Element basis */
        const fe::basis::Basis *basis;

        /* Points (mid-surface) */
        const Float *xpts;

        /* Shell thickness */
        Float thick = 1e-3;

        /* Director vectors */
        geo::Vec3 v1;
        geo::Vec3 v2;
        geo::Vec3 v3;

        /* Dedfined as Gi = 0.5 * thick * vi */
        geo::Vec3 G1;
        geo::Vec3 G2;

        /* Number of quad points */
        int n_qpts;

        /* Quadrature points */
        Float qpts[SFEM_MAX_QUAD_POINTS][3];

        /* Quadrature weights */
        Float qWts[SFEM_MAX_QUAD_POINTS];

        /* Current integration/evaluation point */
        const Float *pt;

        /* Jacobian determinant (Natural to Physical) */
        Float J;

        /* Natural to Physical Jacobian */
        Float dXdxi[3][3];

        /* Physical to Natural Jacobian */
        Float dxidX[3][3];

        /* Shape function */
        Float N[SFEM_MAX_CELL_NODES];

        /* Shape function gradient (w.r.t) natural coordinates */
        Float dNdxi[SFEM_MAX_CELL_NODES][3];

        /* Shape function gradient (w.r.t physical coordinates) */
        Float dNdX[SFEM_MAX_CELL_NODES][3];

        ShellElement(const fe::basis::Basis *basis, const geo::Shape *shape, Float thick, const Float xpts[]);
        void GetQuadrature();
        void GetDirctors();
        void ComputeTransform(const Float pt[]);
        void PhysicalToNaturalJac();
        void NaturalToPhysicalJac();
        void ComputeShapeGrad();
    };

    struct Shell : Elasticity
    {
        /* Shell finite element*/
        ShellElement *shell;

        /* Thickness */
        Float thick;

        Shell(Float E, Float nu, Float thick);
        std::pair<int, int> ShapeDerivMatrixSize() override;
        void ShapeDerivMatrix(Float B[]) override;
        void Constitutive(Float D[]) override;
        void Integrate(fe::FiniteElement *fe, Float kloc[]) override;
    };
}