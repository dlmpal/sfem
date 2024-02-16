#include "finite_element.h"
#include "../common/error.h"
#include "../common/mat_ops.h"

namespace sfem::fe
{
    //=============================================================================
    FiniteElement::FiniteElement(const mesh::Cell &cell, int dim, int n_vars, const std::array<Float, SFEM_MAX_CELL_NODES * 3> &xpts)
        : cell(cell), dim(dim), n_vars(n_vars), xpts(xpts)
    {
        shape = geo::CreateShape(cell);
        basis = basis::CreateBasis(cell);
        n_nodes = basis->GetNumNodes();
        n_dof = n_nodes * n_vars;

        if (dim == basis->GetDim())
        {
            type = ElementType::VOLUME_ELEMENT;
        }
        else
        {
            type = (ElementType)basis->GetDim();
        }

        _J = -1;
        _N.fill(0.0);
        for (auto i = 0; i < n_nodes; i++)
        {
            _dNdX[i].fill(0.0);
            _dNdxi[i].fill(0.0);
        }
    }
    //=============================================================================
    FiniteElement::ElementType FiniteElement::GetType() const
    {
        return type;
    }
    //=============================================================================
    int FiniteElement::GetDim() const
    {
        return dim;
    }
    //=============================================================================
    int FiniteElement::GetNumVars() const
    {
        return n_vars;
    }
    //=============================================================================
    int FiniteElement::GetNumNodes() const
    {
        return n_nodes;
    }
    //=============================================================================
    int FiniteElement::GetNumDof() const
    {
        return n_dof;
    }
    //=============================================================================
    const mesh::Cell &FiniteElement::GetCell() const
    {
        return cell;
    }
    //=============================================================================
    const geo::Shape *FiniteElement::GetShape() const
    {
        return shape;
    }
    //=============================================================================
    const basis::Basis *FiniteElement::GetBasis() const
    {
        return basis;
    }
    //=============================================================================
    const std::array<Float, SFEM_MAX_CELL_NODES * 3> &FiniteElement::GetXpts() const
    {
        return xpts;
    }
    //=============================================================================
    void FiniteElement::ComputeTransform(const std::array<Float, 3> &pt)
    {
        // Evaluate the shape function and its gradients w.r.t natural coords
        std::tie(_N, _dNdxi) = basis->ComputeShapeAndGrad(pt);

        // Natural to Physical (direct)
        Float dXdxi[3][3] = {0};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < n_nodes; k++)
                    dXdxi[i][j] += _dNdxi[k][j] * xpts[k * 3 + i];

        // Physical to Natural (inverse)
        Float dxidX[3][3] = {0};
        _J = InvertJac(dXdxi, dxidX);
        if (_J < 0)
        {
            error::NegativeJacobianError(cell.id, __FILE__, __LINE__);
        }

        // Shape function gradient w.r.t to physical coords
        for (int i = 0; i < n_nodes; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                _dNdX[i][j] = 0.0;
                for (int k = 0; k < 3; k++)
                    _dNdX[i][j] += dxidX[k][j] * _dNdxi[i][k];
            }
        }
    }
    //=============================================================================
    Float FiniteElement::InvertJac(const Float dXdxi[][3], Float dxidX[][3])
    {
        switch (type)
        {
        case ElementType::VOLUME_ELEMENT:
        {
            return mat_ops::inv(basis->GetDim(), dXdxi, dxidX);
        }
        case ElementType::SURFACE_ELEMENT:
        {
            geo::Vec3 j1(dXdxi[0][0], dXdxi[1][0], dXdxi[2][0]);
            geo::Vec3 j2(dXdxi[0][1], dXdxi[1][1], dXdxi[2][1]);
            return j1.CrossProduct(j2).Norm();
        }
        case ElementType::EDGE_ELEMENT:
        {
            return geo::Vec3(xpts[1 * 3 + 0] - xpts[0 * 3 + 0],
                             xpts[1 * 3 + 1] - xpts[0 * 3 + 1],
                             xpts[1 * 3 + 2] - xpts[0 * 3 + 2])
                .Norm();
        }
        // ElementType::POINT_ELEMENT
        default:
            return 1.0;
        }
    }
    //=============================================================================
    Float FiniteElement::J() const
    {
        return _J;
    }
    //=============================================================================
    const std::array<Float, SFEM_MAX_CELL_NODES> FiniteElement::N() const
    {
        return _N;
    }
    //=============================================================================
    const std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES> FiniteElement::dNdX() const
    {
        return _dNdX;
    }

    //=============================================================================
    // void FiniteElement::ComputeFieldGrad(const Float pt[], const Float value_nodes[], Float grad[][3])
    // {
    //     ComputeTransform(pt);
    //     for (auto i = 0; i < n_vars; i++)
    //         for (auto j = 0; j < 3; j++)
    //         {
    //             grad[i][j] = 0;
    //             for (auto k = 0; k < n_nodes; k++)
    //                 grad[i][j] += _dNdX[k][j] * value_nodes[k * n_vars + i];
    //         }
    // }
    //=============================================================================
    // void FiniteElement::InterpValues(const Float pt[], const Float value_nodes[], Float value[])
    // {
    //     ComputeTransform(pt);
    //     for (auto i = 0; i < n_vars; i++)
    //     {
    //         value[i] = 0.0;
    //         for (auto j = 0; j < n_nodes; j++)
    //             value[i] += _N[j] * value_nodes[j * n_vars + i];
    //     }
    // }
}