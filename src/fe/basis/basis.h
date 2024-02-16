#pragma once

#include "../../mesh/cell.h"

namespace sfem::fe::basis
{
    /// @brief Abstract nodal finite flement basis
    class Basis
    {
    public:
        /// @brief Destructor
        virtual ~Basis() = 0;

        /// @brief Get the dimension
        int GetDim() const;

        /// @brief Get the number of nodes
        int GetNumNodes() const;

        /// @brief Get the number of quadrature points
        int GetNumQuadPts() const;

        /// @brief Get the quadrature weight
        Float GetQuadWeight(int npt) const;

        /// @brief Get the quadrature point coordiniates
        const std::array<Float, 3> &GetQuadPt(int npt) const;

        /// @brief Evaluate shape function and its gradient @ pt
        virtual std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
        ComputeShapeAndGrad(const std::array<Float, 3> &pt) const = 0;

    protected:
        /// @brief Basis dimension
        int dim;

        /// @brief Basis number of nodes
        int n_nodes;

        /// @brief Natural coordintes of reference element nodes
        std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES> xi;

        /// @brief Number of quadrature points
        int n_qpts;

        /// @brief Quadrature weights
        std::array<Float, SFEM_MAX_QUAD_POINTS> qwts;

        /// @brief Quadrature point coordinates
        std::array<std::array<Float, 3>, SFEM_MAX_QUAD_POINTS> qpts;
    };

    /// @brief Point
    struct PointBasis : public Basis
    {
    public:
        PointBasis();
        std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
        ComputeShapeAndGrad(const std::array<Float, 3> &pt) const override;
    };

    /// @brief Linear 2-node line segment
    class L1Basis : public Basis
    {
    public:
        L1Basis();
        std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
        ComputeShapeAndGrad(const std::array<Float, 3> &pt) const override;
    };

    /// @brief  Quadratic 3-node line segment
    class L2Basis : public Basis
    {
    public:
        L2Basis();
        std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
        ComputeShapeAndGrad(const std::array<Float, 3> &pt) const override;
    };

    /// @brief Linear 3-node triangle
    class Tri1Basis : public Basis
    {
    public:
        Tri1Basis();
        std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
        ComputeShapeAndGrad(const std::array<Float, 3> &pt) const override;
    };

    /// @brief  Quadratic 6-node triangle
    class Tri2Basis : public Basis
    {
    public:
        Tri2Basis();
        std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
        ComputeShapeAndGrad(const std::array<Float, 3> &pt) const override;
    };

    /// @brief Linear 4-node tetrahedron
    class Tet1Basis : public Basis
    {
    public:
        Tet1Basis();
        std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
        ComputeShapeAndGrad(const std::array<Float, 3> &pt) const override;
    };

    /// @brief  Quadratic 10-node tetrahedron
    class Tet2Basis : public Basis
    {
    public:
        Tet2Basis();
        std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
        ComputeShapeAndGrad(const std::array<Float, 3> &pt) const override;
    };

    /// @brief
    Basis *CreateBasis(const mesh::Cell &cell);
}