#pragma once

#include "../geo/sfem_geo.h"
#include "basis/basis.h"

namespace sfem::fe
{
    /// @brief Nodal finite flement computations
    class FiniteElement
    {
    public:
        enum class ElementType : int
        {
            POINT_ELEMENT = 0,
            EDGE_ELEMENT = 1,
            SURFACE_ELEMENT = 2,
            VOLUME_ELEMENT = 3
        };

        /// @brief Create a FiniteElement
        /// @param cell Corresponding mesh Cell
        /// @param dim Physical dimension of the mesh
        /// @param n_vars Number of variables per node
        /// @param xpts Cell nodal positions
        FiniteElement(const mesh::Cell &cell, int dim, int n_vars, const std::array<Float, SFEM_MAX_CELL_NODES * 3> &xpts);

        /// @brief Get the element type
        ElementType GetType() const;

        /// @brief Get the element's dimension
        int GetDim() const;

        /// @brief Get the number of variables per node
        int GetNumVars() const;

        /// @brief Get the number of nodes
        int GetNumNodes() const;

        /// @brief Get the number of DOF
        int GetNumDof() const;

        /// @brief Get a reference to the nodal positions vector
        const std::array<Float, SFEM_MAX_CELL_NODES * 3> &GetXpts() const;

        /// @brief Get a reference to the underlying Cell
        const mesh::Cell &GetCell() const;

        /// @brief Get a pointer to the underlying Shape
        const geo::Shape *GetShape() const;

        /// @brief Get a pointer to the underlying Basis
        const basis::Basis *GetBasis() const;

        /// @brief Evaluate the coordinate transformation @ pt */
        void ComputeTransform(const std::array<Float, 3> &pt);

        /// @brief Get the Jacobian determinant (Natural to Physical)
        /// for the latest evaluation point
        Float J() const;

        /// @brief Get the shape function values for the latest evaluation point
        const std::array<Float, SFEM_MAX_CELL_NODES> N() const;

        /// @brief Get the shape function gradient (w.r.t physical coords) for the latest evaluation point
        const std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES> dNdX() const;

        /* Evaluate the field gradient @ pt, given the field values at the nodes */
        // void
        // ComputeFieldGrad(const Float pt[], const Float value_nodes[], Float grad[][3]);

        // /* Interpolate the field value @ pt, given the field value @ element nodes */
        // void InterpValues(const Float pt[], const Float value_nodes[], Float value[]);

    private:
        /// @brief Invert the Natural-to-Physical Jacobian
        Float InvertJac(const Float dXdxi[][3], Float dxidX[][3]);

        /// @brief Element type
        ElementType type;

        /// @brief Element dimension
        int dim;

        /// @brief Number of variables per node
        int n_vars;

        /// @brief Number of nodes
        int n_nodes;

        /// @brief Number of DOF
        int n_dof;

        /// @brief Nodal positions
        std::array<Float, SFEM_MAX_CELL_NODES * 3> xpts;

        /// @brief Underlying mesh Cell
        mesh::Cell cell;

        /// @brief Underlying geometrical Shape
        geo::Shape *shape;

        /// @brief Finite Element basis
        basis::Basis *basis;

        /// @brief Jacobian determinant (Natural to Physical)
        Float _J;

        /// @brief Shape function
        std::array<Float, SFEM_MAX_CELL_NODES> _N;

        /// @brief Shape function gradient (w.r.t natural coords)
        std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES> _dNdxi;

        /// @brief Shape function gradient (w.r.t physical coords)
        std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES> _dNdX;
    };
}