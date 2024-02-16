#pragma once

#include "vec.h"
#include "../mesh/cell.h"

namespace sfem::geo
{
    /// @brief Abstract shape class
    class Shape
    {
    public:
        /// @brief Compute the normal at a given face.
        /// @note Set the face index to -1 to compute the outward normal (only for 1D and 2D shapes).
        virtual Vec3 FaceNormal(int f_idx, const Float xpts[]) const = 0;

        /// @brief Compute the tangent for a given face
        virtual Vec3 FaceTangent(int f_idx, const Float xpts[]) const = 0;
    };

    /// p
    class Point : public Shape
    {
    public:
        Vec3 FaceNormal(int f_idx, const Float xpts[]) const;
        Vec3 FaceTangent(int f_idx, const Float xpts[]) const;
    };

    /// p1---p2
    class Line : public Shape
    {
    public:
        Vec3 FaceNormal(int f_idx, const Float xpts[]) const;
        Vec3 FaceTangent(int f_idx, const Float xpts[]) const;
    };

    //      P3
    //     /  \
    // f3 /    \ f2
    //   /      \
    //  p1_ _ _ _p2
    //      f1
    class Triangle : public Shape
    {
    public:
        Vec3 FaceNormal(int f_idx, const Float xpts[]) const;
        Vec3 FaceTangent(int f_idx, const Float xpts[]) const;
    };

    //       f3
    //    p4_ _ _ p3
    //    |      |
    // f4 |      | f2
    //    |      |
    //   p1_ _ _ p2
    //       f1
    class Quad : public Shape
    {
    public:
        Vec3 FaceNormal(int f_idx, const Float xpts[]) const;
        Vec3 FaceTangent(int f_idx, const Float xpts[]) const;
    };

    ///
    class Tet : public Shape
    {
    public:
        Vec3 FaceNormal(int f_idx, const Float xpts[]) const;
        Vec3 FaceTangent(int f_idx, const Float xpts[]) const;
    };

    ///
    class Hex : public Shape
    {
    public:
        Vec3 FaceNormal(int f_idx, const Float xpts[]) const;
        Vec3 FaceTangent(int f_idx, const Float xpts[]) const;
    };

    /// @brief
    Shape *CreateShape(const mesh::Cell &cell);
}