#pragma once

#include "../common/config.h"

namespace sfem::geo
{
    /// @brief 3D vector class. It can also be used for 2D vector math,
    /// by ignoring the z component.
    struct Vec3
    {
        // Default constructor
        Vec3() = default;

        /// @brief Create a Vec3 equal to [x, y, z]
        Vec3(Float x, Float y, Float z);

        /// @brief Create a Vec3 from (x1, y1, z1) to (x2, y2, z2)
        Vec3(Float x1, Float x2, Float y1, Float y2, Float z1, Float z2);

        /// @brief Add two vectors
        Vec3 operator+(const Vec3 &other) const;

        /// @brief Subtruct two vectors
        Vec3 operator-(const Vec3 &other) const;

        /// @brief Inner product between two vectors
        Float operator*(const Vec3 &other) const;

        /// @brief Return a vector scaled by a
        Vec3 operator*(Float a) const;

        /// @brief Vector euclidean norm
        Float Norm() const;

        /// @brief 3D vector cross product
        Vec3 CrossProduct(const Vec3 &other) const;

        /// @brief Return the normalized vector
        Vec3 Normalize() const;

        /// @brief Return a unit vector in the direction of Vec3
        Vec3 UnitTangent() const;

        /// @brief Returns a unit vector, normal to the original, in the 2D sense
        Vec3 UnitNormal() const;

        /// @brief Vector xyz components
        Float x[3] = {0};
    };
}