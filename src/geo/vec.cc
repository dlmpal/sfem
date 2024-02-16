#include "vec.h"

namespace sfem::geo
{
    //=============================================================================
    Vec3::Vec3(Float x, Float y, Float z)
    {
        this->x[0] = x;
        this->x[1] = y;
        this->x[2] = z;
    }
    //=============================================================================
    Vec3::Vec3(Float x1, Float y1, Float z1, Float x2, Float y2, Float z2)
    {
        x[0] = x2 - x1;
        x[1] = y2 - y1;
        x[2] = z2 - z1;
    }
    //=============================================================================
    Vec3 Vec3::operator+(const Vec3 &other) const
    {
        return Vec3{this->x[0] + other.x[0], this->x[1] + other.x[1], this->x[2] + other.x[2]};
    }
    //=============================================================================
    Vec3 Vec3::operator-(const Vec3 &other) const
    {
        return Vec3{this->x[0] - other.x[0], this->x[1] - other.x[1], this->x[2] - other.x[2]};
    }
    //=============================================================================
    Float Vec3::operator*(const Vec3 &other) const
    {
        return this->x[0] * other.x[0] + this->x[1] * other.x[1] + this->x[2] * other.x[2];
    }
    //=============================================================================
    Vec3 Vec3::operator*(Float a) const
    {
        return Vec3{this->x[0] * a, this->x[1] * a, this->x[2] * a};
    }
    //=============================================================================
    Float Vec3::Norm() const
    {
        return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    }
    //=============================================================================
    Vec3 Vec3::CrossProduct(const Vec3 &other) const
    {
        return Vec3{this->x[1] * other.x[2] - this->x[2] * other.x[1],
                    this->x[2] * other.x[0] - this->x[0] * other.x[2],
                    this->x[0] * other.x[1] - this->x[1] * other.x[0]};
    }
    //=============================================================================
    Vec3 Vec3::Normalize() const
    {
        Float mag = this->Norm();
        return Vec3{x[0] / mag, x[1] / mag, x[2] / mag};
    }
    //=============================================================================
    Vec3 Vec3::UnitTangent() const
    {
        return Normalize();
    }
    //=============================================================================
    Vec3 Vec3::UnitNormal() const
    {
        Float mag = this->Norm();
        Float xn = -x[1] / mag;
        Float yn = x[0] / mag;
        Float zn = x[2] / mag;
        return Vec3{xn, yn, zn};
    }

}