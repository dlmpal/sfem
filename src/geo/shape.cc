#include "shape.h"
#include "../common/error.h"

namespace sfem::geo
{
    //=============================================================================
    Vec3 Point::FaceNormal(int f_idx, const Float xpts[]) const
    {
        return Vec3(1, 1, 1);
    }
    //=============================================================================
    Vec3 Point::FaceTangent(int f_idx, const Float xpts[]) const
    {
        return Vec3(1, 1, 1);
    }
    //=============================================================================
    Vec3 Line::FaceNormal(int f_idx, const Float xpts[]) const
    {
        switch (f_idx)
        {
        case 0:
            return Vec3(xpts[0], xpts[1], xpts[2], xpts[3], xpts[4], xpts[5]).UnitTangent();
        case 1:
            return Vec3(xpts[3], xpts[4], xpts[5], xpts[0], xpts[1], xpts[2]).UnitTangent();
        case -1:
            return Vec3(xpts[0], xpts[1], xpts[2], xpts[3], xpts[4], xpts[5]).UnitNormal();
        default:
            error::InvalidFaceError("Line", f_idx, __FILE__, __LINE__);
        }
    }
    //=============================================================================
    Vec3 Line::FaceTangent(int f_idx, const Float xpts[]) const
    {
        switch (f_idx)
        {
        case 0:
            return Vec3(xpts[0], xpts[1], xpts[2], xpts[3], xpts[4], xpts[5]).UnitTangent();
        case 1:
            return Vec3(xpts[3], xpts[4], xpts[5], xpts[0], xpts[1], xpts[2]).UnitTangent();
        default:
            error::InvalidFaceError("Line", f_idx, __FILE__, __LINE__);
        }
    }
    //=============================================================================
    Vec3 Triangle::FaceNormal(int f_idx, const Float xpts[]) const
    {
        switch (f_idx)
        {
        case 0:
            return Vec3(xpts[0], xpts[1], xpts[2], xpts[3], xpts[4], xpts[5]).UnitNormal();
        case 1:
            return Vec3(xpts[3], xpts[4], xpts[5], xpts[6], xpts[7], xpts[8]).UnitNormal();
        case 2:
            return Vec3(xpts[6], xpts[7], xpts[8], xpts[0], xpts[1], xpts[2]).UnitNormal();
        case -1:
        {
            Vec3 v1(xpts[0], xpts[1], xpts[2], xpts[3], xpts[4], xpts[5]);
            Vec3 v2(xpts[0], xpts[1], xpts[2], xpts[6], xpts[7], xpts[8]);
            return v1.CrossProduct(v2).Normalize();
        }
        default:
            error::InvalidFaceError("Triangle", f_idx, __FILE__, __LINE__);
        }
    }
    //=============================================================================
    Vec3 Triangle::FaceTangent(int f_idx, const Float xpts[]) const
    {
        switch (f_idx)
        {
        case 0:
            return Vec3(xpts[0], xpts[1], xpts[2], xpts[3], xpts[4], xpts[5]).UnitTangent();
        case 1:
            return Vec3(xpts[3], xpts[4], xpts[5], xpts[6], xpts[7], xpts[8]).UnitTangent();
        case 2:
            return Vec3(xpts[6], xpts[7], xpts[8], xpts[0], xpts[1], xpts[2]).UnitTangent();
        default:
            error::InvalidFaceError("Triangle", f_idx, __FILE__, __LINE__);
        }
    }
    //=============================================================================
    Vec3 Quad::FaceNormal(int f_idx, const Float xpts[]) const
    {
        switch (f_idx)
        {
        case 0:
            return Vec3(xpts[0], xpts[1], xpts[2], xpts[3], xpts[4], xpts[5]).UnitNormal();
        case 1:
            return Vec3(xpts[3], xpts[4], xpts[5], xpts[6], xpts[7], xpts[8]).UnitNormal();
        case 2:
            return Vec3(xpts[6], xpts[7], xpts[8], xpts[9], xpts[10], xpts[11]).UnitNormal();
        case 3:
            return Vec3(xpts[9], xpts[10], xpts[11], xpts[0], xpts[1], xpts[2]).UnitNormal();
        case -1:
        {
            Vec3 v1(xpts[0], xpts[1], xpts[2], xpts[3], xpts[4], xpts[5]);
            Vec3 v2(xpts[0], xpts[1], xpts[2], xpts[6], xpts[7], xpts[8]);
            return v1.CrossProduct(v2).Normalize();
        }
        default:
            error::InvalidFaceError("Quad", f_idx, __FILE__, __LINE__);
        }
    }
    //=============================================================================
    Vec3 Quad::FaceTangent(int f_idx, const Float xpts[]) const
    {
        switch (f_idx)
        {
        case 0:
            return Vec3(xpts[0], xpts[1], xpts[2], xpts[3], xpts[4], xpts[5]).UnitTangent();
        case 1:
            return Vec3(xpts[3], xpts[4], xpts[5], xpts[6], xpts[7], xpts[8]).UnitTangent();
        case 2:
            return Vec3(xpts[6], xpts[7], xpts[8], xpts[9], xpts[10], xpts[11]).UnitTangent();
        case 3:
            return Vec3(xpts[9], xpts[10], xpts[11], xpts[0], xpts[1], xpts[2]).UnitTangent();
        default:
            error::InvalidFaceError("Quad", f_idx, __FILE__, __LINE__);
        }
    }
    //=============================================================================
    Vec3 Tet::FaceNormal(int f_idxm, const Float xpts[]) const
    {
    }
    //=============================================================================
    Vec3 Tet::FaceTangent(int f_idx, const Float xpts[]) const
    {
    }
    //=============================================================================
    Vec3 Hex::FaceNormal(int f_idx, const Float xpts[]) const
    {
    }
    //=============================================================================
    Vec3 Hex::FaceTangent(int f_idx, const Float xpts[]) const
    {
    }
    //=============================================================================
    Shape *CreateShape(const mesh::Cell &cell)
    {
        Shape *shape = nullptr;
        switch (cell.type)
        {
        case mesh::CellType::POINT:
            shape = new Point();
            break;
        case mesh::CellType::LINE:
            shape = new Line();
            break;
        case mesh::CellType::TRIANGLE:
            shape = new Triangle();
            break;
        case mesh::CellType::QUAD:
            shape = new Quad();
            break;
        case mesh::CellType::TET:
            shape = new Tet();
            break;
        case mesh::CellType::HEX:
            shape = new Hex();
            break;
        }
        if (shape == nullptr)
        {
            error::InvalidCellError(cell.id, (int)cell.type, cell.order, __FILE__, __LINE__);
        }
        return shape;
    }
}