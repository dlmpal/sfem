#pragma once

#include "config.h"

namespace sfem::mat_ops
{
    inline Float det3x3(const Float m[])
    {

        Float det = m[0] * (m[4] * m[8] - m[7] * m[5]) -
                    m[1] * (m[3] * m[8] - m[6] * m[5]) +
                    m[2] * (m[3] * m[7] - m[4] * m[6]);

        return det;
    }

    inline Float det2x2(const Float m[])
    {

        Float det = m[0] * m[3] - m[1] * m[2];

        return det;
    }

    inline double inv3x3(const Float m[], Float mi[])
    {

        Float di = 1 / det3x3(m);

        mi[0] = di * (m[4] * m[8] - m[5] * m[7]);
        mi[1] = -di * (m[1] * m[8] - m[2] * m[7]);
        mi[2] = di * (m[1] * m[5] - m[2] * m[4]);

        mi[3] = -di * (m[3] * m[8] - m[5] * m[6]);
        mi[4] = di * (m[0] * m[8] - m[2] * m[6]);
        mi[5] = -di * (m[0] * m[5] - m[2] * m[3]);

        mi[6] = di * (m[3] * m[7] - m[4] * m[6]);
        mi[7] = -di * (m[0] * m[7] - m[1] * m[6]);
        mi[8] = di * (m[0] * m[4] - m[1] * m[3]);

        return 1.0 / di;
    }

    inline double inv2x2(const Float m[], Float mi[])
    {

        Float di = 1 / det2x2(m);

        mi[0] = di * m[3];
        mi[1] = -di * m[1];
        mi[2] = -di * m[2];
        mi[3] = di * m[0];

        return 1.0 / di;
    }

    inline double inv(const int r, const Float m[][3], Float mi[][3])
    {
        Float _m[3 * 3];
        Float _mi[3 * 3];
        Float det;

        for (auto i = 0; i < r; i++)
            for (auto j = 0; j < r; j++)
                _m[i * r + j] = m[i][j];

        if (r == 3)
            det = inv3x3(_m, _mi);
        else if (r == 2)
            det = inv2x2(_m, _mi);
        else
        {
            _mi[0] = 1 / _m[0];
            det = _mi[0];
        }

        for (auto i = 0; i < r; i++)
            for (auto j = 0; j < r; j++)
                mi[i][j] = _mi[i * r + j];

        return det;
    }

    inline void mat_set(const int size, const Float val, Float m[])
    {
        for (int i = 0; i < size; i++)
        {
            m[i] = val;
        }
    }

    inline void mat_copy(const int size, const Float source[], Float dest[])
    {
        for (int i = 0; i < size; i++)
        {
            dest[i] = source[i];
        }
    }

    inline void mat_mult_scalar(const int size, const Float val, Float m[])
    {
        for (int i = 0; i < size; i++)
        {
            m[i] *= val;
        }
    }

    inline bool mat_check_equal(const int size, const Float tol, const Float m1[], const Float m2[])
    {
        bool eq = true;
        for (int i = 0; i < size; i++)
        {
            eq = (fabs(m1[i] - m2[i]) < tol);
        }
        return eq;
    }

    inline void mat_print(const int r, const int c, const Float m[])
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                std::cout << m[i * c + j] << " ";
            }
            std::cout << "\n";
        }
    }

    inline void mat_transpose(const int r, const int c, const Float m[], Float mt[])
    {
        for (auto i = 0; i < r; i++)

            for (auto j = 0; j < c; j++)

                mt[j * r + i] = m[i * c + j];
    }

    inline void mat_mult(const int r, const int c, const int rc, const Float m1[], const Float m2[], Float m[])
    {

        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                m[i * c + j] = 0.0;

                for (int k = 0; k < rc; k++)
                {
                    m[i * c + j] += m1[i * rc + k] * m2[k * c + j];
                }
            }
        }
    }

    inline void mat_add(const int size, bool inplace, const Float m1[], Float m2[], Float m3[])
    {
        if (inplace)
        {
            for (int i = 0; i < size; i++)
            {
                m2[i] += m1[i];
            }
        }
        else
        {
            for (int i = 0; i < size; i++)
            {
                m3[i] = m1[i] + m2[i];
            }
        }
    }

    inline Float mat_norm(const int size, const Float m[])
    {
        Float norm = 0.0;
        for (int i = 0; i < size; i++)
        {

            norm += m[i] * m[i];
        }
        return sqrt(norm);
    }
}
