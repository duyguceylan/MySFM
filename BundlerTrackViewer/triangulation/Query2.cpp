// Geometric Tools, LLC
// Copyright (c) 1998-2011
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

#include "Query2.h"

Query2::Query2 (int numVertices, const Vec2f* vertices) :
    mNumVertices(numVertices),
    mVertices(vertices)
{
    if(mNumVertices < 0 || mVertices == NULL)
		printf("Query2::Invalid inputs\n");
}

Query2::~Query2 ()
{
}

Query::Type Query2::GetType () const
{
    return Query::QT_REAL;
}

inline int Query2::GetNumVertices () const
{
    return mNumVertices;
}

inline const Vec2f* Query2::GetVertices () const
{
    return mVertices;
}

int Query2::ToLine (int i, int v0, int v1) const
{
    return ToLine(mVertices[i], v0, v1);
}

int Query2::ToLine (const Vec2f& test, int v0, int v1) const
{
    bool positive = Sort(v0, v1);

    const Vec2f& vec0 = mVertices[v0];
    const Vec2f& vec1 = mVertices[v1];

    float x0 = test[0] - vec0[0];
    float y0 = test[1] - vec0[1];
    float x1 = vec1[0] - vec0[0];
    float y1 = vec1[1] - vec0[1];

    float det = Det2(x0, y0, x1, y1);
    if (!positive)
    {
        det = -det;
    }

    return (det > 0.0 ? +1 : (det < 0.0 ? -1 : 0));
}

int Query2::ToTriangle (int i, int v0, int v1, int v2) const
{
    return ToTriangle(mVertices[i], v0, v1, v2);
}

int Query2::ToTriangle (const Vec2f& test, int v0, int v1, int v2) const
{
    int sign0 = ToLine(test, v1, v2);
    if (sign0 > 0)
    {
        return +1;
    }

    int sign1 = ToLine(test, v0, v2);
    if (sign1 < 0)
    {
        return +1;
    }

    int sign2 = ToLine(test, v0, v1);
    if (sign2 > 0)
    {
        return +1;
    }

    return ((sign0 && sign1 && sign2) ? -1 : 0);
}

int Query2::ToCircumcircle (int i, int v0, int v1, int v2) const
{
    return ToCircumcircle(mVertices[i], v0, v1, v2);
}

int Query2::ToCircumcircle (const Vec2f& test, int v0, int v1, int v2) const
{
    bool positive = Sort(v0, v1, v2);

    const Vec2f& vec0 = mVertices[v0];
    const Vec2f& vec1 = mVertices[v1];
    const Vec2f& vec2 = mVertices[v2];

    float s0x = vec0[0] + test[0];
    float d0x = vec0[0] - test[0];
    float s0y = vec0[1] + test[1];
    float d0y = vec0[1] - test[1];
    float s1x = vec1[0] + test[0];
    float d1x = vec1[0] - test[0];
    float s1y = vec1[1] + test[1];
    float d1y = vec1[1] - test[1];
    float s2x = vec2[0] + test[0];
    float d2x = vec2[0] - test[0];
    float s2y = vec2[1] + test[1];
    float d2y = vec2[1] - test[1];
    float z0 = s0x*d0x + s0y*d0y;
    float z1 = s1x*d1x + s1y*d1y;
    float z2 = s2x*d2x + s2y*d2y;

    float det = Det3(d0x, d0y, z0, d1x, d1y, z1, d2x, d2y, z2);
    if (!positive)
    {
        det = -det;
    }

    return (det < 0.0 ? 1 : (det > 0.0 ? -1 : 0));
}

float Query2::Dot (float x0, float y0, float x1, float y1)
{
    return x0*x1 + y0*y1;
}

float Query2::Det2 (float x0, float y0, float x1, float y1)
{
    return x0*y1 - x1*y0;
}

float Query2::Det3 (float x0, float y0, float z0, float x1, float y1,
    float z1, float x2, float y2, float z2)
{
    float c00 = y1*z2 - y2*z1;
    float c01 = y2*z0 - y0*z2;
    float c02 = y0*z1 - y1*z0;
    return x0*c00 + x1*c01 + x2*c02;
}
//----------------------------------------------------------------------------