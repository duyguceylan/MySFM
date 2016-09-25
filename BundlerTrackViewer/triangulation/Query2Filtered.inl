// Geometric Tools, LLC
// Copyright (c) 1998-2011
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

Query2Filtered::Query2Filtered (int numVertices,
    const Vec2f* vertices, float uncertainty) :
    Query2(numVertices, vertices),
    //mRQuery(numVertices, vertices),
    mUncertainty(uncertainty)
{
    if(0.0 <= mUncertainty && mUncertainty <= 1.0)
       printf("Query2Filtered::Invalid uncertainty\n");
}

Query2Filtered::~Query2Filtered ()
{
}

Query::Type Query2Filtered::GetType () const
{
    return Query::QT_FILTERED;
}

int Query2Filtered::ToLine (const Vec2f& test, int v0, int v1)
    const
{
    const Vec2f& vec0 = mVertices[v0];
    const Vec2f& vec1 = mVertices[v1];

    float x0 = test[0] - vec0[0];
    float y0 = test[1] - vec0[1];
    float x1 = vec1[0] - vec0[0];
    float y1 = vec1[1] - vec0[1];

    float len0 = sqrt(x0*x0 + y0*y0);
    float len1 = sqrt(x1*x1 + y1*y1);
    float scaledUncertainty = mUncertainty*len0*len1;

    float det = Det2(x0, y0, x1, y1);
    if (abs(det) >= scaledUncertainty)
    {
        return (det > 0.0 ? +1 : (det < 0.0 ? -1 : 0));
    }

    //return mRQuery.ToLine(test, v0, v1);
}

int Query2Filtered::ToCircumcircle (const Vec2f& test, int v0,
    int v1, int v2) const
{
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

    float len0 = sqrt(d0x*d0x + d0y*d0y + z0*z0);
    float len1 = sqrt(d1x*d1x + d1y*d1y + z1*z1);
    float len2 = sqrt(d2x*d2x + d2y*d2y + z2*z2);
    float scaledUncertainty = mUncertainty*len0*len1*len2;

    float det = Det3(d0x, d0y, z0, d1x, d1y, z1, d2x, d2y, z2);
    if (abs(det) >= scaledUncertainty)
    {
        return (det < 0.0 ? 1 : (det > 0.0 ? -1 : 0));
    }

    //return mRQuery.ToCircumcircle(test, v0, v1, v2);
}
//----------------------------------------------------------------------------