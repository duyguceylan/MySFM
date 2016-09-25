/*
 *  Query2.h
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 8/20/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

// Geometric Tools, LLC
// Copyright (c) 1998-2011
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)
#ifndef _QUERY2_H
#define _QUERY2_H

#include "Query.h"

class Query2 : public Query
{
	public:
		// The base class handles floating-point queries.
		Query2 (int numVertices, const Vec2f* vertices);
		virtual ~Query2 ();
		
		// Run-time type information.
		virtual Query::Type GetType () const;
		
		// Member access.
		inline int GetNumVertices () const;
		inline const Vec2f* GetVertices () const;
		
		// Queries about the relation of a point to various geometric objects.
		
		// Returns:
		//   +1, on right of line
		//   -1, on left of line
		//    0, on the line
		virtual int ToLine (int i, int v0, int v1) const;
		virtual int ToLine (const Vec2f& test, int v0, int v1) const;
		
		// Returns:
		//   +1, outside triangle
		//   -1, inside triangle
		//    0, on triangle
		virtual int ToTriangle (int i, int v0, int v1, int v2) const;
		virtual int ToTriangle (const Vec2f& test, int v0, int v1,
								int v2) const;
		
		// Returns:
		//   +1, outside circumcircle of triangle
		//   -1, inside circumcircle of triangle
		//    0, on circumcircle of triangle
		virtual int ToCircumcircle (int i, int v0, int v1, int v2) const;
		virtual int ToCircumcircle (const Vec2f& test, int v0, int v1,
									int v2) const;
		
		// Helper functions.
		static float Dot (float x0, float y0, float x1, float y1);
		static float Det2 (float x0, float y0, float x1, float y1);
		static float Det3 (float x0, float y0, float z0, float x1, float y1, float z1,
						  float x2, float y2, float z2);
		
	protected:
		// Input points.
		int mNumVertices;
		const Vec2f* mVertices;
		
};
	
//#include "Query2.inl"

#endif