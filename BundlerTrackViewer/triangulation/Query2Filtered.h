/*
 *  Query2Filtered.h
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
// File Version: 5.0.1 (2011/07/23)

#ifndef _QUERY2FILTERED_H
#define _QUERY2FILTERED_H

#include "Query2Rational.h"

class Query2Filtered : public Query2
{
	public:
		// The base class handles floating-point queries.  Each query involves
		// comparing a determinant to zero.  If the determinant is sufficiently
		// close to zero, numerical round-off errors may cause the determinant
		// sign to be misclassified.  To avoid this, the query is repeated with
		// exact rational arithmetic.  You specify the closeness to zero for the
		// switch to rational arithmetic via 'uncertainty', a value in the
		// interval [0,1].  The uncertainty of 0 causes the class to behave
		// as if it were Query2.  The uncertainty of 1 causes the class to
		// behave as if it were Query2Rational.
		Query2Filtered (int numVertices, const Vec2f* vertices,
						float uncertainty);
		virtual ~Query2Filtered ();
		
		// Run-time type information.
		virtual Query::Type GetType () const;
		
		// Queries about the relation of a point to various geometric objects.
		
		virtual int ToLine (const Vec2f& test, int v0, int v1) const;
		
		virtual int ToCircumcircle (const Vec2f& test, int v0, int v1,
									int v2) const;
		
	private:
		using Query2::mVertices;
		using Query2::Sort;
		using Query2::Det2;
		using Query2::Det3;
		
		//Query2Rational mRQuery;
		float mUncertainty;
};
	
#include "Query2Filtered.inl"

#endif