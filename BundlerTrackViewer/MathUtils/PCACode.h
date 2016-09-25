/*
 *  PCA.h
 *  FunctionalModeling
 *
 *  Created by Duygu Ceylan on 9/29/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */


template<class Type, int Dim>
PCA<Type,Dim>::PCA()
{
	clear();
}

/**
 * Reset PCA
 */
template<class Type, int Dim>
void PCA<Type,Dim>::clear()
{
	points.clear();
	weights.clear();
	// Compute covariance matrix
	for ( int i = 0; i < Dim; i++ ){
		centroid[i] = 0;
		for (int j = 0; j < Dim; j++ ){
			C[i][j] = 0;
		}
	}
}

/**
 * Add new point
 */
template<class Type, int Dim>
void PCA<Type,Dim>::addPoint( Point p )
{
	centroid = ( centroid * float( points.size() ) + p ) / float( points.size() + 1 );
	points.push_back( p );
}

/**
 * Add new point
 */
template<class Type, int Dim>
void PCA<Type,Dim>::addPoint( Point p, Type w )
{
	weights.push_back(w);
	points.push_back( p );
}

/**
 * Solve eigen system
 */
template<class Type, int Dim>
void PCA<Type,Dim>::analyze( Point& eigenValues, Matrix& eigenVectors, Point& c )
{
	if( weights.size() > 0)
	{
		if( weights.size() != points.size() )
		{
			cout<< "Don't mix weighted and unweigted addPoint() methods!" << endl;
			return;
		}
		
		// estimate centroid
		float summedWeight = 0;
		for( unsigned i=0;i<points.size();i++ )
		{
			centroid += points[i]*weights[i];
			summedWeight += weights[i];
		}
		if( summedWeight > 1e-20f )
			centroid /= summedWeight;
		
		// compute weighted covariance
		for ( int n = 0; n < points.size(); n++ ){
			for ( int i = 0; i < Dim; i++ ){
				for ( int j = 0; j <= i; j++ ){
					C[i][j] += ( points[n][i]-centroid[i]) * (points[n][j]-centroid[j])*weights[n];
					C[j][i] = C[i][j];
				}
			}
		}
		if( summedWeight > 1e-20f )
			C /= summedWeight;
		
	}
	else
	{
		for (int n = 0; n < points.size(); n++ ){
			for ( int i = 0; i < Dim; i++ ){
				for ( int j = 0; j <= i; j++ ){
					C[i][j] += ( points[n][i]-centroid[i]) * (points[n][j]-centroid[j]);
					C[j][i] = C[i][j];
				}
			}
		}
	}
	
	/// Normalize covariance matrix
	//if ( points.size() > 0 ) C /= float( points.size() );
	
	// Analyze eigenstructure
	int noRows = Dim; 
	MathUtils::eigenValueAnalysis(&(C[0][0]), &noRows, &(eigenValues[0]), &(eigenVectors[0][0]));
	c = centroid;
}
