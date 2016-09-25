/*
 *  GraphWrapper.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 2/14/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "GraphWrapper.h"
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

using namespace boost;
typedef boost::minstd_rand base_generator_type;

GraphWrapper::GraphWrapper(int noVertices_)
{
	noVertices = noVertices_;
	graph = Graph(noVertices);	
}

void GraphWrapper::addEdge(int nodeIndex1, int nodeIndex2, float weight)
{
	add_edge(nodeIndex1, nodeIndex2, EdgeWeight(weight), graph);
}

bool GraphWrapper::doesEdgeExist(int nodeIndex1, int nodeIndex2, edge_descriptor &ed)
{
	std::pair<edge_descriptor, bool> ep = edge(vertex(nodeIndex1, graph), vertex(nodeIndex2,graph), graph);
	if(ep.second)
	{
		ed = ep.first;
		return true;
	}
	else
	{
		ep = edge(vertex(nodeIndex2, graph), vertex(nodeIndex1,graph), graph);
		if(ep.second)
		{
			ed = ep.first;
			return true;
		}
		else
			return false;
	}
}

void GraphWrapper::uniformlySampleSpanningTree()
{
	base_generator_type generator(42u);
	static_property_map<double> s(0.0);
	
	IndexMap index_map = get(vertex_index, graph);
    predecessors_t predecessors (num_vertices (graph), graph_traits <Graph>::null_vertex ());
    predecessor_map_t predecessor_map (predecessors.begin (), index_map);
	
    /// Initialize predecessor map
	vertex_iter vertex_iter, vertex_end;
    for (tie (vertex_iter, vertex_end) = vertices (graph); vertex_iter != vertex_end; ++vertex_iter)
    {
		put (predecessor_map, *vertex_iter, *vertex_iter);
    }
	
	std::vector<default_color_type> colors(num_vertices(graph)); 
	make_iterator_property_map(&colors[0], get(vertex_index, graph));
	
	random_spanning_tree(graph, generator, *vertices(graph).first, predecessor_map, s, &colors[0]);
}

void GraphWrapper::findAll3CycleInferences(vector<CycleInference> &inferences)
{
	for(int i=0; i<noVertices; i++)
	{
		for(int j=i+1; j<noVertices; j++)
		{
			CycleInference ci;
			
			//form image pair i-j
			Edge e;
			e.first = i;
			e.second = j;
			ci.mainEdge = e;
			
			for(int k=0; k<noVertices; k++)
			{
				if(k==i || k==j)
					continue;
				//check if i-k and k-j exists
				edge_descriptor ed1, ed2;
				if(doesEdgeExist(i, k, ed1) && doesEdgeExist(k, j, ed2))
				{
					vector<Edge> path;
					Edge firstEdge, secondEdge;
					firstEdge.first = i; firstEdge.second = k;
					secondEdge.first = k; secondEdge.second = j;
					path.push_back(firstEdge);
					path.push_back(secondEdge);
					ci.paths.push_back(path);
				}
			}
			
			if(ci.paths.size() > 0)
				inferences.push_back(ci);
		}
	}
}

void GraphWrapper::findAll3Cycles(vector<Cycle> &cycles)
{
	for(int i=0; i<noVertices; i++)
	{
		for(int j=i+1; j<noVertices; j++)
		{
			edge_descriptor ed1;
			if(doesEdgeExist(i, j, ed1))
			{
			
				CycleInference ci;
			
				//form image pair i-j
				Edge e;
				e.first = i;
				e.second = j;
				ci.mainEdge = e;
			
				for(int k=0; k<noVertices; k++)
				{
					if(k==i || k==j)
						continue;
					//check if i-k and k-j exists
					edge_descriptor ed1, ed2;
					if(doesEdgeExist(i, k, ed1) && doesEdgeExist(k, j, ed2) && doesEdgeExist(j, i, ed2))
					{
						bool add=true;
						for(int m=0; m<cycles.size(); m++)
						{
							if(find(cycles[m].path.begin(), cycles[m].path.end(), i) != cycles[m].path.end() &&
							   find(cycles[m].path.begin(), cycles[m].path.end(), k) != cycles[m].path.end() &&
							   find(cycles[m].path.begin(), cycles[m].path.end(), j) != cycles[m].path.end())
							{
								add=false;
								break;
							}
						}
					
						if(add)
						{
							Cycle c;
							c.path.push_back(i);
							c.path.push_back(k);
							c.path.push_back(j);
							cycles.push_back(c);
						}
					}
				}
			}
		}
	}
}

void GraphWrapper::findMinSpanningTree(vector<int> &sourceVertices, vector<int> &targetVertices)
{
	property_map < Graph, edge_weight_t >::type weight = get(edge_weight, graph);
	std::vector < edge_descriptor > spanning_tree;
	
	kruskal_minimum_spanning_tree(graph, std::back_inserter(spanning_tree));
	
	for (std::vector < edge_descriptor >::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei) 
	{
		sourceVertices.push_back(source(*ei, graph));
		targetVertices.push_back(target(*ei, graph));
	}
}

void GraphWrapper::printGraph()
{
	std::pair<vertex_iter, vertex_iter> vp;
	IndexMap index = get(vertex_index, graph);
	EdgeWeightMap edgeWeightMap = get(edge_weight, graph);
	
	std::cout << "Vertices:\n";
	for (vp = vertices(graph); vp.first != vp.second; ++vp.first)
		std::cout << index[*vp.first] <<  " ";
	std::cout << std::endl;
	std::cout << "Edges:\n";
	
    edge_iter ei, ei_end;
	
    for (tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei)
        std::cout << "(" << index[source(*ei, graph)] << "," << index[target(*ei, graph)] << ") has weight " << 
				get(edge_weight, graph, *ei) << "\n";
    std::cout << std::endl;
	
	uniformlySampleSpanningTree();
}