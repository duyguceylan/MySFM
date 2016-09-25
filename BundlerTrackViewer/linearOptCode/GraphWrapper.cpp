#include "GraphWrapper.h"
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

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
