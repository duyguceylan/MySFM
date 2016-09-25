#include <vector>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#ifndef _GRAPH_WRAPPER_H
#define _GRAPH_WRAPPER_H

using namespace boost;
using namespace std;

typedef property<edge_weight_t, float> EdgeWeight;
typedef std::pair<int, int> Edge;
typedef adjacency_list<listS, vecS, undirectedS, no_property, EdgeWeight > Graph;
typedef property_map<Graph, vertex_index_t>::type IndexMap;
typedef property_map<Graph, edge_weight_t>::type EdgeWeightMap;

typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
typedef graph_traits<Graph>::edge_descriptor edge_descriptor;
typedef graph_traits<Graph>::vertex_iterator vertex_iter;
typedef graph_traits<Graph>::edge_iterator edge_iter;

typedef std::vector <vertex_descriptor> predecessors_t;
typedef iterator_property_map <predecessors_t::iterator, IndexMap, vertex_descriptor, vertex_descriptor&> predecessor_map_t;

struct CycleInference
{
	Edge mainEdge;
	vector<vector<Edge> > paths;
};

struct Cycle
{
	vector<int> path;
};

class GraphWrapper
{
private:
	Graph graph;
	int noVertices;
	
public:
	GraphWrapper(int noVertices_);
	void addEdge(int nodeIndex1, int nodeIndex2, float weight);
	bool doesEdgeExist(int nodeIndex1, int nodeIndex2, edge_descriptor &ed); 
	
	void findAll3Cycles(vector<Cycle> &cycles);
	
	void findMinSpanningTree(vector<int> &sourceVertices, vector<int> &targetVertices);
	
};
