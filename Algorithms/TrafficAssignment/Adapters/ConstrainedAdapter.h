#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>
#include "DataStructures/Graph/Graph.h"
#include <lemon/dijkstra.h>
#include <lemon/static_graph.h>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/r_c_shortest_paths.hpp>
#include <iostream>
#include <map>
#include <limits>

using namespace lemon;
using namespace boost;

struct VertProp 
{
	VertProp(int num=0, double min_dist=0, double max_dist=0) : num(num), min_distance(min_dist), max_distance(max_dist) { }
    int num;				// id
    double min_distance;	// minimal traversal distance
    double max_distance;	// maximal traversal distance
};

struct EdgeProp
{
	EdgeProp(int n, double t, double d) : num(n), time(t), distance(d)  { }
				
    int num; // id
    double time; // traversal time
    double distance; // traversal distance
};

typedef adjacency_list<vecS, vecS, directedS, VertProp, EdgeProp> BoostGraph;
typedef graph_traits<BoostGraph>::edge_descriptor Edge;
typedef graph_traits<BoostGraph>::vertex_descriptor Vertex;

// data structures for shortest path problem with distance windows (spptw)
// ResourceContainer model 
struct ResourceContainer
{
	ResourceContainer(int time = 0, int distance = 0) : time(time), distance(distance) {}
	
    ResourceContainer& operator=(const ResourceContainer& other)
    {
        if (this == &other)
            return *this;
		
        this->time = other.time;
		this->distance = other.distance;

		return *this;
    }
		double time;
    double distance;
};

bool operator==(const ResourceContainer& res_cont_1, const ResourceContainer& res_cont_2)
{
    return (res_cont_1.time == res_cont_2.time && res_cont_1.distance == res_cont_2.distance);
}

bool operator<(const ResourceContainer& res_cont_1, const ResourceContainer& res_cont_2)
{
    if (res_cont_1.time > res_cont_2.time)
        return false;
    if (res_cont_1.time == res_cont_2.time)
        return res_cont_1.distance < res_cont_2.distance;
    return true;
}

// ResourceExtensionFunction model
class ResourceExtensionFunction
{
public:
inline bool operator()(const BoostGraph& g, ResourceContainer& new_cont, const ResourceContainer& old_cont, Edge ed) const
{
	const EdgeProp& arc_prop = get(edge_bundle, g)[ed];
	const VertProp& vert_prop = get(vertex_bundle, g)[target(ed, g)];
	new_cont.time = old_cont.time + arc_prop.time;
	double& i_distance = new_cont.distance;
	i_distance = old_cont.distance + arc_prop.distance;
	i_distance < vert_prop.min_distance ? i_distance = vert_prop.min_distance : 0;
	return i_distance <= vert_prop.max_distance ? true : false;
}
};

// DominanceFunction model
class DominanceFunction
{
public:
    inline bool operator()(const ResourceContainer& res_cont_1,
						   const ResourceContainer& res_cont_2) const
    {
        // must be "<=" here, not "<"
        return res_cont_1.time <= res_cont_2.time
            && res_cont_1.distance <= res_cont_2.distance;
    }
};
// end data structures for shortest path problem with distance windows (spptw)

// The search algorithm using the graph and possibly auxiliary data to compute shortest paths.
class ConstrainedAdapter {
public:
	// Constructs a query algorithm instance working on the specified data.
	explicit ConstrainedAdapter(Graph& graph) : graph(graph), weights(graph.getWeights()), lengthMap(edgeLengths, lemonGraph), dijkstra(lemonGraph, lengthMap), normalDistanceMultiplier(graph.noramlDistanceMultiplier()) { }
	
	// Computes shortest paths from source to target
	void run(const int source_id, const int target_id, std::list<int>& path) {
		path.clear();

		double st_distance;
		std::pair<int,int> st_pair = std::make_pair(source_id, target_id);

		// find normal distance between source and target
		if (distances.find(st_pair) != distances.end())
			st_distance = distances[st_pair];
		else
		{
			Node s = lemonGraph.node(source_id);
			Node t = lemonGraph.node(target_id);
			
			dijkstra.run(s,t);
			
			st_distance = dijkstra.dist(t);		
			distances[st_pair] = st_distance;
		}

		// set distance window for the target vertex
		boostGraph[target_id].max_distance = normalDistanceMultiplier * st_distance;

		Vertex s = boost_vertices[source_id];
		Vertex t = boost_vertices[target_id];

		std::vector<std::vector<Edge>> opt_solutions;
		std::vector<ResourceContainer> pareto_opt_rcs;

		r_c_shortest_paths(boostGraph, get(&VertProp::num, boostGraph),
						   get(&EdgeProp::num, boostGraph), s, t, opt_solutions,
						   pareto_opt_rcs, ResourceContainer(0, 0), ResourceExtensionFunction(),
						   DominanceFunction(),
						   std::allocator<r_c_shortest_paths_label<BoostGraph, ResourceContainer>>(),
						   default_r_c_shortest_paths_visitor());

		// find the shortest (time-wise) path of all the pareto optimal paths
		double min_time = pareto_opt_rcs[0].time;
		int min_index = 0;

		for (int i = 1; i < opt_solutions.size(); ++i)
		{
			if (pareto_opt_rcs[i].time <  min_time)
			{
				min_time = pareto_opt_rcs[i].time;
				min_index = i;
			}
		}
		
		// extract path
		for (int j = opt_solutions[min_index].size() - 1; j >= 0; --j)
			path.push_back(edge_map[std::make_pair(source(opt_solutions[min_index][j], boostGraph), target(opt_solutions[min_index][j], boostGraph))]);

		// revert distance window for the target vertex
		boostGraph[target_id].max_distance = std::numeric_limits<double>::max();
	}
	
	void preprocess(){
		// prepare distances map
		std::vector<std::pair<int,int>> edges(graph.numEdges());
		edgeLengths = std::vector<double>(graph.numEdges());
		
		for (int e = 0; e < graph.numEdges(); ++e){
			edges[e] = std::make_pair(graph.tail(e), graph.head(e));
			edgeLengths[e] = graph.length(e);
		}
		
		lemonGraph.build(graph.numVertices(), edges.begin(), edges.end());
		dijkstra = Dijkstra<LemonGraph, LengthMap>(lemonGraph, lengthMap);
	}

	void customize(){
		// construct the boost graph
		boostGraph.clear();
		
		boost_vertices.clear();
		for (int v = 0; v < graph.numVertices(); v++)
			boost_vertices.push_back(add_vertex(VertProp(v, 0.0, std::numeric_limits<double>::max()), boostGraph));
		
		for (int e = 0; e < graph.numEdges(); e++)
		{
			add_edge(graph.tail(e), graph.head(e), EdgeProp(e, weights[e], graph.length(e)), boostGraph);
			edge_map[std::make_pair(graph.tail(e), graph.head(e))] = e;
		}
	}
	
private:
	using LemonGraph = StaticDigraph;
	using Node = LemonGraph::Node;
	using Arc = LemonGraph::Arc;
	using NodeIt = LemonGraph::NodeIt;
	using ArcIt = LemonGraph::ArcIt;

	template<typename Type>
	using NodeMap = LemonGraph::NodeMap<Type>;

	template<typename Type>
	using ArcMap = LemonGraph::ArcMap<Type>;
	
	struct LengthMap 
	{
		typedef double Value;
		LengthMap(std::vector<double>& lengths, LemonGraph& lg) : lengths(lengths), lg(lg) { }
		
		double operator[](Arc e) const
		{
			return lengths[lg.index(e)];
		}

		std::vector<double>& lengths;
		LemonGraph& lg;
	};

	Graph& graph;									// Input graph
	std::vector<double>& weights;					// Specifies edge travel time for search
	std::map<std::pair<int,int>,double> distances;	// Specifies normal distances used as constraints
	BoostGraph boostGraph;							// Graph for constrained search
	LemonGraph lemonGraph;							// Graph for computing the normal distance
	LengthMap lengthMap;
	Dijkstra<LemonGraph, LengthMap> dijkstra;		// Dijkstra search for normal distance
	std::map<std::pair<int,int>,int> edge_map;		// Maps from pair edges to edge index
	double normalDistanceMultiplier;
	std::vector<double> boost_vertices;
	std::vector<double> edgeLengths;
	
};
