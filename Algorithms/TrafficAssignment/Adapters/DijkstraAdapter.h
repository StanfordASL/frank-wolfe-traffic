#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>
#include "DataStructures/Graph/Graph.h"
#include <lemon/dijkstra.h>
#include <lemon/static_graph.h>

using namespace lemon;

// The search algorithm using the graph and possibly auxiliary data to compute shortest paths.
class DijkstraAdapter {
public:
	// Constructs a query algorithm instance working on the specified data.
	explicit DijkstraAdapter(Graph& graph) : graph(graph), weightMap(graph.getWeights(), lemonGraph),
											 dijkstra(lemonGraph, weightMap){ }

	// Computes shortest paths from each source to its target simultaneously.
	void run(const int source, const int target, std::list<int>& path) {
		path.clear();

		Node s = lemonGraph.node(source);
		Node t = lemonGraph.node(target);

		if (source != currentSource){
			currentSource = source;

			dijkstra.run(s);
		}
		
		Node node = t;
		while (node != s)
		{
			Arc in_arc =  dijkstra.predMap()[node];
			node = lemonGraph.source(in_arc);
			path.push_front(lemonGraph.id(in_arc));
		}

		assert(!path.empty());		
	}

	void preprocess()
		{
			std::vector<std::pair<int,int>> edges(graph.numEdges());
			for (int e = 0; e < graph.numEdges(); ++e)
				edges[e] = std::make_pair(graph.tail(e), graph.head(e));
		
			lemonGraph.build(graph.numVertices(), edges.begin(), edges.end());
	
			/*
			  // weightMap = WeightMap(graph.getWeights(), lemonGraph);
			// start and target nodes
			Node s = g.node(0);
			Node t = g.node(3);
		
			dij.run(s,t);

			// extract shortest path (in reverse order)
			Node node = t;
			while (node != s)
			{
			Arc in_arc =  dij.predMap()[c_n];
			node = g.source(in_arc);
		
			cout << g.id(node) << endl;
			}*/
		}

	void customize()
		{
			//dijkstra = Dijkstra<LemonGraph, WeightMap>(lemongGraph, weightMap);
			//dijkstra.init();
			//dijkstra = Dijkstra<LemonGraph, WeightMap>(lemonGraph, weightMap);
			dijkstra.lengthMap(weightMap);
			int currentSource = 0;
			dijkstra.run(lemonGraph.node(currentSource));
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
	
	struct WeightMap 
	{
		typedef double Value;
		WeightMap(std::vector<double>& weights, LemonGraph& lg) : weights(weights), lg(lg) { }
		
		double operator[](Arc e) const
			{
				return weights[lg.index(e)];
			}

		std::vector<double>& weights;
		LemonGraph& lg;
	};
	
	Graph& graph;           // The input graph.
	LemonGraph lemonGraph;  // The graph used for dijkstra search
	WeightMap weightMap;	// Specifies edge weights for Dijkstra search
	Dijkstra<LemonGraph, WeightMap> dijkstra; // Dijkstra search
	int currentSource;
};


/*
// LEMON DEMO CODE

// TEST BEGIN
std::vector<std::pair<int,int>> edges;
edges.push_back(std::make_pair(0,1));
edges.push_back(std::make_pair(0,2));
edges.push_back(std::make_pair(1,2));	
edges.push_back(std::make_pair(1,3));
edges.push_back(std::make_pair(2,3));
edges.push_back(std::make_pair(3,0));

lemonGraph.clear();
		
lemonGraph.build(4, edges.begin(), edges.end());
//WeightMap weightMap(graph.getWeights());
ArcMap<double> length(lemonGraph, 1.0);
length[lemonGraph.arc(1)] = 0.1;
length[lemonGraph.arc(5)] = 8.0;
length[lemonGraph.arc(2)] = 5.0;
length[lemonGraph.arc(4)] = 0.1;

Dijkstra<LemonGraph, ArcMap<double>> dijkstra(lemonGraph, length);

// start and target nodes
Node s = lemonGraph.node(0);
Node t = lemonGraph.node(3);
		
dijkstra.run(s,t);

// extract shortest path (in reverse order)
Node node = t;
while (node != s)
{
Arc in_arc =  dijkstra.predMap()[node];
node = lemonGraph.source(in_arc);
		
std::cout << lemonGraph.id(node) << std::endl;
}
*/
