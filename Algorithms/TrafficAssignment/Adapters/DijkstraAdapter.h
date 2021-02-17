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

		//assert(source != target);

		/*if (source == target)
		  std::cout << "origin and destination are equal (" << source << "\n";*/

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

		//assert(!path.empty()); // Graph not connected!
	}

	void preprocess(){
		std::vector<std::pair<int,int>> edges(graph.numEdges());
		for (int e = 0; e < graph.numEdges(); ++e)
			edges[e] = std::make_pair(graph.tail(e), graph.head(e));
		
		lemonGraph.build(graph.numVertices(), edges.begin(), edges.end());
	}

	void customize(){
		dijkstra.lengthMap(weightMap);
		currentSource = 0;
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
