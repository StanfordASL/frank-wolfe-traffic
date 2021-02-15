#pragma once

#include <list>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <vector>

#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Stats/TrafficAssignment/AllOrNothingAssignmentStats.h"
#include "Tools/Simd/AlignedVector.h"
#include "Tools/Timer.h"

// Implementation of an iterative all-or-nothing traffic assignment. Each OD-pair is processed in
// turn and the corresponding OD-flow (in our case always a single flow unit) is assigned to each
// edge on the shortest path between O and D. Other O-D paths are not assigned any flow. The
// procedure can be used with different shortest-path algorithms.
template <typename ShortestPathAlgoT>
class AllOrNothingAssignment {
public:
	// Constructs an all-or-nothing assignment instance.
	AllOrNothingAssignment(Graph& graph,
						   const std::vector<ClusteredOriginDestination>& odPairs,
						   const bool verbose = true)
		: stats(odPairs.size()),
		  shortestPathAlgo(graph),
		  inputGraph(graph),
		  odPairs(odPairs),
		  verbose(verbose)
		{
			Timer timer;
			shortestPathAlgo.preprocess();
			stats.totalPreprocessingTime = timer.elapsed();
			stats.lastRoutingTime = stats.totalPreprocessingTime;
			stats.totalRoutingTime = stats.totalPreprocessingTime;
			if (verbose) std::cout << "  Prepro: " << stats.totalPreprocessingTime << "ms" << std::endl;
			paths = std::vector<std::list<int>>(odPairs.size(), std::list<int>());
		
		}

	// Assigns all OD-flows to their currently shortest paths.
	void run() {
		Timer timer;
		++stats.numIterations;
		if (verbose) std::cout << "Iteration " << stats.numIterations << ": " << std::flush;
	
		shortestPathAlgo.customize();
		stats.lastCustomizationTime = timer.elapsed();

		timer.restart();

		// assign initial flow of 0 to each edge
		trafficFlows.assign(inputGraph.numEdges(), 0);
		stats.startIteration();

		// find shortest path between each OD pair and collect flows
		for (int i = 0; i < odPairs.size(); i++)
		{
			shortestPathAlgo.run(odPairs[i].origin, odPairs[i].destination, paths[i]);
			for(const auto& e : paths[i])
				trafficFlows[e] += odPairs[i].volume;
		}
		
		stats.lastQueryTime = timer.elapsed();
		stats.finishIteration();

		if (verbose) {
			std::cout << " done.\n";
			std::cout << "  Checksum: " << stats.lastChecksum;
			std::cout << "  Custom: " << stats.lastCustomizationTime << "ms";
			std::cout << "  Queries: " << stats.lastQueryTime << "ms";
			std::cout << "  Routing: " << stats.lastRoutingTime << "ms\n";
			std::cout << std::flush;
		}
	}

	// Returns the traffic flow on edge e.
	double trafficFlowOn(const int e) const {
		assert(e >= 0); assert(e < inputGraph.numEdges());
		return (double)trafficFlows[e];
	}

	std::vector<std::list<int>>& getPaths()
	{
		return paths;
	}

	AllOrNothingAssignmentStats stats; // Statistics about the execution.

private:
	using ODPairs = std::vector<ClusteredOriginDestination>;

	ShortestPathAlgoT shortestPathAlgo; // Algo computing shortest paths between OD-pairs.
	Graph& inputGraph;					// The input graph.
	const ODPairs& odPairs;             // The OD-pairs to be assigned onto the graph.
	std::vector<int> trafficFlows;			// The traffic flows on the edges.
	std::vector<std::list<int>> paths;	// paths of the individual od pairs
	const bool verbose;                 // Should informative messages be displayed?
	
};
