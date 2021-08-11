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
						   const bool verbose = true, const bool elasticRebalance = false)
		: stats(odPairs.size()),
		  shortestPathAlgo(graph),
		  inputGraph(graph),
		  odPairs(odPairs),
		  verbose(verbose),
		  elasticRebalance(elasticRebalance)
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
		if (elasticRebalance) // comptue for elastic AMoD
		{
			/*
				TODO:
				The current subroutine doesn't aggregate shortest path queries according to origin vertex, which may lead to large computation times.

				A better solution would be to have four different shortest path queries to each origin-destination pair. That is, the OD-pair file would specify for each (virtual) origin-destination the id of the real od-pair, and a number in {0,..,3} specifying the type of query it represents. 
			 */
			
			for (int i = 0; i < odPairs.size(); i++)
			{
				std::list<int> path_od, path_dr, path_or;
				double cost_od, cost_dr, cost_or = 0;
				
				cost_od = shortestPathAlgo.run(odPairs[i].origin, odPairs[i].destination, path_od); // passenger path from new origin to real destination
				cost_dr = shortestPathAlgo.run(odPairs[i].destination, odPairs[i].rebalancer, path_dr); // path for rebalancer

				path_or.push_back(odPairs[i].edge1);
				path_or.push_back(odPairs[i].edge2);

				cost_or = inputGraph.weight(odPairs[i].edge1) + inputGraph.weight(odPairs[i].edge2);
				
				if (cost_od + cost_dr < cost_or) 
				{ // real path used
					paths[i] = path_od;
					paths[i].insert(paths[i].end(), path_dr.begin(), path_dr.end());
				} else 
				{ // virtual path used
					paths[i] = path_or;
				}
				
				for(const auto& e : paths[i])
					trafficFlows[e] += odPairs[i].volume;
			}
		}
		else // compute for classic traffic assignment
		{
			for (int i = 0; i < odPairs.size(); i++)
			{
				shortestPathAlgo.run(odPairs[i].origin, odPairs[i].destination, paths[i]);
				for(const auto& e : paths[i])
					trafficFlows[e] += odPairs[i].volume;
			}
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
	const bool elasticRebalance;		// if true, compute compute AMoD with elastic demand
	
};
