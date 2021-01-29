#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include <vectorclass/vectorclass.h>

#include "Algorithms/TrafficAssignment/AllOrNothingAssignment.h"
#include "Algorithms/TrafficAssignment/UnivariateMinimization.h"
#include "DataStructures/Graph/Attributes/TravelCostAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Stats/TrafficAssignment/FrankWolfeAssignmentStats.h"
#include "Tools/Simd/AlignedVector.h"
#include "Tools/Timer.h"

// A traffic assignment procedure based on the Frank-Wolfe method (also known as convex combinations
// method). At its heart are iterative shortest-paths computations. The algo can be parameterized to
// compute the user equilibrium or system optimum, and to use different travel cost functions and
// shortest-path algorithms.
template <
    template <typename> class ObjFunctionT, template <typename> class TravelCostFunctionT,
    template <typename, typename> class ShortestPathAlgoT, typename InputGraphT>
class FrankWolfeAssignment {
public:
	using InputGraph = InputGraphT;

	// Constructs an assignment procedure based on the Frank-Wolfe method.
	FrankWolfeAssignment(InputGraphT& graph, const std::vector<ClusteredOriginDestination>& odPairs, std::ofstream& csv, std::ofstream& distFile, std::ofstream& patternFile, std::ofstream& pathFile, std::ofstream& weightFile, const bool verbose = true)
		: allOrNothingAssignment(graph, odPairs, verbose),
		  graph(graph),	
		  trafficFlows(graph.numEdges()),
		  pointOfSight(graph.numEdges()),
		  travelCostFunction(graph),
		  objFunction(travelCostFunction),
		  csv(csv),
		  distanceFile(distFile),
		  patternFile(patternFile),
		  pathFile(pathFile),
		  weightFile(weightFile),
		  verbose(verbose) {
		stats.totalRunningTime = allOrNothingAssignment.stats.totalRoutingTime;
	}

	// Assigns all OD-flows onto the input graph.
	void run(const int numIterations = 0, const std::vector<int>& samplingIntervals = {}) {
		assert(numIterations >= 0);
		assert(samplingIntervals.empty() || samplingIntervals[0] > 0);
		for (int i = 1; i < samplingIntervals.size(); ++i) {
			assert(samplingIntervals[i] > 0);
			assert(samplingIntervals[i - 1] % samplingIntervals[i] == 0);
		}
		const AllOrNothingAssignmentStats& substats = allOrNothingAssignment.stats;

		std::vector<double> weights = std::vector<double>(numIterations, 0.0);
		weights[0]=1.0;
		
		Timer timer;
		const int interval = !samplingIntervals.empty() ? samplingIntervals[0] : 1;
		determineInitialSolution(interval);
		paths = allOrNothingAssignment.getPaths();

		stats.lastRunningTime = timer.elapsed();
		stats.lastLineSearchTime = stats.lastRunningTime - substats.lastRoutingTime;
		stats.objFunctionValue = objFunction(trafficFlows);
		stats.finishIteration();

		if (csv.is_open()) {
			csv << substats.numIterations << "," << interval << ",";
			csv << substats.lastCustomizationTime << "," << substats.lastQueryTime << ",";
			csv << stats.lastLineSearchTime << "," << stats.lastRunningTime << ",nan,nan,";
			csv << stats.objFunctionValue << "," << stats.totalTravelCost << ",";
			csv << substats.lastChecksum << std::endl;
		}

		if (distanceFile.is_open())
			for (const auto dist : substats.lastDistances)
				distanceFile << substats.numIterations << ',' << dist << '\n';

		if (pathFile.is_open())
			{
				for (auto i = 0; i < paths.size(); i++)
				{
					pathFile << substats.numIterations << ',' << i;
					for(const auto& e : paths[i])
						pathFile << "," << e;

					pathFile << '\n';
				}
			}
		

		if (verbose) {
			std::cout << "  Line search: " << stats.lastLineSearchTime << "ms";
			std::cout << "  Total: " << stats.lastRunningTime << "ms\n";
			std::cout << "  Objective function value: " << stats.objFunctionValue << "\n";
			std::cout << "  Total travel cost: " << stats.totalTravelCost << "\n";
			std::cout << std::flush;
		}

		// Perform iterations of Frank-Wolfe		
		do {
			stats.startIteration();
			Timer timer;

			// Update travel costs
			updateTravelCosts();

			// Direction finding.
			const int interval = substats.numIterations < samplingIntervals.size() ?  samplingIntervals[substats.numIterations] : 1;
			findDescentDirection(interval);
			paths = allOrNothingAssignment.getPaths();

			const auto tau = findMoveSize();
			moveAlongDescentDirection(tau);

			// update weights vector
			for (auto i = 0; i < substats.numIterations - 1; i++)
				weights[i] = weights[i] * (1.0-tau);

			weights[substats.numIterations-1] = tau;
									
			stats.lastRunningTime = timer.elapsed();
			stats.lastLineSearchTime = stats.lastRunningTime - substats.lastRoutingTime;
			stats.objFunctionValue = objFunction(trafficFlows);
			stats.finishIteration();

			if (csv.is_open()) {
				csv << substats.numIterations << "," << interval << ",";
				csv << substats.lastCustomizationTime << "," << substats.lastQueryTime << ",";
				csv << stats.lastLineSearchTime << "," << stats.lastRunningTime << ",";
				csv << substats.avgChangeInDistances << "," << substats.maxChangeInDistances << ",";
				csv << stats.objFunctionValue << "," << stats.totalTravelCost << ",";
				csv << substats.lastChecksum << std::endl;
			}

			if (distanceFile.is_open())
				for (const auto dist : substats.lastDistances)
					distanceFile << substats.numIterations << ',' << dist << '\n';
			
			if (pathFile.is_open())
			{
				for (auto i = 0; i < paths.size(); i++)
				{
					pathFile << substats.numIterations << ',' << i;
					for(const auto& e : paths[i])
						pathFile << "," << e;

					pathFile << '\n';
				}
			}
			

			if (verbose) {
				std::cout << "  Line search: " << stats.lastLineSearchTime << "ms";
				std::cout << "  Total: " << stats.lastRunningTime << "ms\n";
				std::cout << "  Max change in OD-distances: " << substats.maxChangeInDistances << "\n";
				std::cout << "  Avg change in OD-distances: " << substats.avgChangeInDistances << "\n";
				std::cout << "  Objective function value: " << stats.objFunctionValue << "\n";
				std::cout << "  Total travel cost: " << stats.totalTravelCost << "\n";
				std::cout << std::flush;
			}
		} while ((numIterations > 0 || substats.avgChangeInDistances > 1e-2) &&
				 (numIterations == 0 || substats.numIterations != numIterations));

		if (verbose) {
			std::cout << "Total:\n";
			std::cout << "  Checksum: " << substats.totalChecksum;
			std::cout << "  Prepro: " << substats.totalPreprocessingTime << "ms";
			std::cout << "  Custom: " << substats.totalCustomizationTime << "ms";
			std::cout << "  Queries: " << substats.totalQueryTime << "ms";
			std::cout << "  Routing: " << substats.totalRoutingTime << "ms\n";
			std::cout << "  Line search: " << stats.totalLineSearchTime << "ms";
			std::cout << "  Total: " << stats.totalRunningTime << "ms\n";
			std::cout << std::flush;
		}

		if (patternFile.is_open()) 
			FORALL_EDGES(graph, e)
			{			
				const int tail = graph.edgeTail_z(e);
				const int head = graph.edgeHead(e);
				const auto flow = trafficFlows[e];
			
				patternFile << substats.numIterations << ',' << tail << ',' << head << ',' << graph.travelTime(e) << ',' << travelCostFunction(e, flow) << ',' << graph.capacity(e) << ',' << flow << '\n';
			}

		if (weightFile.is_open())
			for (int i=0; i < numIterations; i++)
				weightFile << i+1 << "," << weights[i] << std::endl;
		
	}

	void determineInitialSolution(const int interval) {
		//int iterations = 7;

		FORALL_EDGES(graph, e)
			graph.travelCost(e) = objFunction.derivative(e, 0);

		allOrNothingAssignment.run(interval);

		/*FORALL_EDGES(graph, e)
			trafficFlows[e] = allOrNothingAssignment.trafficFlowOn(e);

		for (int i = 2; i <= iterations; i++){		
			FORALL_EDGES(graph, e) 
				graph.travelCost(e) = objFunction.derivative(e, trafficFlows[e]);

			findDescentDirection(interval);
				
			const auto tau = findMoveSize();
			moveAlongDescentDirection(tau);		
			}  */

		FORALL_EDGES(graph, e)
		{
			trafficFlows[e] = allOrNothingAssignment.trafficFlowOn(e);
			stats.totalTravelCost += trafficFlows[e] * travelCostFunction(e, trafficFlows[e]);
		}

		// allOrNothingAssignment.stats.numIterations = 1;
	}

	// Updates traversal costs.
	void updateTravelCosts() {
		FORALL_EDGES(graph, e) 
			graph.travelCost(e) = objFunction.derivative(e, trafficFlows[e]);
	}

	// Finds the descent direction.
	void findDescentDirection(const int skipInterval) {
		allOrNothingAssignment.run(skipInterval);

		if (allOrNothingAssignment.stats.numIterations == 2) {
			FORALL_EDGES(graph, e)
				pointOfSight[e] = allOrNothingAssignment.trafficFlowOn(e);
			return;
		}

		auto num = 0.0, den = 0.0;
		FORALL_EDGES(graph, e) {
			const auto residualDirection = pointOfSight[e] - trafficFlows[e];
			const auto secondDerivative = objFunction.secondDerivative(e, trafficFlows[e]);
			const auto fwDirection = allOrNothingAssignment.trafficFlowOn(e) - trafficFlows[e];
			num += residualDirection * secondDerivative * fwDirection;
			den += residualDirection * secondDerivative * (fwDirection - residualDirection);
		}

		const auto alpha = std::min(std::max(0.0, num / den), 1 - 1e-15);
    
		FORALL_EDGES(graph, e)
			pointOfSight[e] = alpha * pointOfSight[e] + (1 - alpha) * allOrNothingAssignment.trafficFlowOn(e);
	}

	// Find the optimal move size.
	double findMoveSize() const {
		return bisectionMethod([this](const double tau) {
								   auto sum = 0.0;
								   FORALL_EDGES(graph, e) {
									   const auto direction = pointOfSight[e] - trafficFlows[e];
									   sum += direction * objFunction.derivative(e, trafficFlows[e] + tau * direction);
								   }
								   return sum;
							   }, 0, 1);
	}

	// Moves along the descent direction.
	void moveAlongDescentDirection(const double tau) {
		FORALL_EDGES(graph, e)
		{
			trafficFlows[e] += tau * (pointOfSight[e] - trafficFlows[e]);
			stats.totalTravelCost += trafficFlows[e] * travelCostFunction(e, trafficFlows[e]);
		}  
	}
	
	// Returns the traffic flow on edge e.
	const double& trafficFlowOn(const int e) const {
		assert(e >= 0); assert(e < graph.numEdges());
		return trafficFlows[e];
	}

	FrankWolfeAssignmentStats stats; // Statistics about the execution.

private:
	using AllOrNothing = AllOrNothingAssignment<ShortestPathAlgoT<InputGraphT, TravelCostAttribute>>;
	using TravelCostFunction = TravelCostFunctionT<InputGraphT>;
	using ObjFunction = ObjFunctionT<TravelCostFunction>;

	AllOrNothing allOrNothingAssignment;   // The all-or-nothing assignment algo used as a subroutine.
	InputGraphT& graph;               // The input graph.
	// InputGraphT graphReversed;        // Reversed input graph (needed to get edge tails)
	AlignedVector<double> trafficFlows;    // The traffic flows on the edges.
	std::vector<double> pointOfSight;            // The point defining the descent direction d = s - x
	TravelCostFunction travelCostFunction; // A functor returning the travel cost on an edge.
	ObjFunction objFunction;               // The objective function to be minimized (UE or SO).
	std::ofstream& csv;                    // The output CSV file containing statistics.
	std::ofstream& distanceFile;           // The output file containing the OD-distances.
	std::ofstream& patternFile;            // The output file containing the flow patterns.
	std::ofstream& pathFile;				// Output file for individual paths
	std::ofstream& weightFile;				// Output file for path weights
	const bool verbose;                    // Should informative messages be displayed?
	std::vector<std::list<int>> paths;	// paths of the individual od pairs
};

// An alias template for a user-equilibrium (UE) traffic assignment.
template <
    template <typename> class TravelCostFunctionT,
    template <typename, typename> class ShortestPathAlgoT, typename InputGraphT>
using UEAssignment =
								 FrankWolfeAssignment<UserEquilibrium, TravelCostFunctionT, ShortestPathAlgoT, InputGraphT>;

// An alias template for a system-optimum (SO) traffic assignment.
template <
    template <typename> class TravelCostFunctionT,
    template <typename, typename> class ShortestPathAlgoT, typename InputGraphT>
using SOAssignment =
								 FrankWolfeAssignment<SystemOptimum, TravelCostFunctionT, ShortestPathAlgoT, InputGraphT>;
