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
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Tools/Timer.h"
#include "Stats/TrafficAssignment/FrankWolfeAssignmentStats.h"

// #define TA_NO_CFW 

// A traffic assignment procedure based on the Frank-Wolfe method (also known as convex combinations
// method). At its heart are iterative shortest-paths computations. The algo can be parameterized to
// compute the user equilibrium or system optimum, and to use different travel cost functions and
// shortest-path algorithms.
template <
    template <typename> class ObjFunctionT, typename TravelCostFunction,
    typename ShortestPathAlgoT>
class FrankWolfeAssignment {
public:

	// Constructs an assignment procedure based on the Frank-Wolfe method.
	FrankWolfeAssignment(Graph& graph, const std::vector<ClusteredOriginDestination>& odPairs, std::ofstream& csv, std::ofstream& patternFile, std::ofstream& pathFile, std::ofstream& weightFile, const bool verbose = true)
		: allOrNothingAssignment(graph, odPairs, verbose),
		  graph(graph),	
		  trafficFlows(graph.numEdges()),
		  pointOfSight(graph.numEdges()),
		  travelCostFunction(graph),
		  objFunction(travelCostFunction),
		  csv(csv),
		  patternFile(patternFile),
		  pathFile(pathFile),
		  weightFile(weightFile),
		  verbose(verbose) {
		stats.totalRunningTime = allOrNothingAssignment.stats.totalRoutingTime;
	}

	// Assigns all OD-flows onto the input graph.
	void run(const int numIterations = 0) {
		assert(numIterations >= 0);
		const AllOrNothingAssignmentStats& substats = allOrNothingAssignment.stats;

		std::vector<double> weights = std::vector<double>(numIterations, 0.0);
		weights[0]=1.0;
		
		Timer timer;
		determineInitialSolution();
		paths = allOrNothingAssignment.getPaths();

		stats.lastRunningTime = timer.elapsed();
		stats.lastLineSearchTime = stats.lastRunningTime - substats.lastRoutingTime;
		stats.objFunctionValue = objFunction(trafficFlows);
		stats.finishIteration();

		if (csv.is_open()) {
			csv << substats.numIterations << "," << substats.lastCustomizationTime << "," << substats.lastQueryTime << ",";
			csv << stats.lastLineSearchTime << "," << stats.lastRunningTime << ",";
			csv << stats.objFunctionValue << "," << stats.totalTravelCost << "," << std::endl;
		}
		
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
			Timer timer;
			stats.startIteration();

			// Update travel costs
			updateTravelCosts();

			// Direction finding.
			findDescentDirection();
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
				csv << substats.numIterations << "," << substats.lastCustomizationTime << "," << substats.lastQueryTime << ",";
				csv << stats.lastLineSearchTime << "," << stats.lastRunningTime << ",";
				csv << stats.objFunctionValue << "," << stats.totalTravelCost << "," << std::endl;
			}	
			
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
				const int tail = graph.tail(e);
				const int head = graph.head(e);
				const auto flow = trafficFlows[e];
			
				patternFile << substats.numIterations << ',' << tail << ',' << head << ',' << graph.freeTravelTime(e) << ',' << travelCostFunction(e, flow) << ',' << graph.capacity(e) << ',' << flow << '\n';
			}

		if (weightFile.is_open())
			for (int i=0; i < numIterations; i++)
				weightFile << i+1 << "," << weights[i] << std::endl;
		
	}

	void determineInitialSolution() {
		FORALL_EDGES(graph, e)
			graph.setWeight(e, objFunction.derivative(e, 0));

		allOrNothingAssignment.run();
		
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
			graph.setWeight(e, objFunction.derivative(e, trafficFlows[e]));
	}

	// Finds the descent direction.
	void findDescentDirection() {
		allOrNothingAssignment.run();
#ifndef TA_NO_CFW
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
#endif
	}

	// Find the optimal move size.
	double findMoveSize() const {
		return bisectionMethod([this](const double tau) {
								   double sum = 0.0;
								   FORALL_EDGES(graph, e) {
#ifndef TA_NO_CFW
									   const auto direction = pointOfSight[e] - trafficFlows[e];
#else
									   const auto direction = allOrNothingAssignment.trafficFlowOn(e) - trafficFlows[e];
#endif
									   sum += direction * objFunction.derivative(e, trafficFlows[e] + tau * direction);
								   }
								   return sum;
							   }, 0, 1);
	}

	// Moves along the descent direction.
	void moveAlongDescentDirection(const double tau) {
		FORALL_EDGES(graph, e)
		{
#ifndef TA_NO_CFW			
			trafficFlows[e] += tau * (pointOfSight[e] - trafficFlows[e]);
#else
			trafficFlows[e] += tau * (allOrNothingAssignment.trafficFlowOn(e) - trafficFlows[e]);
#endif			
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
	using AllOrNothing = AllOrNothingAssignment<ShortestPathAlgoT>;
	using ObjFunction = ObjFunctionT<TravelCostFunction>;

	AllOrNothing allOrNothingAssignment;   // The all-or-nothing assignment algo used as a subroutine.
	Graph& graph;               // The input graph.
	std::vector<double> trafficFlows;    // The traffic flows on the edges.
	std::vector<double> pointOfSight;            // The point defining the descent direction d = s - x
	TravelCostFunction travelCostFunction; // A functor returning the travel cost on an edge.
	ObjFunction objFunction;               // The objective function to be minimized (UE or SO).
	std::ofstream& csv;                    // The output CSV file containing statistics.
	std::ofstream& patternFile;            // The output file containing the flow patterns.
	std::ofstream& pathFile;				// Output file for individual paths
	std::ofstream& weightFile;				// Output file for path weights
	const bool verbose;                    // Should informative messages be displayed?
	std::vector<std::list<int>> paths;	// paths of the individual od pairs
};

// An alias template for a user-equilibrium (UE) traffic assignment.
template <
    typename TravelCostFunction,
    typename ShortestPathAlgoT>
using UEAssignment = FrankWolfeAssignment<UserEquilibrium, TravelCostFunction, ShortestPathAlgoT>;

// An alias template for a system-optimum (SO) traffic assignment.
template <
    typename TravelCostFunction,
    typename ShortestPathAlgoT>
using SOAssignment = FrankWolfeAssignment<SystemOptimum, TravelCostFunction, ShortestPathAlgoT>;
