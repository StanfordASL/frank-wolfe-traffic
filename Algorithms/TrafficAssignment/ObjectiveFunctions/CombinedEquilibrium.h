#pragma once

#include "DataStructures/Graph/Graph.h"
#include "Algorithms/TrafficAssignment/ObjectiveFunctions/SystemOptimum.h"
#include "Algorithms/TrafficAssignment/ObjectiveFunctions/UserEquilibrium.h"

// Represents the combined-equilibrium (CE) objective function. The flow pattern that minimizes the CE
// objective function (while satisfying the flow conservation constraint) that balances between the user-equlibrium and system-optimum objective.
template <typename TravelCostFunctionT>
class CombinedEquilibrium {
public:
	// Constructs an UE objective function.
CombinedEquilibrium(TravelCostFunctionT travelCostFunction, Graph& graph) : travelCostFunction(travelCostFunction), graph(graph), alpha(graph.combinedEquilibriumParameter()), systemOptimumObj(travelCostFunction, graph), userEquilibriumObj(travelCostFunction, graph) {
	}

	// Returns the value of the objective function for the specified edge flows.
	double operator()(const std::vector<double>& flows) const {
		return interpolate(systemOptimumObj(flows), userEquilibriumObj(flows));
	}

	// Returns the weight of edge e, given the flow x on e.
	double derivative(const int e, const double x) const {
		return interpolate(systemOptimumObj.derivative(e, x), userEquilibriumObj.derivative(e, x));
	}

	// Returns the weight of edge e, given the flow x on e.
	double secondDerivative(const int e, const double x) const {
		return interpolate(systemOptimumObj.secondDerivative(e, x), userEquilibriumObj.secondDerivative(e, x));
	}

private:
	double interpolate(const double so_value, const double ue_value) const
	{
		return alpha * so_value + (1 - alpha) * ue_value;
	}

	TravelCostFunctionT travelCostFunction; // A functor returning the travel cost on an edge.
	Graph& graph;
	
	double alpha;

	SystemOptimum<TravelCostFunctionT> systemOptimumObj;
	UserEquilibrium<TravelCostFunctionT> userEquilibriumObj;
};
