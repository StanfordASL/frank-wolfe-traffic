#pragma once

#include "Algorithms/TrafficAssignment/TravelCostFunctions/BprFunction.h"

// The BPR travel cost function, relating the travel time on an edge to the flow on this edge.
class ModifiedBprFunction {
 public:
  // Constructs a BPR function.
	ModifiedBprFunction(const Graph& graph) : graph(graph), bpr(graph) {
	}

  // Returns the travel time on edge e, given the flow x on e.
  double operator()(const int e, const double x) const {
    const double pt = APT * graph.capacity(e); // The point at which we linearize.
    if (x <= pt)
      return bpr(e, x);
    else
      return bpr(e, pt) + bpr.derivative(e, pt) * (x - pt);
  }

  // Returns the derivative of e's travel cost function at x.
  double derivative(const int e, const double x) const {
    const double pt = APT * graph.capacity(e); // The point at which we linearize.
    if (x <= pt)
      return bpr.derivative(e, x);
    else
      return bpr.derivative(e, pt);
  }

	// Returns the derivative of e's travel cost function at x.
  double secondDerivative(const int e, const double x) const {
    const double pt = APT * graph.capacity(e); // The point at which we linearize.
    if (x <= pt)
      return bpr.secondDerivative(e, x);
    else
		return 0;
  }

  // Returns the integral of e's travel cost function from 0 to b.
  double integral(const int e, const double b) const {
    const double pt = APT * graph.capacity(e); // The point at which we linearize.
    if (b <= pt)
      return bpr.integral(e, b);
    else
      return bpr.integral(e, pt) + (b - pt) * (operator()(e, b) + operator()(e, pt)) / 2;
  }


 private:
	const Graph& graph; // The graph on whose edges we operate.
	BprFunction bpr; // The original BPR function.
};
