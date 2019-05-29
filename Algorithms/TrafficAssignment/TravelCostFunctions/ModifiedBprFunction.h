#pragma once

#include <vectorclass/vectorclass.h>
#include "Algorithms/TrafficAssignment/TravelCostFunctions/BprFunction.h"

// The BPR travel cost function, relating the travel time on an edge to the flow on this edge.
template <typename GraphT>
class ModifiedBprFunction {
 public:
  // Constructs a BPR function.
	ModifiedBprFunction(const GraphT& graph) : graph(graph), bpr(graph) {
		std::cout << "here" << std::endl;
       //Lucas
      std::vector<double> vect(graph.numEdges(),0);
      edgeTotalShift=vect;//initialized to zero by default
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

  // Returns the integral of e's travel cost function from 0 to b.
  double integral(const int e, const double b) const {
    const double pt = APT * graph.capacity(e); // The point at which we linearize.
    if (b <= pt)
      return bpr.integral(e, b);
    else
      return bpr.integral(e, pt) + (b - pt) * (operator()(e, b) + operator()(e, pt)) / 2;
  }

  // Returns the travel times on four consecutive edges starting at e, given the flows x on them.
  Vec4d operator()(const int e, const Vec4d& x) const {
    Vec4d pt = APT * to_double(Vec4i().load(&graph.capacity(e)));
    return bpr(e, min(x, pt)) +
        bpr.derivative(e, pt) * (max(x, pt) - pt);
  }

  // Returns the derivative of e's travel cost function at x.
  Vec4d derivative(const int e, const Vec4d& x) const {
    Vec4d pt =APT * to_double(Vec4i().load(&graph.capacity(e)));
    return bpr.derivative(e, min(x, pt));
  }
    
    //Lucas
    void setEdgeShift(std::vector<double> inputVectorShift){//Accessor to edit the vectorshift
        edgeTotalShift = inputVectorShift;
    }

 private:
	const GraphT& graph; // The graph on whose edges we operate.
	BprFunction<GraphT> bpr; // The original BPR function.
    
    std::vector<double> edgeTotalShift;//Lucas
};
