#pragma once

#include <vectorclass/vectorclass.h>
#include "Algorithms/TrafficAssignment/TravelCostFunctions/BprFunction.h"

// The BPR travel cost function, relating the travel time on an edge to the flow on this edge.
template <typename GraphTT>
class CustomBprFunction {
 public:
  // Constructs a BPR function.
	CustomBprFunction(const GraphTT& graph) : graph(graph), bpr(graph) {
		exo = graph.exogenous();
		dummy_id = graph.dummyId();
		exo_v = Vec4d(exo);
	}

  // Returns the travel time on edge e, given the flow x on e.
  double operator()(const int e, const double x) const {
	  double x_new;
	  const double cap = graph.capacity(e);
	  
	  if (graph.edgeReal(e))
		  x_new = x + exo*cap;
	  else
		  x_new = x;
	  const double pt = APT * cap;
	  if (x <= pt)
		  return bpr(e, x_new);
	  else
		  return bpr(e, pt) + bpr.derivative(e, pt) * (x_new - pt);
  }

  // Returns the derivative of e's travel cost function at x.
  double derivative(const int e, const double x) const {
	  double x_new;
	  const double cap = graph.capacity(e);
	  if (graph.edgeReal(e))
		  x_new = x + exo*cap;
	  else
		  x_new = x;

	  const double pt = APT * graph.capacity(e); // The point at which we linearize.
	  if (x_new <= pt)
		  return bpr.derivative(e, x_new);
	  else
		  return bpr.derivative(e, pt);
  }

  // Returns the integral of e's travel cost function from 0 to b.
  double integral(const int e, const double b) const {
	double b_new;
	const double cap = graph.capacity(e);
	
	if (graph.edgeReal(e))
	  b_new = b + exo*cap;
	else
	  b_new = b;
	
    const double pt = APT * graph.capacity(e); // The point at which we linearize.
	double int_old, int_new;
	
    if (b_new <= pt)
      int_new = bpr.integral(e, b_new);
    else
      int_new = bpr.integral(e, pt) + (b_new - pt) * (operator()(e, b_new) + operator()(e, pt)) / 2;

	if (b <= pt)
		int_old =  bpr.integral(e, b);
	else
		int_old = bpr.integral(e, pt) + (b - pt) * (operator()(e, b) + operator()(e, pt)) / 2;

	return int_new - int_old;		
  }

  // Returns the travel times on four consecutive edges starting at e, given the flows x on them.
  Vec4d operator()(const int e, const Vec4d& x) const {
	  const double cap = graph.capacity(e);
	  const Vec4d cap_v(cap);
	  Vec4d edgeReal;
	  edgeReal.load(&graph.edgeReal(e));
	  
	  Vec4d x_new;
	  x_new = x + exo_v*cap_v*edgeReal;
	  
    Vec4d pt = APT * to_double(Vec4i().load(&graph.capacity(e)));
    return bpr(e, min(x_new, pt)) +
        bpr.derivative(e, pt) * (max(x_new, pt) - pt);
  }

  // Returns the derivative of e's travel cost function at x.
  Vec4d derivative(const int e, const Vec4d& x) const{
	  const double cap = graph.capacity(e);
	  const Vec4d cap_v(cap);
	  Vec4d edgeReal;
	  edgeReal.load(&graph.edgeReal(e));
	  
	  Vec4d x_new;
  	  x_new = x + (exo_v*cap_v)*edgeReal;
	  
	  Vec4d pt = APT * to_double(Vec4i().load(&graph.capacity(e)));
	  return bpr.derivative(e, min(x_new, pt));
  }

 private:
	const GraphTT& graph; // The graph on whose edges we operate.
	BprFunction<GraphTT> bpr; // The original BPR function.
	double exo;
	int dummy_id;
	Vec4d exo_v;
};
