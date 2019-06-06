#pragma once

#include <vectorclass/vectorclass.h>

#define APT 3.0 // Global linearization point
#define XTH 1.2
#define XEND 2

// The BPR travel cost function, relating the travel time on an edge to the flow on this edge.
template <typename GraphT>
class BprFunction {
 public:
  // Constructs a BPR function.
  BprFunction(const GraphT& graph) : graph(graph) {
      //Lucas
  	  std::vector<double> vect(graph.numEdges(),0);
      edgeTotalShift=vect;//initialized to zero by default
      std::vector<double> vect2(graph.numEdges(),-1);
      edgeRebalancers=vect2;
      AlignedVector<double> vect3(graph.numEdges(),0);
      trafficFlows=vect3;
  }

  // Returns the travel time on edge e, given the flow x on e.
  double operator()(const int e, const double x) const {
    const double tmp = x / graph.capacity(e);
    return graph.travelTime(e) * (1 + 0.15 * tmp * tmp * tmp * tmp);
  }

  // Returns the derivative of e's travel cost function at x.
  double derivative(const int e, const double x) const {
    const double tmp = x / graph.capacity(e);
    return graph.travelTime(e) * 0.15 * 4 * tmp * tmp * tmp / graph.capacity(e);
  }

  // Returns the antiderivative of e's travel cost function at x.
  double antiderivative(const int e, const double x) const {
    const double tmp = x / graph.capacity(e);
    return graph.travelTime(e) * (x + 0.15 * x * tmp * tmp * tmp * tmp / (4 + 1));
  }

  // Returns the integral of e's travel cost function from 0 to b.
  double integral(const int e, const double b) const {
    return antiderivative(e, b) - antiderivative(e, 0);
  }

  // Returns the travel times on four consecutive edges starting at e, given the flows x on them.
  Vec4d operator()(const int e, const Vec4d& x) const {
    Vec4d time = to_double(Vec4i().load(&graph.travelTime(e)));
    Vec4d capacity = to_double(Vec4i().load(&graph.capacity(e)));
    return time * (1 + 0.15 * pow_const(x / capacity, 4));
  }

  // Returns the derivative of e's travel cost function at x.
  Vec4d derivative(const int e, const Vec4d& x) const {
    Vec4d time = to_double(Vec4i().load(&graph.travelTime(e)));
    Vec4d capacity = to_double(Vec4i().load(&graph.capacity(e)));
    return time * 0.15 * 4 * pow_const(x / capacity, 4 - 1) / capacity;
  }
    
    //Lucas
    void setEdgeShift(std::vector<double> inputVectorShift){//Accessor to edit the vectorshift
        edgeTotalShift = inputVectorShift;
    }

    
    double getEdgeShift(const int e){
      return edgeTotalShift[e];
    }
    
    void setEdgeRebalancers(std::vector<double> inputEdgeRebalancers){//Accessor to edit the edgerebalancers
        edgeRebalancers = inputEdgeRebalancers;
    }
    
    void updateTrafficFlows(AlignedVector<double> inputTrafficFlows){
        trafficFlows=inputTrafficFlows;
    }

 private:
  const GraphT& graph; // The graph on whose edges we operate.
    
  std::vector<double> edgeTotalShift;//Lucas
    std::vector<double> edgeRebalancers;
    AlignedVector<double> trafficFlows;
};
