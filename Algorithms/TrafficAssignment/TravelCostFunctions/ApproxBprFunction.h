#pragma once

#include <vectorclass/vectorclass.h>
#include "Algorithms/TrafficAssignment/TravelCostFunctions/BprFunction.h"
#include "Algorithms/TrafficAssignment/TravelCostFunctions/CustomBprFunction.h"

// The BPR travel cost function, relating the travel time on an edge to the flow on this edge.
template <typename GraphT>
class ApproxBprFunction {
 public:
  // Constructs a BPR function.
	ApproxBprFunction(const GraphT& graph) : graph(graph), bpr(graph), cbpr(graph) {
		exo = graph.exogenous();
		dummy_id = graph.dummyId();
		exo_v = Vec4d(exo);
        
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
	  double x_new;
	  const double cap = graph.capacity(e);
	  
	  if (graph.edgeReal(e))
	  {
		  x_new = x + exo*cap;
		  const double pt = XTH * cap;
		  if (x_new <= pt)
			  return bpr(e, exo*cap);
		  else
		  {
			  double b_slope = (bpr(e, 2*cap) - bpr(e, pt)) / (2*cap - pt);
			  double cost =  b_slope * (x_new - pt) +  bpr(e, exo*cap);
			  return cost;
		  }
	  }
	  else
	  {
		  x_new = x;
		  const double pt = APT * cap;
		  if (x_new <= pt)
			  return bpr(e, x_new);
		  else
			  return bpr(e, pt) + bpr.derivative(e, pt) * (x_new - pt);
	  }
  }

  // Returns the derivative of e's travel cost function at x.
  double derivative(const int e, const double x) const {
	  double x_new;
	  const double cap = graph.capacity(e);
	  if (graph.edgeReal(e))
	  {
		  x_new = x + exo*cap;
		  const double pt =XTH * cap;
		  if (x <= pt)
			  return 0.0;
		  else
		  {
			  double b_slope = (bpr(e, 2*cap) - bpr(e, pt)) / (2*cap - pt);
			  return b_slope;
		  }
	  }
	  else
	  {
		  x_new = x;
		  const double pt = APT * graph.capacity(e); // The point at which we linearize.
		  if (x_new <= pt)
			  return bpr.derivative(e, x_new);
		  else
			  return bpr.derivative(e, pt);
	  }	  
  }

  // Returns the integral of e's travel cost function from 0 to b.
  double integral(const int e, const double b) const {
	  // this is not supposed to happen
	  std::cout << "INTEGRAL OF APPROXBPR IS USED " << b << e << std::endl;	  
	  return 0.0;
  }

  // Returns the travel times on four consecutive edges starting at e, given the flows x on them.
  Vec4d operator()(const int e, const Vec4d& x) const {
	  const double cap = graph.capacity(e);
	  const Vec4d cap_v(cap);
	  Vec4d edgeReal, edgeFake;
	  edgeReal.load(&graph.edgeReal(e));
	  edgeFake.load(&graph.edgeFake(e));
	  Vec4d val_fake = cbpr(e,x)*edgeFake;

	  Vec4d x_new;
	  x_new = x + exo_v*cap_v;

	  Vec4d pt_real = XTH * cap_v;
	  Vec4d slope = (bpr(e,2*cap_v)-bpr(e,pt_real))/(2*cap-pt_real);
	  Vec4d val_real = (bpr(e, exo_v*cap_v) + slope * (max(x_new, pt_real) - pt_real)) * edgeReal;
	  
	  return val_fake + val_real;
  }

  // Returns the derivative of e's travel cost function at x.
  Vec4d derivative(const int e, const Vec4d& x) const{
	  const double cap = graph.capacity(e);
	  const Vec4d cap_v(cap);
	  Vec4d edgeReal, edgeFake;
	  edgeReal.load(&graph.edgeReal(e));
	  edgeFake.load(&graph.edgeFake(e));
	  Vec4d val_fake = cbpr.derivative(e,x)*edgeFake;

	  Vec4d x_new;
	  x_new = x + exo_v*cap_v;
	  Vec4d pt_real = XTH * cap_v;
	  Vec4d slope = (bpr(e,2*cap_v)-bpr(e,pt_real))/(2*cap-pt_real);
	  Vec4d val_real = slope *  ((min(x_new,pt_real)-x_new) / (pt_real-x_new)) * edgeReal;

	  return val_fake + val_real; 
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
	BprFunction<GraphT> bpr; // The original BPR function.
	CustomBprFunction<GraphT> cbpr; // The original BPR function.
	double exo;
	int dummy_id;
	Vec4d exo_v;
    
    std::vector<double> edgeTotalShift;//Lucas
    std::vector<double> edgeRebalancers;
    AlignedVector<double> trafficFlows;
};
