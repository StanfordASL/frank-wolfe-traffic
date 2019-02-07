#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>

#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Labels/SimdLabelSet.h"
#include "Tools/Simd/AlignedVector.h"

namespace trafficassignment {

// An adapter that makes Dijkstra's algorithm usable in the all-or-nothing assignment procedure.
template <typename InputGraphT, typename WeightT>
class DijkstraAdapter {
 public:
#if TA_LOG_K < 2 || defined(TA_NO_SIMD_SEARCH)
  using LabelSet = BasicLabelSet<TA_LOG_K, ParentInfo::FULL_PARENT_INFO>;
#else
  using LabelSet = SimdLabelSet<TA_LOG_K, ParentInfo::FULL_PARENT_INFO>;
#endif
  using InputGraph = InputGraphT;

  static constexpr int K = LabelSet::K; // The number of simultaneous shortest-path computations.

  // The search algorithm using the graph and possibly auxiliary data to compute shortest paths.
  // Multiple instances can work on the same data concurrently.
  class QueryAlgo {
   public:
    // Constructs a query algorithm instance working on the specified data.
	  QueryAlgo(const InputGraph& inputGraph, AlignedVector<int>& flowsOnForwardEdges, const int odNum)
        : search(inputGraph),
          flowsOnForwardEdges(flowsOnForwardEdges),
          localFlowsOnForwardEdges(flowsOnForwardEdges.size()) {
      assert(inputGraph.numEdges() == flowsOnForwardEdges.size());
	  (void)odNum;
    }

    // Computes shortest paths from each source to its target simultaneously.
	  void run(std::array<int, K>& sources, std::array<int, K>& targets, const int k, const bool consider_loss, const int first_k) {
		// loss is not implemented here
		assert(!consider_loss);
		(void)consider_loss;
		(void)first_k;
		
		
		// Run a centralized Dijkstra search.
		search.run(sources, targets);

		// Assign flow to the edges on the computed paths.
		for (auto i = 0; i < k; ++i) {
			for (const auto e : search.getReverseEdgePath(targets[i], i)) {
				assert(e >= 0); assert(e < localFlowsOnForwardEdges.size());
				++localFlowsOnForwardEdges[e];
			}
		}
    }

    // Returns the length of the i-th shortest path.
    int getDistance(const int dst, const int i) {
      return search.getDistance(dst, i);
    }

    // Adds the local flow counters to the global ones. Must be synchronized externally.
    void addLocalToGlobalFlows(const bool consider_loss = false) {
		// loss is not implemented here
		assert(!consider_loss);
		if (consider_loss)
			std::cout << "This is not supposed to happen" << std::endl;
		
		for (auto e = 0; e < flowsOnForwardEdges.size(); ++e)
			flowsOnForwardEdges[e] += localFlowsOnForwardEdges[e];
    }

   private:
    using Dijkstra = StandardDijkstra<InputGraph, WeightT, LabelSet>;

    Dijkstra search;                           // The Dijkstra search.
    AlignedVector<int>& flowsOnForwardEdges;   // The flows in the forward graph.
    std::vector<int> localFlowsOnForwardEdges; // The local flows in the forward graph.
  };

  // Constructs an adapter for Dijkstra's algorithm.
	explicit DijkstraAdapter(const InputGraph& inputGraph, const int odNum)
		: inputGraph(inputGraph), flowsOnForwardEdges(inputGraph.numEdges()), odNum(odNum) {}

  // Invoked before the first iteration.
  void preprocess() { /* do nothing */ }

  // Invoked before each iteration.
  void customize() {
    std::fill(flowsOnForwardEdges.begin(), flowsOnForwardEdges.end(), 0);
  }

  // Returns an instance of the query algorithm.
  QueryAlgo getQueryAlgoInstance() {
	  return {inputGraph, flowsOnForwardEdges, odNum};
  }

  // Propagates the flows on the edges in the search graphs to the edges in the input graph.
  void propagateFlowsToInputEdges(AlignedVector<int>& flowsOnInputEdges) {
    assert(flowsOnInputEdges.size() == flowsOnInputEdges.size());
    flowsOnInputEdges.swap(flowsOnForwardEdges);
  }

 private:
  const InputGraph& inputGraph;           // The input graph.
  AlignedVector<int> flowsOnForwardEdges; // The flows on the edges in the forward graph.
	const int odNum;
	
};

}
