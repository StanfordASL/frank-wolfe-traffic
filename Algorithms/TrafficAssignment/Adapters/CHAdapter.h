#pragma once

#include <array>
#include <cassert>
#include <vector>

#include "Algorithms/CH/CH.h"
#include "Algorithms/CH/CHQuery.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Labels/SimdLabelSet.h"
#include "Tools/Simd/AlignedVector.h"

namespace trafficassignment {

// An adapter that makes CHs usable in the all-or-nothing assignment procedure.
template <typename InputGraphT, typename WeightT>
class CHAdapter {
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
    QueryAlgo(
        const CH& ch, AlignedVector<int>& flowsOnUpEdges, AlignedVector<int>& flowsOnDownEdges, const int odNum)
        : ch(ch),
          search(ch),
          flowsOnUpEdges(flowsOnUpEdges),
          flowsOnDownEdges(flowsOnDownEdges),
          localFlowsOnUpEdges(flowsOnUpEdges.size()),
          localFlowsOnDownEdges(flowsOnDownEdges.size()), odPathsUp(odNum), odPathsDown(odNum), odLengths(odNum){
      assert(ch.upwardGraph().numEdges() == flowsOnUpEdges.size());
      assert(ch.downwardGraph().numEdges() == flowsOnDownEdges.size());
    }

    // Computes shortest paths from each source to its target simultaneously.
	  void run(std::array<int, K>& sources, std::array<int, K>& targets, const int k, const bool consider_loss, const int first_k) {
		// Run a centralized CH search.
		for (auto i = 0; i < K; ++i) {
			sources[i] = ch.rank(sources[i]);
			targets[i] = ch.rank(targets[i]);
		}
		search.run(sources, targets);
		
		// Assign flow to the edges on the computed paths.
		if (!consider_loss){
			for (auto i = 0; i < k; ++i) {				
				for (const auto e : search.getUpEdgePath(i)) {
					assert(e >= 0); assert(e < localFlowsOnUpEdges.size());
					++localFlowsOnUpEdges[e];
				}
				for (const auto e : search.getDownEdgePath(i)) {
					assert(e >= 0); assert(e < localFlowsOnDownEdges.size());
					++localFlowsOnDownEdges[e];
				}
			}
		} else {
			for (auto i = 0; i < k; ++i) {
				for (const auto e : search.getUpEdgePath(i)) 
					odPathsUp[i+first_k].push_back(e);
				for (const auto e : search.getDownEdgePath(i))
					odPathsDown[i+first_k].push_back(e);
				
				odLengths[i+first_k] = search.getDistance(i);
			}
		}
	}

    // Returns the length of the i-th shortest path.
    int getDistance(const int /*dst*/, const int i) {
      return search.getDistance(i);
    }

    // Adds the local flow counters to the global ones. Must be synchronized externally.
    void addLocalToGlobalFlows(const bool consider_loss = false) {
		if (consider_loss){
			int N = odLengths.size();
			int n = N / 3;
			std::cout << "N=" << N << std::endl;
			assert(N % 3 == 0);
			
			std::vector<int> to_process;
			for (auto i = 0; i < n; i++){
				int cost_rebalance = odLengths[i] + odLengths[i + n];
				int cost_loss = odLengths[i + 2*n];
				
				if (cost_rebalance < cost_loss){
					to_process.push_back(i);
					to_process.push_back(i+n);
				}
				else
					to_process.push_back(i+2*n);
			}

			for (auto it = to_process.begin(); it != to_process.end(); it++){
				for (const auto e : odPathsUp[*it]) {
					assert(e >= 0); assert(e < localFlowsOnUpEdges.size());
					++localFlowsOnUpEdges[e];
				}
				for (const auto e : odPathsDown[*it]) {
					assert(e >= 0); assert(e < localFlowsOnDownEdges.size());
					++localFlowsOnDownEdges[e];
				}
			}
		}
		
		FORALL_EDGES(ch.upwardGraph(), e)
			flowsOnUpEdges[e] += localFlowsOnUpEdges[e];
		FORALL_EDGES(ch.downwardGraph(), e)
			flowsOnDownEdges[e] += localFlowsOnDownEdges[e];
    }

   private:
	  const CH& ch;                           // The CH rebuilt in each iteration.
	  StandardCHQuery<LabelSet> search;       // The CH search.
	  AlignedVector<int>& flowsOnUpEdges;     // The flows in the upward graph.
	  AlignedVector<int>& flowsOnDownEdges;   // The flows in the downward graph.
	  std::vector<int> localFlowsOnUpEdges;   // The local flows in the upward graph.
	  std::vector<int> localFlowsOnDownEdges; // The local flows in the downward graph.
	  std::vector<std::vector<int>> odPathsUp;  // The actual paths (used only for LOSS)
	  std::vector<std::vector<int>> odPathsDown;  // The actual paths (used only for LOSS)
	  std::vector<int> odLengths;          // Lengths of paths
  };

  // Constructs an adapter for CHs.
	explicit CHAdapter(const InputGraph& inputGraph, const int odNum) : inputGraph(inputGraph), odNum(odNum) {}

  // Invoked before the first iteration.
  void preprocess() { /* do nothing */ }

  // Invoked before each iteration.
  void customize() {
    ch.preprocess<WeightT>(inputGraph);
    flowsOnUpEdges.assign(ch.upwardGraph().numEdges(), 0);
    flowsOnDownEdges.assign(ch.downwardGraph().numEdges(), 0);
  }

  // Returns an instance of the query algorithm.
  QueryAlgo getQueryAlgoInstance() {
	  return {ch, flowsOnUpEdges, flowsOnDownEdges, odNum};
  }

  // Propagates the flows on the edges in the search graphs to the edges in the input graph.
  void propagateFlowsToInputEdges(AlignedVector<int>& flowsOnInputEdges) {
    const auto& upGraph = ch.upwardGraph();
    const auto& downGraph = ch.downwardGraph();
    for (auto u = inputGraph.numVertices() - 1; u >= 0; --u) {
      FORALL_INCIDENT_EDGES(upGraph, u, e)
        if (upGraph.unpackingInfo(e).second == INVALID_EDGE) {
          flowsOnInputEdges[upGraph.unpackingInfo(e).first] = flowsOnUpEdges[e];
        } else {
          flowsOnDownEdges[upGraph.unpackingInfo(e).first] += flowsOnUpEdges[e];
          flowsOnUpEdges[upGraph.unpackingInfo(e).second] += flowsOnUpEdges[e];
        }
      FORALL_INCIDENT_EDGES(downGraph, u, e)
        if (downGraph.unpackingInfo(e).second == INVALID_EDGE) {
          flowsOnInputEdges[downGraph.unpackingInfo(e).first] = flowsOnDownEdges[e];
        } else {
          flowsOnDownEdges[downGraph.unpackingInfo(e).first] += flowsOnDownEdges[e];
          flowsOnUpEdges[downGraph.unpackingInfo(e).second] += flowsOnDownEdges[e];
        }
    }
  }

 private:
  const InputGraph& inputGraph; // The input graph.
  CH ch;                        // The weighted CH rebuilt in each iteration.

  AlignedVector<int> flowsOnUpEdges;   // The flows on the edges in the upward graph.
  AlignedVector<int> flowsOnDownEdges; // The flows on the edges in the downward graph.
	const int odNum;
};

}
