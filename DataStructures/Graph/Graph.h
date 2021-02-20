#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <ostream>
#include <iostream>
#include <cstdio>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <csv.h>

class Graph
{
public:
	// Constructs a graph from csv edge file
	Graph(const std::string& filename, const double ceParameter) : vertexNum(0), ceParameter(ceParameter) {
		vertexNum = 0;
		readFrom(filename);
	}
							  
	// Returns the number of vertices in the graph
	int numVertices() const {
		return vertexNum;
	}

	// Returns the number of edges in the graph.
	int numEdges() const {
		return edgeHead.size();
	}
	
	// Returns the head vertex of edge e.
	int head(const int e) const {
		assert(e >= 0);
		assert(e < edgeHead.size()); 
		return edgeHead[e];
	}

	// Returns the head vertex of edge e.
	int tail(const int e) const {
		assert(e >= 0);
		assert(e < edgeTail.size()); 
		return edgeTail[e];
	}

	// Returns the capacity of edge e.
	int capacity(const int e) const {
		assert(e >= 0);
		assert(e < edgeCapacity.size()); 
		return edgeCapacity[e];
	}

	// Returns the length of edge e.
	int length(const int e) const {
		assert(e >= 0);
		assert(e < edgeLength.size()); 
		return edgeLength[e];
	}

	// Returns the speed of edge e.
	int speed(const int e) const {
		assert(e >= 0);
		assert(e < edgeSpeed.size()); 
		return edgeSpeed[e];
	}

	double weight(const int e) const 
	{
		assert(e >= 0);
		assert(e < edgeWeight.size()); 
		return edgeWeight[e];
	}
	
		
	// Returns the travel time of edge e at free flow.
	double freeTravelTime(const int e) const {
		assert(e >= 0);
		assert(e < edgeFreeTravelTime.size()); 
		return edgeFreeTravelTime[e];
	}

	// Sets edge weight, which will be used in the AoN routine
	void setWeight(const int e, const double v)
	{
		assert(e >= 0);
		assert(e < edgeWeight.size());
		edgeWeight[e] = v;
	}

	std::vector<double>& getWeights()
	{
		return edgeWeight;
	}

	double combinedEquilibriumParameter()
	{
		return ceParameter;
	}
	
private:
	template <int numFields>
	using CsvDialect = io::CSVReader<numFields>;
	
	void readFrom(const std::string& filename) {		
		CsvDialect<5> edgeFile(filename);            // The CSV file containing the edge records.
		edgeFile.read_header(io::ignore_extra_column, "edge_tail", "edge_head", "length", "capacity", "speed");

		int tail, head, length, capacity, speed;
		
		while(edgeFile.read_row(tail, head, length, capacity, speed))
		{
			edgeTail.push_back(tail);
			edgeHead.push_back(head);
			edgeCapacity.push_back(capacity);
			edgeLength.push_back(length);
			edgeSpeed.push_back(speed);

			assert(tail >= 0);
			assert(head >= 0);
			assert(capacity > 0);
			assert(length > 0);
			assert(speed > 0);

			// update vertex number
			if (tail + 1 > vertexNum || head + 1 > vertexNum)
				vertexNum = std::max(tail, head) + 1;

			// compute free flow travel time in minutes
			double freeFlowTravelTime = 60 * 60 * ((double) length / 1000.0) / ((double) speed);
			edgeFreeTravelTime.push_back(freeFlowTravelTime);
			edgeWeight.push_back(freeFlowTravelTime); // initial edge weight
		}
	}
	
	int vertexNum;
	std::vector<int> edgeTail;
	std::vector<int> edgeHead;
	std::vector<int> edgeCapacity; // vehicles / h
	std::vector<int> edgeLength; // length (m)
	std::vector<int> edgeSpeed; // travel time in free flow (k/h)
	std::vector<double> edgeFreeTravelTime; // hours
	std::vector<double> edgeWeight; // current travel time, in hours

	double ceParameter; // paramter for combined equilibrium calculation
};

// Iteration macros for conveniently looping through vertices or edges of a graph.
#define FORALL_VERTICES(G, u) for (int u = 0; u < G.numVertices(); ++u)
#define FORALL_EDGES(G, e) for (int e = 0; e < G.numEdges(); ++e)
// #define FORALL_INCIDENT_EDGES(G, u, e) for (int e = G.firstEdge(u); e < G.lastEdge(u); ++e)
