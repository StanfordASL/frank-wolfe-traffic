#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <string>
#include <unordered_map>

#include <csv.h>

#include "DataStructures/Geometry/CoordinateConversion.h"
#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Graph/Attributes/CapacityAttribute.h"
#include "DataStructures/Graph/Attributes/CoordinateAttribute.h"
#include "DataStructures/Graph/Attributes/FreeFlowSpeedAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/NumLanesAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Attributes/VertexIdAttribute.h"
#include "Tools/LexicalCast.h"
#include "Tools/StringHelpers.h"


/*
  Input description:
  - analysisPeriod: unit is 1h
  - length: in meters
  - capacity: cars per hour
  - speed: speed in free flow (km/h)
*/

// An importer to read graphs in Visum network file format. The transport system (e.g. pedestrian,
// bicycle, car) whose network is to be imported must be supplied in the constructor. First, the
// Graph class repeatedly calls nextVertex to read the next vertex from disk and fetches various
// vertex attributes. Then, it repeatedly calls nextEdge to read the next edge from disk and
// fetches various edge attributes.
class CsvImporter {
 public:
  // Constructs an importer to read the specified system's network.
  CsvImporter(const std::string& filename,
			  const int analysisPeriod)
      : vertexFile(filename + "/vertices.csv"),
        edgeFile(filename + "/edges.csv"),
        analysisPeriod(analysisPeriod),
        nextVertexId(0)
		{
			assert(analysisPeriod > 0);
			std::cout << "Analysis period:" << analysisPeriod << std::endl;
			
		}

	// Copy constructor.
  CsvImporter(const CsvImporter& im)
      : vertexFile(im.vertexFile.get_truncated_file_name()),
        edgeFile(im.edgeFile.get_truncated_file_name()),
        analysisPeriod(im.analysisPeriod),
        nextVertexId(0) {}

	/* Opens the input file(s) and reads the header line(s).
	   Filename is given as input, although not used in practice, to keep the
	   interface  of CsvImporter idencital to the other importers.*/
  void init(const std::string& dummy_filename){
	  // used to eliminate the unused paramter warning
	  std::string dummy_string = dummy_filename;
	  dummy_string = "vert_id";
    vertexFile.read_header(io::ignore_extra_column, dummy_string, "xcoord", "ycoord");
    edgeFile.read_header(
        io::ignore_extra_column, "edge_tail", "edge_head", "length", "capacity", "speed");    
  }

  // Returns the number of vertices in the graph, or 0 if the number is not yet known.
  int numVertices() const {
    return 0;
  }

  // Returns the number of edges in the graph, or 0 if the number is not yet known.
  int numEdges() const {
    return 0;
  }

  // Reads the next vertex from disk. Returns false if there are no more vertices.
  bool nextVertex() {
	  double x,y;
	  
	  if (!vertexFile.read_row(currentVertex.id, x, y))
		  return false;

	  
	  /*const int lat = LatLng::PRECISION / 100000 * int(x);
		const int lng = LatLng::PRECISION / 100000 * int(y);*/
	  assert(origToNewIds.find(currentVertex.id) == origToNewIds.end());
	  origToNewIds[currentVertex.id] = nextVertexId++;
	  // currentVertex.latLng = LatLng(lat,lng);
	  currentVertex.latLng = LatLng(x,y);
	  return true;
  }

  // Returns the ID of the current vertex. Vertices must have sequential IDs from 0 to n − 1.
  int vertexId() const {
    return nextVertexId - 1;
  }

  // Reads the next edge from disk. Returns false if there are no more edges.
  bool nextEdge() {
    char *lengthField, *speedField;
      
	if (!edgeFile.read_row(currentEdge.edgeTail, currentEdge.edgeHead, lengthField, currentEdge.capacity, speedField))
        return false;

	//std::cout << currentEdge.edgeTail << ',' << currentEdge.edgeHead << ',' <<  lengthField << ',' << currentEdge.capacity << ',' << speedField << std::endl;
      
	assert(currentEdge.capacity >= 0);

	// assert(endsWith(speedField, "km/h"));
	currentEdge.freeFlowSpeed = lexicalCast<int>(speedField);	  
	assert(currentEdge.freeFlowSpeed >= 0);
	
    assert(origToNewIds.find(currentEdge.edgeTail) != origToNewIds.end());
    assert(origToNewIds.find(currentEdge.edgeHead) != origToNewIds.end());
    currentEdge.edgeTail = origToNewIds.at(currentEdge.edgeTail);
    currentEdge.edgeHead = origToNewIds.at(currentEdge.edgeHead);

    currentEdge.length = std::round(lexicalCast<double>(lengthField));
    assert(currentEdge.length >= 0);
    return true;
  }

  // Returns the tail vertex of the current edge.
  int edgeTail() const {
    return currentEdge.edgeTail;
  }

  // Returns the head vertex of the current edge.
  int edgeHead() const {
    return currentEdge.edgeHead;
  }

  // Returns the value of the specified attribute for the current vertex/edge, or the attribute's
  // default value if the attribute is not part of the file format.
  template <typename Attr>
  typename Attr::Type getValue() const {
    return Attr::defaultValue();
  }

  // Closes the input file(s).
  void close() { /* do nothing */ }

 private:
  // The CSV dialect used by Visum.
  template <int numFields>
  using CsvDialect = io::CSVReader<numFields>;

  // A vertex record in Visum network file format.
  struct VertexRecord {
    int id;
    Point coordinate;
    LatLng latLng;
  };

  // An edge record in Visum network file format.
  struct EdgeRecord {
    int edgeTail;
    int edgeHead;
    int length;
    int capacity;
    int freeFlowSpeed;
  };

  using IdMap = std::unordered_map<int, int>;

  CsvDialect<3> vertexFile;          // The CSV file containing the vertex records.
  CsvDialect<5> edgeFile;            // The CSV file containing the edge records.
  const int analysisPeriod;          // The analysis period in hours (capacity is in vehicles/AP).

  IdMap origToNewIds; // A map from original vertex IDs to new sequential IDs.
  int nextVertexId;   // The next free vertex ID.

  VertexRecord currentVertex; // The vertex record read by the last call of nextVertex.
  EdgeRecord currentEdge;     // The edge record read by the last call of nextEdge.
};

// Returns the value of the coordinate attribute for the current vertex.
template <>
inline CoordinateAttribute::Type CsvImporter::getValue<CoordinateAttribute>() const {
  return currentVertex.coordinate;
}

// Returns the value of the LatLng attribute for the current vertex.
template <>
inline LatLngAttribute::Type CsvImporter::getValue<LatLngAttribute>() const {
  return currentVertex.latLng;
}

// Returns the value of the capacity attribute for the current edge.
template <>
inline CapacityAttribute::Type CsvImporter::getValue<CapacityAttribute>() const {
  return std::round(1.0 * currentEdge.capacity / analysisPeriod);
}

// Returns the value of the free-flow speed attribute for the current edge.
template <>
inline FreeFlowSpeedAttribute::Type CsvImporter::getValue<FreeFlowSpeedAttribute>() const {
  return currentEdge.freeFlowSpeed;
}

// Returns the value of the length attribute for the current edge.
template <>
inline LengthAttribute::Type CsvImporter::getValue<LengthAttribute>() const {
  return currentEdge.length;
}

/*// Returns the value of the number of lanes attribute for the current edge.
template <>
inline NumLanesAttribute::Type CsvImporter::getValue<NumLanesAttribute>() const {
  return currentEdge.numLanes;
  }*/

// Returns the value of the travel time attribute for the current edge.
// This is the time it takes to traverse edge in free flow
template <>
inline TravelTimeAttribute::Type CsvImporter::getValue<TravelTimeAttribute>() const {
  return std::round(360.0 * currentEdge.length / currentEdge.freeFlowSpeed);
}

// Returns the value of the vertex ID attribute for the current vertex.
template <>
inline VertexIdAttribute::Type CsvImporter::getValue<VertexIdAttribute>() const {
  return currentVertex.id;
}
