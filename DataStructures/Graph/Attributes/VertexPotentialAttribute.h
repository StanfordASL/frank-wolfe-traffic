#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"
#include "Tools/Constants.h"

// An attribute associating an ID with each vertex of a graph.
class VertexPotentialAttribute : public AbstractAttribute<double> {
 public:
  // Returns the attribute's default value.
  static Type defaultValue() {
    return 0;
  }

  // Returns the ID of vertex v.
  const Type& vertexPotential(const int v) const {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

  // Returns a reference to the ID of vertex v.
  Type& vertexPotential(const int v) {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

 protected:
  static constexpr const char* NAME = "vertex_potential"; // The attribute's unique name.
};
