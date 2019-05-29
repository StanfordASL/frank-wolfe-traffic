#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"
#include "Tools/Constants.h"

// An attribute associating a travel cost with each edge of a graph.
class EdgeNegativeShiftAttribute : public AbstractAttribute<double> {
 public:
  // Returns the attribute's default value.
  static Type defaultValue() {
    return 0.0;
  }

  // Returns the travel cost on edge e.
  const Type& edgeNegativeShift(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the travel cost on edge e.
  Type& edgeNegativeShift(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

 protected:
  static constexpr const char* NAME = "edgeNegativeShift"; // The attribute's unique name.
};
