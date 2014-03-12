#ifndef PGMLINK_TRACK_FEATURES_H
#define PGMLINK_TRACK_FEATURES_H

#include <vector>

#include "pgmlink/traxels.h"

namespace pgmlink {

class FeatureAggregator {
 public:
  virtual ~FeatureAggregator() = 0; // rein virtuell
  virtual feature_array vector_valued(const feature_arrays features) const;
  virtual feature_type scalar_valued(const feature_arrays features) const;
  virtual const std::string& name() const;

  static const std::string name_;
}; // Class FeatureAggregator

} // Namespace pgmlink

#endif // PGMLINK_TRACK_FEATURES_H
