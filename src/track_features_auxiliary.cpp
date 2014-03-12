#include "pgmlink/track_features_auxiliary.h"

namespace pgmlink{
////
//// Class FeatureAggregator
////
const std::string FeatureAggregator::name_ = "";

feature_array FeatureAggregator::vector_valued(
  const feature_arrays features
) const {
  throw std::runtime_error(
    "FeatureAggregator \"" + name() + "\" has no vector valued method"
  );
  feature_array ret;
  return ret;
}

feature_type FeatureAggregator::scalar_valued(
  const feature_arrays features
) const {
  throw std::runtime_error(
    "FeatureAggregator \"" + name() + "\" has no scalar valued method"
  );
  feature_type ret;
  return ret;
}

const std::string& FeatureAggregator::name() const {
  return FeatureAggregator::name_;
}

} // namespace pgmlink
