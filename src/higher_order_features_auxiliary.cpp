#include <cassert>

#include "pgmlink/higher_order_features_auxiliary.h"

namespace pgmlink {

////
//// class TrackFeaturesIdentity
////
const std::string TrackFeaturesIdentity::name_ = "TrackFeaturesIdentity";

TrackFeaturesIdentity::TrackFeaturesIdentity(
  const std::vector<std::string>& feature_names
) : feature_names_(feature_names) {
}

const std::string& TrackFeaturesIdentity::name() const {
  return name_;
}

const feature_arrays TrackFeaturesIdentity::operator()(
  const Track& track
) const {
  feature_arrays feature_matrix;
  // iterate over all traxel
  for(
    Traxelvector::const_iterator traxel_it = track.traxels_.begin();
    traxel_it != track.traxels_.end();
    traxel_it++
  ) {
    // fetch the features which are stored in a map
    FeatureMap feature_map = traxel_it->features;
    feature_array feature_vector;
    // iterate over all feature names
    for(
      std::vector<std::string>::const_iterator fname_it = feature_names_.begin();
      fname_it != feature_names_.end();
      fname_it++
    ) {
      // assert if the feature name exists
      FeatureMap::const_iterator f = feature_map.find(*fname_it);
      assert(f != feature_map.end());
      // append the features to the feature vector
      feature_vector.insert(
        feature_vector.end(),
        (f->second).begin(),
        (f->second).end()
      );
    }
    // append feature_vector to the feature matrix
    feature_matrix.push_back(feature_vector);
  }
  return feature_matrix;
}

} // namespace pgmlink
