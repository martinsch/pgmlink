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


////
//// class TrackFeaturesDiff
////
const std::string TrackFeaturesDiff::name_ = "TrackFeaturesDiff";

TrackFeaturesDiff::TrackFeaturesDiff(
  const std::vector<std::string>& feature_names
) : feature_names_(feature_names) {
}

const std::string& TrackFeaturesDiff::name() const {
  return name_;
}

const feature_arrays TrackFeaturesDiff::operator()(
  const Track& track
) const {
  feature_arrays feature_matrix;
  // iterate over all traxel
  assert(track.traxels_.size() >= 2);
  Traxelvector::const_iterator t1_it, t2_it;
  t1_it = track.traxels_.begin()+1;
  t2_it = track.traxels_.begin();
  for(; t1_it != track.traxels_.end(); t1_it++, t2_it++) {
    // fetch the features which are stored in a map
    const FeatureMap& fmap1 = t1_it->features;
    const FeatureMap& fmap2 = t2_it->features;
    feature_array feature_vector;
    // iterate over all feature names
    for(
      std::vector<std::string>::const_iterator fname_it = feature_names_.begin();
      fname_it != feature_names_.end();
      fname_it++
    ) {
      // check if the feature name exists
      FeatureMap::const_iterator fmap1_it = fmap1.find(*fname_it);
      FeatureMap::const_iterator fmap2_it = fmap2.find(*fname_it);
      assert(fmap1_it != fmap1.end());
      assert(fmap2_it != fmap2.end());
      // check if the feature vectors have the same size
      assert((fmap1_it->second).size() == (fmap2_it->second).size());
      // append the differences to the feature vector
      feature_array::const_iterator farray1_it, farray2_it;
      farray1_it = (fmap1_it->second).begin();
      farray2_it = (fmap2_it->second).begin();
      for (
        ;
        farray1_it != (fmap1_it->second).end();
        farray1_it++, farray2_it++
      ) {
        feature_vector.push_back(*farray1_it - *farray2_it);
      }
    }
    // append feature_vector to the feature matrix
    feature_matrix.push_back(feature_vector);
  }
  return feature_matrix;
}

} // namespace pgmlink
