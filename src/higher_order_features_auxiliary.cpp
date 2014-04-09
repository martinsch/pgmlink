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

TrackFeaturesIdentity::TrackFeaturesIdentity(
  const std::string& feature_name
) {
  feature_names_.resize(1);
  feature_names_[0] = feature_name;
}

const std::string& TrackFeaturesIdentity::name() const {
  return name_;
}

feature_arrays TrackFeaturesIdentity::operator()(
  const Track& track
) {
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

TrackFeaturesDiff::TrackFeaturesDiff(
  const std::string& feature_name
) {
  feature_names_.resize(1);
  feature_names_[0] = feature_name;
}

const std::string& TrackFeaturesDiff::name() const {
  return name_;
}

feature_arrays TrackFeaturesDiff::operator()(
  const Track& track
) {
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

////
//// class TrackFeaturesCurvature
////
const std::string TrackFeaturesCurvature::name_ = "TrackFeaturesCurvature";

TrackFeaturesCurvature::TrackFeaturesCurvature(
  const std::vector<std::string>& feature_names
) : feature_names_(feature_names) {
}

TrackFeaturesCurvature::TrackFeaturesCurvature(
  const std::string& feature_name
) {
  feature_names_.resize(1);
  feature_names_[0] = feature_name;
}

const std::string& TrackFeaturesCurvature::name() const {
  return name_;
}

feature_arrays TrackFeaturesCurvature::operator()(
  const Track& track
) {
  feature_arrays feature_matrix;
  // iterate over all traxel
  assert(track.traxels_.size() >= 3);
  Traxelvector::const_iterator t1_it, t2_it, t3_it;
  t1_it = track.traxels_.begin()+2;
  t2_it = track.traxels_.begin()+1;
  t3_it = track.traxels_.begin();
  for(; t1_it != track.traxels_.end(); t1_it++, t2_it++, t3_it++) {
    // fetch the features which are stored in a map
    const FeatureMap& fmap1 = t1_it->features;
    const FeatureMap& fmap2 = t2_it->features;
    const FeatureMap& fmap3 = t3_it->features;
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
      FeatureMap::const_iterator fmap3_it = fmap3.find(*fname_it);
      assert(fmap1_it != fmap1.end());
      assert(fmap2_it != fmap2.end());
      assert(fmap3_it != fmap3.end());
      // check if the feature vectors have the same size
      assert((fmap1_it->second).size() == (fmap2_it->second).size());
      assert((fmap1_it->second).size() == (fmap3_it->second).size());
      // append the curvature to the feature vector
      feature_array::const_iterator farray1_it, farray2_it, farray3_it;
      farray1_it = (fmap1_it->second).begin();
      farray2_it = (fmap2_it->second).begin();
      farray3_it = (fmap3_it->second).begin();
      for (
        ;
        farray1_it != (fmap1_it->second).end();
        farray1_it++, farray2_it++, farray3_it++
      ) {
        feature_vector.push_back(*farray1_it - 2 * *farray2_it + *farray3_it);
      }
    }
    // append feature_vector to the feature matrix
    feature_matrix.push_back(feature_vector);
  }
  return feature_matrix;
}

////
//// class OutlierCountAggregator
////
const std::string OutlierCountAggregator::name_ = "OutlierCountAggregator";

const std::string& OutlierCountAggregator::name() const {
  return name_;
}

feature_type OutlierCountAggregator::operator()(
  const feature_arrays& features
) {
  std::vector<size_t> outlier_ids = outlier_calculator_->calculate(features);
  return static_cast<feature_type>(outlier_ids.size());
}

////
//// class OutlierBadnessAggregator
////
const std::string OutlierBadnessAggregator::name_ = "OutlierBadnessAggregator";

const std::string& OutlierBadnessAggregator::name() const {
  return name_;
}

feature_type OutlierBadnessAggregator::operator()(
  const feature_arrays& features
) {
  outlier_calculator_->calculate(features);
  feature_array outlier_badness = outlier_calculator_->get_measures();
  feature_type ret = 0.0;
  for (
    feature_array::const_iterator f_it = outlier_badness.begin();
    f_it != outlier_badness.end();
    f_it++
  ) {
    ret = *f_it > ret ? *f_it : ret;
  }
  return ret;
}


////
//// function to_arma_matrix
////
arma::Mat<feature_type> to_arma_matrix(const feature_arrays& features) {
  assert(features.empty() == false);
  size_t cols = features.size();
  size_t rows = features[0].size();
  arma::Mat<feature_type> ret(rows, cols);

  typename std::vector<feature_array>::const_iterator feature_array_it;
  feature_array_it = features.begin();
  for (
    size_t j=0;
    feature_array_it != features.end();
    feature_array_it++, j++
  ) {
    assert(feature_array_it->size() == rows);
    arma::Col<feature_type> column(*feature_array_it);
    ret.col(j) = column;
  }
  return ret;
}

////
//// class OutlierCalculator
////
const std::string OutlierCalculator::name_ = "OutlierCalculator";

const std::string& OutlierCalculator::name() const {
  return name_;
}

const feature_array& OutlierCalculator::get_measures() const {
  throw std::runtime_error(
    "OutlierCalculator \"" + name() + "\"doesn't provide a measure"
  );
  return measures_;
}

////
//// class MVNOutlierCalculator
////
const std::string MVNOutlierCalculator::name_ = "MVNOutlierCalculator";

MVNOutlierCalculator::MVNOutlierCalculator(const feature_type sigma_threshold) {
  sigma_threshold_ = sigma_threshold;
}

const std::string& MVNOutlierCalculator::name() const {
  return name_;
}

const feature_array& MVNOutlierCalculator::get_measures() const {
  return measures_;
}

const arma::Mat<feature_type>& MVNOutlierCalculator::get_covariance() const {
  return covariance_;
}

const arma::Mat<feature_type>& MVNOutlierCalculator::get_inverse_covariance() const {
  return inv_covariance_;
}

const arma::Col<feature_type>& MVNOutlierCalculator::get_mean() const {
  return mean_;
}

const std::vector<size_t>& MVNOutlierCalculator::calculate(
  const feature_arrays& features
) {
  if (features.size() <= features[0].size()) {
    measures_.clear();
    outlier_ids_.clear();
    mean_.clear();
    covariance_.clear();
    inv_covariance_.clear();
  } else {
    // Get covariance and inverse covariance matrix
    arma::Mat<feature_type> features_mat(to_arma_matrix(features));
    arma::Mat<feature_type> features_mat_t(trans(features_mat));
    try {
      covariance_ = arma::cov(features_mat_t);
      inv_covariance_ = arma::inv_sympd(covariance_);

      LOG(logDEBUG4) << "In MVNOutlierCalculator: covariance matrix";
      LOG(logDEBUG4) << covariance_;
      LOG(logDEBUG4) << "In MVNOutlierCalculator: inverse covariance matrix";
      LOG(logDEBUG4) << inv_covariance_;

      // Get mean values
      mean_ = arma::mean(features_mat, 1);
  
      // Calculate the outliers
      outlier_ids_.clear();
      measures_.clear();
      feature_arrays::const_iterator features_it = features.begin();
      for(size_t id=0; features_it != features.end(); features_it++, id++) {
        arma::Col<feature_type> diff_vector(*features_it);
        diff_vector -= mean_;
        feature_type norm_residual = arma::dot(diff_vector, inv_covariance_*diff_vector);
        measures_.push_back(norm_residual);
        if (norm_residual >= sigma_threshold_) {
          outlier_ids_.push_back(id);
        }
      }
    } catch (std::exception& exception) {
      LOG(logDEBUG) << "In MVNOutlierCalculator: Too few data to calculate outliers" << std::endl;
    }
  } // else
  return outlier_ids_;
}

} // namespace pgmlink
