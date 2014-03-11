#include "pgmlink/outlier_detection.h"

namespace pgmlink {

//
// function to_arma_mat
//
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

//
// Class OutlierCalculator
//
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

//
// Class MVNOutlierCalculator
//
const std::string MVNOutlierCalculator::name_ = "MVNOutlierCalculator";

const feature_type MVNOutlierCalculator::sigma_threshold_ = 3.0;

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
    std::cerr << "Too few data to calculate outliers" << std::endl;
  } else {
    // Get covariance and inverse covariance matrix
    arma::Mat<feature_type> features_mat(to_arma_matrix(features));
    arma::Mat<feature_type> features_mat_t(trans(features_mat));
    covariance_ = arma::cov(features_mat_t);
    inv_covariance_ = arma::inv_sympd(covariance_);
  
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
  } // else
  return outlier_ids_;
}

//
// Class Outlier
//
Outlier::Outlier(
  boost::shared_ptr<TrackFeatureExtractor> extractor,
  boost::shared_ptr<OutlierCalculator> outlier_calculator
) : extractor_(extractor), outlier_calculator_(outlier_calculator) {
}

const std::vector<size_t> Outlier::calculate(const Track& track) {
  feature_arrays extracted_features = extractor_->extract(track);
  return outlier_calculator_->calculate(extracted_features);
}

const feature_array& Outlier::get_measures() const {
  return outlier_calculator_->get_measures();
}

} // Namespace pgmlink
