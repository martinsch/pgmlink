#include "pgmlink/outlier_detection.h"

namespace pgmlink {
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
