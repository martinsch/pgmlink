#ifndef PGMLINK_HIGHER_ORDER_FEATURES_AUXILIARY_H
#define PGMLINK_HIGHER_ORDER_FEATURES_AUXILIARY_H

#include <armadillo>

#include "pgmlink/higher_order_features.h" // for Track, Tracks, ...

namespace pgmlink{

/*=============================================================================
  pure virtual functors
=============================================================================*/

////
//// class TrackFeatureExtractor
////
class TrackFeatureExtractor {
 public:
  TrackFeatureExtractor() {};
  virtual ~TrackFeatureExtractor() {};
  virtual const std::string& name() const = 0;
  virtual const feature_arrays operator()(const Track& track) const = 0;
};

////
//// class FeatureAggregator
////
class FeatureAggregator {
 public:
  FeatureAggregator() {};
  virtual ~FeatureAggregator() {};
  virtual const std::string& name() const = 0;
  virtual feature_type operator()(
    const feature_arrays& features
  ) const = 0;
};

/*=============================================================================
  specific functors
=============================================================================*/
////
//// class TrackFeaturesIdentity
////
class TrackFeaturesIdentity : public TrackFeatureExtractor {
 public:
  TrackFeaturesIdentity(const std::vector<std::string>& feature_names);
  TrackFeaturesIdentity(const std::string& feature_name);
  virtual ~TrackFeaturesIdentity() {};
  virtual const std::string& name() const;
  virtual const feature_arrays operator()(const Track& track) const;
 protected:
  static const std::string name_;
  std::vector<std::string> feature_names_;
}; // class TrackFeaturesIdentity

////
//// class TrackFeaturesDiff
////
class TrackFeaturesDiff : public TrackFeatureExtractor {
 public:
  TrackFeaturesDiff(const std::vector<std::string>& feature_names);
  TrackFeaturesDiff(const std::string& feature_name);
  virtual ~TrackFeaturesDiff() {};
  virtual const std::string& name() const;
  virtual const feature_arrays operator()(const Track& track) const;
 protected:
  static const std::string name_;
  std::vector<std::string> feature_names_;
};

////
//// class TrackFeaturesCurvature
////
class TrackFeaturesCurvature : public TrackFeatureExtractor {
 public:
  TrackFeaturesCurvature(const std::vector<std::string>& feature_names);
  TrackFeaturesCurvature(const std::string& feature_name);
  virtual ~TrackFeaturesCurvature() {};
  virtual const std::string& name() const;
  virtual const feature_arrays operator()(const Track& track) const;
 protected:
  static const std::string name_;
  std::vector<std::string> feature_names_;
};

////
//// class OutlierCountAggregator
////
template<typename OutlierCalculator_T>
class OutlierCountAggregator : public FeatureAggregator {
 public:
  OutlierCountAggregator() : outlier_calculator_() {};
  virtual ~OutlierCountAggregator() {};
  virtual const std::string& name() const;
  virtual feature_type operator()(
    const feature_arrays& features
  ) const;
 protected:
  static const std::string name_;
  OutlierCalculator_T outlier_calculator_;
};

////
//// class OutlierBadnessAggregator
////
template<typename OutlierCalculator_T>
class OutlierBadnessAggregator : public FeatureAggregator {
 public:
  OutlierBadnessAggregator() : outlier_calculator_() {};
  virtual ~OutlierBadnessAggregator() {};
  virtual const std::string& name() const;
  virtual feature_type operator()(
    const feature_arrays& features
  ) const;
 protected:
  static const std::string name_;
  OutlierCalculator_T outlier_calculator_;
};

////
//// class OutlierCalculator
////
class OutlierCalculator {
  public:
    OutlierCalculator() {};
    ~OutlierCalculator() {};
    virtual const std::vector<size_t>& calculate(
      const feature_arrays& features
    ) = 0;
    virtual const feature_array& get_measures() const;
    virtual const std::string& name() const;
  protected:
    static const std::string name_;
    feature_array measures_;
}; // class OutlierCalculator

////
//// class MVNOutlierCalculator
////
/* the template parameter "sigma_threshold" will be scaled with factor 1/1000.
A sigma_threshold of 3000 corresponds to an actual threshold of 3.000 */
template<int sigma_threshold = 3000>
class MVNOutlierCalculator : public OutlierCalculator {
  public:
    MVNOutlierCalculator();
    ~MVNOutlierCalculator() {};
    const std::vector<size_t>& calculate(const feature_arrays& features);
    const feature_array& get_measures() const;
    const arma::Mat<feature_type>& get_covariance() const;
    const arma::Mat<feature_type>& get_inverse_covariance() const;
    const arma::Col<feature_type>& get_mean() const;
    virtual const std::string& name() const;
  protected:
    feature_type sigma_threshold_;
    static const std::string name_;
    feature_array measures_;
    std::vector<size_t> outlier_ids_;
    arma::Col<feature_type> mean_;
    arma::Mat<feature_type> covariance_;
    arma::Mat<feature_type> inv_covariance_;
}; // class MVNOutlierCalculator

} // namespace pgmlink

#endif
