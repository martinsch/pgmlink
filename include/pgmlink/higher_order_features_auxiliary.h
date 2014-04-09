#ifndef PGMLINK_HIGHER_ORDER_FEATURES_AUXILIARY_H
#define PGMLINK_HIGHER_ORDER_FEATURES_AUXILIARY_H

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

} // namespace pgmlink

#endif
