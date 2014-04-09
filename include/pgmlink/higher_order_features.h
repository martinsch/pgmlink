/**
* @file
* @ingroup tracking
* @brief calculation of higher order features
*
* This file provides an interface to calculate higher order features from a
* tracking.
*/

#ifndef PGMLINK_HIGHER_ORDER_FEATURES_H
#define PGMLINK_HIGHER_ORDER_FEATURES_H

// stl
#include <vector>

// pgmlink
#include "pgmlink/traxels.h" /* for traxels */
#include "pgmlink/hypotheses.h" /* for hypotheses graph */

// boost
#include <boost/serialization/serialization.hpp> /* for serialization */

// armadillo
#include <armadillo> /* for MVNOutlierCalculator */

namespace pgmlink {

// forward declaration of class Track
class Track;

/*=============================================================================
 type definitions
=============================================================================*/
typedef std::vector<Traxel> Traxelvector;
typedef std::vector<Track> Trackvector;

/*=============================================================================
 class definitions
=============================================================================*/
////
//// class Track
////
class Track {
 public:
  Track();
  Track(
    const size_t id,
    const size_t time_start,
    const Traxelvector& traxels,
    const size_t parent_id = 0,
    const size_t left_child_id = 0,
    const size_t right_child_id = 0
  );
  void set_id(const size_t id);
  void set_time_start(const size_t time_start);
  void set_parent_id(const size_t parent_id);
  void set_child_ids(const size_t left, const size_t right);
  size_t get_id() const;
  size_t get_time_start() const;
  size_t get_parent_id() const;
  const std::vector<size_t>& get_child_ids() const;
  size_t get_length() const;

  size_t id_;
  size_t time_start_;
  Traxelvector traxels_;
  size_t parent_id_;
  std::vector<size_t> child_ids_;
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & archive, const unsigned int version) {
    (void) version;
    archive & id_;
    archive & time_start_;
    archive & traxels_;
    archive & parent_id_;
    archive & child_ids_;
  }
}; // class Track

////
//// class Tracking
////
class Tracking{
 public:
  Tracking();
  Tracking(const HypothesesGraph& graph, const size_t index);

  Trackvector tracks_;
  size_t index_;
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & archive, const unsigned int version) {
    (void) version;
    archive & tracks_;
    archive & index_;
  }
}; // class Tracking

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
  virtual feature_arrays operator()(const Track& track) = 0;
};

////
//// class FeatureAggregator
////
class FeatureAggregator {
 public:
  FeatureAggregator() {};
  virtual ~FeatureAggregator() {};
  virtual const std::string& name() const = 0;
  virtual feature_type operator()(const feature_arrays& features) = 0;
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
  virtual feature_arrays operator()(const Track& track);
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
  virtual feature_arrays operator()(const Track& track);
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
  virtual feature_arrays operator()(const Track& track);
 protected:
  static const std::string name_;
  std::vector<std::string> feature_names_;
};

////
//// class MVNOutlierCalculator
////
/* the template parameter "sigma_threshold" will be scaled with factor 1/1000.
A sigma_threshold of 3000 corresponds to an actual threshold of 3.000 */
class MVNOutlierCalculator : public OutlierCalculator {
  public:
    MVNOutlierCalculator(const feature_type sigma_threshold = 3.0);
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

////
//// class OutlierCountAggregator
////
class OutlierCountAggregator : public FeatureAggregator {
 public:
  OutlierCountAggregator(
    OutlierCalculator* outlier_calculator = new MVNOutlierCalculator()
  ) : outlier_calculator_(outlier_calculator) {};
  virtual ~OutlierCountAggregator() {};
  virtual const std::string& name() const;
  virtual feature_type operator()(
    const feature_arrays& features
  );
 protected:
  static const std::string name_;
  OutlierCalculator* outlier_calculator_;
};

////
//// class OutlierBadnessAggregator
////
class OutlierBadnessAggregator : public FeatureAggregator {
 public:
  OutlierBadnessAggregator(
    OutlierCalculator* outlier_calculator = new MVNOutlierCalculator()
  ) : outlier_calculator_(outlier_calculator) {};
  virtual ~OutlierBadnessAggregator() {};
  virtual const std::string& name() const;
  virtual feature_type operator()(
    const feature_arrays& features
  );
 protected:
  static const std::string name_;
  OutlierCalculator* outlier_calculator_;
};

} // namespace pgmlink

#endif // PGMLINK_HIGHER_ORDER_FEATURES_H
