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

// vigra
#include <vigra/multi_array.hxx> /* for the feature extractors */

namespace pgmlink {

// forward declaration of class Track
class Track;

/*=============================================================================
 type definitions
=============================================================================*/
typedef std::vector<Traxel> Traxelvector;
typedef std::vector<Track> Trackvector;
typedef std::vector<HypothesesGraph::Node> Nodevector;

typedef feature_type FeatureScalar;
typedef vigra::MultiArray<1, feature_type> FeatureVector;
typedef vigra::MultiArray<2, feature_type> FeatureMatrix;

typedef vigra::MultiArrayView<1, feature_type> FeatureVectorView;

/*=============================================================================
 function definitions
=============================================================================*/
void set_solution(HypothesesGraph& graph);

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
    Track* const parent = NULL,
    Track* const left_child = NULL,
    Track* const right_child = NULL
  );
  void set_id(const size_t id);
  void set_time_start(const size_t time_start);
  void set_parent(Track* const parent);
  void set_child(Track* const child);
  void set_children(Track* const left_child, Track* const right_child);
  size_t get_id() const;
  size_t get_time_start() const;
  Track* get_parent() const;
  const std::vector<Track*>& get_children() const;
  size_t get_length() const;

  size_t id_;
  size_t time_start_;
  Traxelvector traxels_;
  Track* parent_;
  std::vector<Track*> children_;
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & archive, const unsigned int version) {
    (void) version;
    archive & id_;
    archive & time_start_;
    archive & traxels_;
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
  pure virtual classes
=============================================================================*/
//// TODO: deprecated / use FeatureExtractor
//// class TrackFeatureExtractor
////
class TrackFeatureExtractor {
 public:
  TrackFeatureExtractor() {};
  virtual ~TrackFeatureExtractor() {};
  virtual const std::string& name() const = 0;
  virtual const feature_arrays& operator()(const Track& track) = 0;
};

////
//// class FeatureAggregator
////
class FeatureAggregator {
 public:
  FeatureAggregator() {};
  virtual ~FeatureAggregator() {};
  virtual const std::string& name() const = 0;
  virtual const feature_type& operator()(const feature_arrays& features) = 0;
};

////
//// class SubsetsOfInterest
////
class SubsetsOfInterest {
 public:
  SubsetsOfInterest() {};
  virtual ~SubsetsOfInterest() {};
  virtual const std::string& name() const = 0;
  virtual const std::vector<Nodevector>& operator()(
    const HypothesesGraph& graph
  ) = 0;
};

////
//// class SubsetFeatureExtractor
////
class SubsetFeatureExtractor {
 public:
  SubsetFeatureExtractor() {};
  virtual ~SubsetFeatureExtractor() {};
  virtual const std::string& name() const = 0;
  virtual const FeatureMatrix& extract_matrix(
    const Nodevector& nodevector,
    const HypothesesGraph& graph
  );
  virtual const FeatureVector& extract_vector(
    const Nodevector& nodevector,
    const HypothesesGraph& graph
  );
  virtual const FeatureScalar& extract_scalar(
    const Nodevector& nodevector,
    const HypothesesGraph& graph
  );
};

////
//// class SubsetFeatureAggregator
////
class SubsetFeatureAggregator {
 public:
  SubsetFeatureAggregator() {};
  virtual ~SubsetFeatureAggregator() {};
  virtual const std::string& name() const = 0;
  virtual const feature_array& operator()(
    const feature_arrays& features
  ) = 0;
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
//// class TrackValue
////
class TrackValue {
 public:
  TrackValue(
    TrackFeatureExtractor* track_feature_extractor,
    FeatureAggregator* feature_aggregator
  ) :
    track_feature_extractor_(track_feature_extractor),
    feature_aggregator_(feature_aggregator)
  {};
  const feature_type& operator()(const Track& track);
 protected:
  TrackFeatureExtractor* track_feature_extractor_;
  FeatureAggregator* feature_aggregator_;
  feature_type ret_;
};

////
//// class TrackingValue
////
class TrackingValue {
 public:
  TrackingValue(
    SubsetsOfInterest* subsets_of_interest,
    SubsetFeatureExtractor* subset_feature_extractor,
    SubsetFeatureAggregator* subset_feature_aggregator,
    FeatureAggregator* feature_aggregator
  );
  const feature_type& operator()(const Tracking& tracking);
 protected:
  SubsetsOfInterest* subsets_of_interest_;
  SubsetFeatureExtractor* subset_feature_extractor_;
  SubsetFeatureAggregator* subset_feature_aggregator_;
  FeatureAggregator* feature_aggregator_;
  feature_type ret_;
};

////
//// class SubsetFeaturesIdentity
////
class SubsetFeaturesIdentity : public SubsetFeatureExtractor {
 public:
  SubsetFeaturesIdentity(const std::vector<std::string>& feature_names);
  SubsetFeaturesIdentity(const std::string& feature_name);
  virtual ~SubsetFeaturesIdentity() {};
  virtual const std::string& name() const;
  virtual const FeatureMatrix& extract_matrix(
    const Nodevector& nodevector,
    const HypothesesGraph& graph
  );
 protected:
  static const std::string name_;
  std::vector<std::string> feature_names_;
  FeatureMatrix ret_matrix_;
}; // class SubsetFeaturesIdentity

////
//// class TrackFeaturesDiff
////
class TrackFeaturesDiff : public TrackFeatureExtractor {
 public:
  TrackFeaturesDiff(const std::vector<std::string>& feature_names);
  TrackFeaturesDiff(const std::string& feature_name);
  virtual ~TrackFeaturesDiff() {};
  virtual const std::string& name() const;
  virtual const feature_arrays& operator()(const Track& track);
 protected:
  static const std::string name_;
  std::vector<std::string> feature_names_;
  feature_arrays ret_;
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
  virtual const feature_arrays& operator()(const Track& track);
 protected:
  static const std::string name_;
  std::vector<std::string> feature_names_;
  feature_arrays ret_;
};

////
//// class SumAggregator
////
class SumAggregator : public FeatureAggregator {
 public:
  SumAggregator() {};
  virtual ~SumAggregator() {};
  virtual const std::string& name() const;
  virtual const feature_type& operator()(const feature_arrays& features);
 protected:
  static const std::string name_;
  feature_type ret_;
};

////
//// class TrackSubsets
////
class TrackSubsets : public SubsetsOfInterest {
 public:
  TrackSubsets() {};
  virtual ~TrackSubsets() {};
  virtual const std::string& name() const;
  virtual const std::vector<Nodevector>& operator()(
    const HypothesesGraph& graph
  );
 protected:
  static const std::string name_;
  std::vector<Nodevector> ret_;
};

////
//// class DivisionSubsets
////
class DivisionSubsets : public SubsetsOfInterest {
 public:
  DivisionSubsets() {};
  virtual ~DivisionSubsets() {};
  virtual const std::string& name() const;
  virtual const std::vector<Nodevector>& operator()(
    const HypothesesGraph& graph
  );
 protected:
  static const std::string name_;
  std::vector<Nodevector> ret_;
};

////
//// class SubsetFeatureExtractorFromFE
////
class SubsetFeatureExtractorFromFE : public SubsetFeatureExtractor {
 public:
  SubsetFeatureExtractorFromFE(
    TrackFeatureExtractor* track_feature_extractor
  ) : track_feature_extractor_(track_feature_extractor) {};
  virtual ~SubsetFeatureExtractorFromFE() {};
  virtual const std::string& name() const;
  virtual const feature_arrays& operator()(
    const Trackvector& tracks
  );
 protected:
  static const std::string name_;
  TrackFeatureExtractor* track_feature_extractor_;
  feature_arrays ret_;
};

////
//// class DivisionFeatureExtractor
////
// class DivisionFeatureExtractor : public SubsetFeatureExtractor {
//  public:
//   DivisionFeatureExtractor(const std::string& feature_name, size_t depth=1);
//   DivisionFeatureExtractor(
//     const std::vector<std::string>& feature_names,
//     size_t depth=1
//   );
//   virtual ~DivisionFeatureExtractor() {};
//   virtual const std::string& name() const;
//   virtual const feature_arrays& operator()(
//     const Trackvector& tracks
//   );
//  protected:
//   static const std::string name_;
//   TrackFeaturesIdentity features_identity_;
//   size_t depth_;
//   feature_arrays ret_;
// };

////
//// class SubsetFeatureAggregatorFromFA
////
class SubsetAggregatorFromFA : public SubsetFeatureAggregator {
 public:
  SubsetAggregatorFromFA(
    FeatureAggregator* feature_aggregator
  ) : feature_aggregator_(feature_aggregator) {};
  virtual ~SubsetAggregatorFromFA() {};
  virtual const std::string& name() const;
  virtual const feature_array& operator()(const feature_arrays& features);
 protected:
  static const std::string name_;
  FeatureAggregator* feature_aggregator_;
  feature_array ret_;
};

////
//// class ChildRatioAggregator
////
class ChildRatioAggregator : public SubsetFeatureAggregator {
 public:
  ChildRatioAggregator(const size_t depth=1) : depth_(depth) {};
  virtual ~ChildRatioAggregator() {};
  virtual const std::string& name() const;
  virtual const feature_array& operator()(const feature_arrays& features);
 protected:
  static const std::string name_;
  size_t depth_;
  feature_array ret_;
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
  virtual const feature_type& operator()(
    const feature_arrays& features
  );
 protected:
  static const std::string name_;
  OutlierCalculator* outlier_calculator_;
  feature_type ret_;
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
  virtual const feature_type& operator()(
    const feature_arrays& features
  );
 protected:
  static const std::string name_;
  OutlierCalculator* outlier_calculator_;
  feature_type ret_;
};

} // namespace pgmlink

#endif // PGMLINK_HIGHER_ORDER_FEATURES_H
