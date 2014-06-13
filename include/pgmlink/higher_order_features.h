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
#include <boost/shared_ptr.hpp> /* for shared_ptr */

// vigra
#include <vigra/multi_array.hxx> /* for the feature extractors */

namespace pgmlink {

/*=============================================================================
 type definitions
=============================================================================*/
typedef std::vector<HypothesesGraph::Node> Nodevector;
typedef std::vector<const Traxel*> ConstTraxelRefVector;

typedef feature_type FeatureScalar;
typedef vigra::MultiArray<1, feature_type> FeatureVector;
typedef vigra::MultiArray<2, feature_type> FeatureMatrix;

typedef vigra::MultiArrayView<1, feature_type> FeatureVectorView;
typedef vigra::MultiArrayView<2, feature_type> FeatureMatrixView;

/*=============================================================================
 functions
=============================================================================*/
void set_solution(HypothesesGraph& graph, const size_t solution_index);

/*=============================================================================
  pure virtual classes
=============================================================================*/
////
//// class TraxelsOfInterest
////
class TraxelsOfInterest {
 public:
  TraxelsOfInterest() {};
  virtual ~TraxelsOfInterest() {};
  virtual const std::string& name() const = 0;
  virtual const std::vector<ConstTraxelRefVector>& operator()(
    const HypothesesGraph& graph
  ) = 0;
};

////
//// class TraxelsFeatureExtractor
////
class TraxelsFeatureExtractor {
 public:
  TraxelsFeatureExtractor() {};
  virtual ~TraxelsFeatureExtractor() {};
  virtual const std::string& name() const = 0;
  virtual void extract(
    const ConstTraxelRefVector& traxelrefs,
    FeatureMatrix& feature_matrix
  ) const = 0;
  virtual FeatureMatrix extract(
    const ConstTraxelRefVector& traxelrefs
  ) const;
};

////
//// class TraxelsFeatureCalculator
////
class TraxelsFeatureCalculator {
 public:
  TraxelsFeatureCalculator() {};
  virtual ~ TraxelsFeatureCalculator() {};
  virtual const std::string& name() const = 0;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const = 0;
  virtual FeatureMatrix calculate(
    const FeatureMatrix& feature_matrix
  ) const;
};

/*=============================================================================
  specific classes
=============================================================================*/
////
//// class GraphFeatureCalculator
////
class GraphFeatureCalculator {
 public:
  GraphFeatureCalculator(
    boost::shared_ptr<TraxelsOfInterest> subsets_extractor_ptr,
    boost::shared_ptr<TraxelsFeatureExtractor> feature_extractor_ptr,
    boost::shared_ptr<TraxelsFeatureCalculator> feature_calculator_ptr
  ) : 
    subsets_extractor_ptr_(subsets_extractor_ptr),
    feature_extractor_ptr_(feature_extractor_ptr),
    feature_calculator_ptr_(feature_calculator_ptr) {};
  virtual ~GraphFeatureCalculator() {}
  virtual const FeatureVector& calculate_vector(
    const HypothesesGraph& graph
  );
 protected:
  boost::shared_ptr<TraxelsOfInterest> subsets_extractor_ptr_;
  boost::shared_ptr<TraxelsFeatureExtractor> feature_extractor_ptr_;
  boost::shared_ptr<TraxelsFeatureCalculator> feature_calculator_ptr_;
  FeatureVector ret_vector_;
};

////
//// class TraxelsFeaturesIdentity
////
class TraxelsFeaturesIdentity : public TraxelsFeatureExtractor {
 public:
  TraxelsFeaturesIdentity(const std::vector<std::string>& feature_names);
  TraxelsFeaturesIdentity(const std::string& feature_name);
  virtual ~TraxelsFeaturesIdentity() {};
  virtual const std::string& name() const;
  virtual void extract(
    const ConstTraxelRefVector& traxelrefs,
    FeatureMatrix& feature_matrix
  ) const;
 protected:
  static const std::string name_;
  std::vector<std::string> feature_names_;
}; // class TraxelsFeaturesIdentity

////
//// class TrackTraxels
////
class TrackTraxels : public TraxelsOfInterest {
 public:
  TrackTraxels() {};
  virtual ~TrackTraxels() {};
  virtual const std::string& name() const;
  virtual const std::vector<ConstTraxelRefVector>& operator()(
    const HypothesesGraph& graph
  );
 protected:
  static const std::string name_;
  std::vector<ConstTraxelRefVector> ret_;
};

////
//// class DivisionTraxels
////
class DivisionTraxels : public TraxelsOfInterest {
 public:
  DivisionTraxels(size_t depth = 1) : depth_(depth) {};
  virtual ~DivisionTraxels() {};
  virtual const std::string& name() const;
  virtual const std::vector<ConstTraxelRefVector>& operator()(
    const HypothesesGraph& graph
  );
  virtual const std::vector<ConstTraxelRefVector>& operator()(
    const HypothesesGraph& graph,
    size_t depth
  );
 protected:
  const std::vector<ConstTraxelRefVector>& from_tracklet_graph(
    const HypothesesGraph& graph,
    size_t depth
  );
  const std::vector<ConstTraxelRefVector>& from_traxel_graph(
    const HypothesesGraph& graph,
    size_t depth
  );
  bool get_children_to_depth(
    const HypothesesGraph::Node& node,
    const HypothesesGraph& graph,
    size_t depth,
    ConstTraxelRefVector& traxelrefs
  );
  bool get_parents_to_depth(
    const HypothesesGraph::Node& node,
    const HypothesesGraph& graph,
    size_t depth,
    ConstTraxelRefVector& traxelrefs
  );
  static const std::string name_;
  std::vector<ConstTraxelRefVector> ret_;
  size_t depth_;
};

////
//// class CompositionCalculator
////
class CompositionCalculator : public TraxelsFeatureCalculator {
 public:
  CompositionCalculator(
    boost::shared_ptr<TraxelsFeatureCalculator> first_calculator_ptr,
    boost::shared_ptr<TraxelsFeatureCalculator> second_calculator_ptr
  ) :
    first_calculator_ptr_(first_calculator_ptr),
    second_calculator_ptr_(second_calculator_ptr) {};
  virtual ~CompositionCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  static const std::string name_;
  boost::shared_ptr<TraxelsFeatureCalculator> first_calculator_ptr_;
  boost::shared_ptr<TraxelsFeatureCalculator> second_calculator_ptr_;
};

////
//// class SumCalculator
////
class SumCalculator : public TraxelsFeatureCalculator {
 public:
  SumCalculator() {};
  virtual ~SumCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_vector
  ) const;
 protected:
  static const std::string name_;
};

////
//// class DiffCalculator
////
class DiffCalculator : public TraxelsFeatureCalculator {
 public:
  DiffCalculator() {};
  virtual ~DiffCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  static const std::string name_;
};

////
//// class CurveCalculator
////
class CurveCalculator : public TraxelsFeatureCalculator {
 public:
  CurveCalculator() {};
  virtual ~CurveCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  static const std::string name_;
};

////
//// class MinCalculator
////
class MinCalculator : public TraxelsFeatureCalculator {
 public:
  MinCalculator() {};
  virtual ~MinCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  static const std::string name_;
};

////
//// class MaxCalculator
////
class MaxCalculator : public TraxelsFeatureCalculator {
 public:
  MaxCalculator() {};
  virtual ~MaxCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  static const std::string name_;
};

////
//// class MeanCalculator
////
class MeanCalculator : public TraxelsFeatureCalculator {
 public:
  MeanCalculator() {};
  virtual ~MeanCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  static const std::string name_;
};

////
//// class SquaredNormCalculator
////
class SquaredNormCalculator : public TraxelsFeatureCalculator {
 public:
  SquaredNormCalculator() {};
  virtual ~SquaredNormCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  static const std::string name_;
};

////
//// class SquaredDiffCalculator
////
class SquaredDiffCalculator : public TraxelsFeatureCalculator {
 public:
  SquaredDiffCalculator() {};
  virtual ~SquaredDiffCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  DiffCalculator diff_calculator_;
  SquaredNormCalculator squared_norm_calculator_;
  static const std::string name_;
};

////
//// class DiffusionCalculator
////
class DiffusionCalculator : public TraxelsFeatureCalculator {
 public:
  DiffusionCalculator() {};
  virtual ~DiffusionCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  SquaredDiffCalculator squared_diff_calculator_;
  MeanCalculator mean_calculator_;
  static const std::string name_;
};

////
//// class ChildParentDiffCalculator
////
class ChildParentDiffCalculator : public TraxelsFeatureCalculator {
 public:
  ChildParentDiffCalculator() {};
  virtual ~ChildParentDiffCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix,
    size_t depth
  ) const;
 protected:
  static const std::string name_;
};

////
//// class DotProductCalculator
////
class DotProductCalculator : public TraxelsFeatureCalculator {
 public:
  DotProductCalculator() {};
  virtual ~DotProductCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  static const std::string name_;
};

////
//// class ChildDeceleration
////
class ChildDeceleration : public TraxelsFeatureCalculator {
 public:
  ChildDeceleration() {};
  virtual ~ChildDeceleration() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix,
    size_t depth
  ) const;
 protected:
  static const std::string name_;
};

////
//// class MVNOutlierCalculator
////
class MVNOutlierCalculator : public TraxelsFeatureCalculator {
 public:
  MVNOutlierCalculator() {};
  virtual ~MVNOutlierCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix,
    const FeatureScalar& sigma_threshold
  ) const;
  virtual void calculate_inverse_covariance_matrix(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
  virtual void calculate_outlier_badness(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  MeanCalculator mean_calculator_;
  static const std::string name_;
};

} // namespace pgmlink

#endif // PGMLINK_HIGHER_ORDER_FEATURES_H
