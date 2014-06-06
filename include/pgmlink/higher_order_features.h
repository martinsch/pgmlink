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
//// class SubsetsOfInterest
////
class SubsetsOfInterest {
 public:
  SubsetsOfInterest() {};
  virtual ~SubsetsOfInterest() {};
  virtual const std::string& name() const = 0;
  virtual const std::vector<ConstTraxelRefVector>& operator()(
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
    const ConstTraxelRefVector& traxelrefs
  );
  virtual const FeatureVector& extract_vector(
    const ConstTraxelRefVector& traxelrefs
  );
  virtual FeatureScalar extract_scalar(
    const ConstTraxelRefVector& traxelrefs
  );
};

////
//// class SubsetFeatureCalculator
////
class SubsetFeatureCalculator {
 public:
  SubsetFeatureCalculator() {};
  virtual ~ SubsetFeatureCalculator() {};
  virtual const std::string& name() const = 0;
  virtual const FeatureMatrix& calculate_matrix(
    const FeatureMatrix& feature_matrix
  );
  virtual const FeatureVector& calculate_vector(
    const FeatureMatrix& feature_matrix
  );
  virtual FeatureScalar calculate_scalar(
    const FeatureMatrix& feature_matrix
  );
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
    boost::shared_ptr<SubsetsOfInterest> subsets_extractor_ptr,
    boost::shared_ptr<SubsetFeatureExtractor> feature_extractor_ptr,
    boost::shared_ptr<SubsetFeatureCalculator> feature_calculator_ptr
  ) : 
    subsets_extractor_ptr_(subsets_extractor_ptr),
    feature_extractor_ptr_(feature_extractor_ptr),
    feature_calculator_ptr_(feature_calculator_ptr) {};
  virtual ~GraphFeatureCalculator() {}
  virtual const FeatureVector& calculate_vector(
    const HypothesesGraph& graph
  );
 protected:
  boost::shared_ptr<SubsetsOfInterest> subsets_extractor_ptr_;
  boost::shared_ptr<SubsetFeatureExtractor> feature_extractor_ptr_;
  boost::shared_ptr<SubsetFeatureCalculator> feature_calculator_ptr_;
  FeatureVector ret_vector_;
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
    const ConstTraxelRefVector& traxelrefs
  );
 protected:
  static const std::string name_;
  std::vector<std::string> feature_names_;
  FeatureMatrix ret_matrix_;
}; // class SubsetFeaturesIdentity

////
//// class TrackSubsets
////
class TrackSubsets : public SubsetsOfInterest {
 public:
  TrackSubsets() {};
  virtual ~TrackSubsets() {};
  virtual const std::string& name() const;
  virtual const std::vector<ConstTraxelRefVector>& operator()(
    const HypothesesGraph& graph
  );
 protected:
  static const std::string name_;
  std::vector<ConstTraxelRefVector> ret_;
};

////
//// class DivisionSubsets
////
class DivisionSubsets : public SubsetsOfInterest {
 public:
  DivisionSubsets(size_t depth = 1) : depth_(depth) {};
  virtual ~DivisionSubsets() {};
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
class CompositionCalculator : public SubsetFeatureCalculator {
 public:
  CompositionCalculator(
    boost::shared_ptr<SubsetFeatureCalculator> first_calculator_ptr,
    boost::shared_ptr<SubsetFeatureCalculator> second_calculator_ptr
  ) :
    first_calculator_ptr_(first_calculator_ptr),
    second_calculator_ptr_(second_calculator_ptr) {};
  virtual ~CompositionCalculator() {};
  virtual const std::string& name() const;
  virtual const FeatureMatrix& calculate_matrix(
    const FeatureMatrix& feature_matrix
  );
  virtual const FeatureVector& calculate_vector(
    const FeatureMatrix& feature_matrix
  );
  FeatureScalar calculate_scalar(
    const FeatureMatrix& feature_matrix
  );
 protected:
  static const std::string name_;
  boost::shared_ptr<SubsetFeatureCalculator> first_calculator_ptr_;
  boost::shared_ptr<SubsetFeatureCalculator> second_calculator_ptr_;
};

////
//// class SumCalculator
////
class SumCalculator : public SubsetFeatureCalculator {
 public:
  SumCalculator() {};
  virtual ~SumCalculator() {};
  virtual const std::string& name() const;
  virtual const FeatureVector& calculate_vector(
    const FeatureMatrix& feature_matrix
  );
  virtual FeatureScalar calculate_scalar(
    const FeatureMatrix& feature_matrix
  );
 protected:
  static const std::string name_;
  FeatureVector ret_;
};

////
//// class DiffCalculator
////
class DiffCalculator : public SubsetFeatureCalculator {
 public:
  DiffCalculator() {};
  virtual ~DiffCalculator() {};
  virtual const std::string& name() const;
  virtual const FeatureMatrix& calculate_matrix(
    const FeatureMatrix& feature_matrix
  );
 protected:
  static const std::string name_;
  FeatureMatrix ret_;
};

////
//// class CurveCalculator
////
class CurveCalculator : public SubsetFeatureCalculator {
 public:
  CurveCalculator() {};
  virtual ~CurveCalculator() {};
  virtual const std::string& name() const;
  virtual const FeatureMatrix& calculate_matrix(
    const FeatureMatrix& feature_matrix
  );
 protected:
  static const std::string name_;
  FeatureMatrix ret_;
};

////
//// class MinCalculator
////
class MinCalculator : public SubsetFeatureCalculator {
 public:
  MinCalculator() {};
  virtual ~MinCalculator() {};
  virtual const std::string& name() const;
  virtual const FeatureVector& calculate_vector(
    const FeatureMatrix& feature_matrix
  );
  FeatureScalar calculate_scalar(
    const FeatureMatrix& feature_matrix
  );
 protected:
  static const std::string name_;
  FeatureVector ret_;
};

////
//// class MaxCalculator
////
class MaxCalculator : public SubsetFeatureCalculator {
 public:
  MaxCalculator() {};
  virtual ~MaxCalculator() {};
  virtual const std::string& name() const;
  virtual const FeatureVector& calculate_vector(
    const FeatureMatrix& feature_matrix
  );
  FeatureScalar calculate_scalar(
    const FeatureMatrix& feature_matrix
  );
 protected:
  static const std::string name_;
  FeatureVector ret_;
};

////
//// class MeanCalculator
////
class MeanCalculator : public SubsetFeatureCalculator {
 public:
  MeanCalculator() {};
  virtual ~MeanCalculator() {};
  virtual const std::string& name() const;
  virtual const FeatureVector& calculate_vector(
    const FeatureMatrix& feature_matrix
  );
  FeatureScalar calculate_scalar(
    const FeatureMatrix& feature_matrix
  );
 protected:
  static const std::string name_;
  FeatureVector ret_;
};

////
//// class SquaredNormCalculator
////
class SquaredNormCalculator : public SubsetFeatureCalculator {
 public:
  SquaredNormCalculator() {};
  virtual ~SquaredNormCalculator() {};
  virtual const std::string& name() const;
  virtual const FeatureVector& calculate_vector(
    const FeatureMatrix& feature_matrix
  );
 protected:
  static const std::string name_;
  FeatureVector ret_;
};

////
//// class SquaredDiffCalculator
////
class SquaredDiffCalculator : public SubsetFeatureCalculator {
 public:
  SquaredDiffCalculator() {};
  virtual ~SquaredDiffCalculator() {};
  virtual const std::string& name() const;
  virtual const FeatureMatrix& calculate_matrix(
    const FeatureMatrix& feature_matrix
  );
  virtual const FeatureVector& calculate_vector(
    const FeatureMatrix& feature_matrix
  );
 protected:
  DiffCalculator diff_calculator_;
  SquaredNormCalculator squared_norm_calculator_;
  static const std::string name_;
  FeatureMatrix ret_matrix_;
  FeatureVector ret_vector_;
};

////
//// class DiffusionCalculator
////
class DiffusionCalculator : public SubsetFeatureCalculator {
 public:
  DiffusionCalculator() {};
  virtual ~DiffusionCalculator() {};
  virtual const std::string& name() const;
  virtual FeatureScalar calculate_scalar(
    const FeatureMatrix& feature_matrix
  );
 protected:
  SquaredDiffCalculator squared_diff_calculator_;
  static const std::string name_;
};

////
//// class ChildParentDiffCalculator
////
class ChildParentDiffCalculator : public SubsetFeatureCalculator {
 public:
  ChildParentDiffCalculator() {};
  virtual ~ChildParentDiffCalculator() {};
  virtual const std::string& name() const;
  virtual const FeatureMatrix& calculate_matrix(
    const FeatureMatrix& feature_matrix
  );
  virtual const FeatureMatrix& calculate_matrix(
    const FeatureMatrix& feature_matrix,
    size_t division_depth
  );
 protected:
  static const std::string name_;
  FeatureMatrix ret_;
};

////
//// class DotProductCalculator
////
class DotProductCalculator : public SubsetFeatureCalculator {
 public:
  DotProductCalculator() {};
  virtual ~DotProductCalculator() {};
  virtual const std::string& name() const;
  virtual FeatureScalar calculate_scalar(
    const FeatureMatrix& feature_matrix
  );
 protected:
  static const std::string name_;
};

////
//// class ChildDeceleration
////
class ChildDeceleration : public SubsetFeatureCalculator {
 public:
  ChildDeceleration() {};
  virtual ~ChildDeceleration() {};
  virtual const std::string& name() const;
  virtual const FeatureVector& calculate_vector(
    const FeatureMatrix& feature_matrix
  );
  virtual const FeatureVector& calculate_vector(
    const FeatureMatrix& feature_matrix,
    size_t depth
  );
  virtual FeatureScalar calculate_scalar(
    const FeatureMatrix& feature_matrix
  );
  virtual FeatureScalar calculate_scalar(
    const FeatureMatrix& feature_matrix,
    size_t depth
  );
 protected:
  static const std::string name_;
  FeatureVector ret_vector_;
};

////
//// class MVNOutlierCalculator
////
class MVNOutlierCalculator : public SubsetFeatureCalculator {
 public:
  MVNOutlierCalculator() {};
  virtual ~MVNOutlierCalculator() {};
  virtual const std::string& name() const;
  virtual const FeatureMatrix& calculate_matrix(
    const FeatureMatrix& feature_matrix
  );
  virtual const FeatureVector& calculate_vector(
    const FeatureMatrix& feature_matrix
  );
  virtual FeatureScalar calculate_scalar(
    const FeatureMatrix& feature_matrix
  );
  virtual FeatureScalar calculate_scalar(
    const FeatureMatrix& feature_matrix,
    FeatureScalar sigma_threshold
  );
 protected:
  static const std::string name_;
  FeatureMatrix ret_matrix_;
  FeatureVector ret_vector_;
};

} // namespace pgmlink

#endif // PGMLINK_HIGHER_ORDER_FEATURES_H
