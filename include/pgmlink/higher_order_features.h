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

/*=============================================================================
 type definitions
=============================================================================*/
typedef std::vector<HypothesesGraph::Node> Nodevector;
typedef std::vector<const Traxel*> ConstTraxelRefVector;

typedef feature_type FeatureScalar;
typedef vigra::MultiArray<1, feature_type> FeatureVector;
typedef vigra::MultiArray<2, feature_type> FeatureMatrix;

typedef vigra::MultiArrayView<1, feature_type> FeatureVectorView;

/*=============================================================================
 functions
=============================================================================*/
void set_solution(HypothesesGraph& graph, const size_t solution_index);

/*=============================================================================
  pure virtual classes
=============================================================================*/
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

/*
TODO clean file from needless code
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
*/

/*=============================================================================
  specific classes
=============================================================================*/
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

/*
////
//// class MVNOutlierCalculator
////
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
*/

} // namespace pgmlink

#endif // PGMLINK_HIGHER_ORDER_FEATURES_H
