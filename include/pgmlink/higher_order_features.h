/**
\file
\brief calculation of higher order features

This file provides an interface to calculate higher order features of a
tracking. Workflow of the calculation of higher order feautres of the MAP
solution:
- Use an implementation of <TT>pgmlink::TraxelsOfInterest</TT> to get all
interesting subsets of traxels, e.g. a track.
- Use an implementation of the <TT>pgmlink::TraxelsFeatureExtractor</TT> to
extract features like positions or sizes of the involved traxels.
- Compose the implementations of <TT>pgmlink::TraxelsFeatureCalculator</TT> to
calculate any interesting new higher order feature

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
/**
\brief write the solution stored in the property maps <TT>node_active_count</TT>
  and <TT>arc_active_count</TT> into the <TT>node_active</TT> and 
  <TT>arc_active</TT>.
*/
void set_solution(HypothesesGraph& graph, const size_t solution_index);

/*=============================================================================
  pure virtual classes
=============================================================================*/
////
//// class TraxelsOfInterest
////
/**
\brief The class <TT>TraxelsOfInterest</TT> extracts subsets of the graph that
  are of interest, e.g. tracks or division.

The <TT>TraxelsOfInterest</TT> is a pure virtual class. Its ()-operator returns
a vector of vectors of traxel references. One vector of traxel references
represents one subset of traxels of the whole graph over which the higher order
features will be calculated.

Interesting subsets of traxels might be the traxels in one track, or the traxels
involved in a cell division.
*/
class TraxelsOfInterest {
 public:
  TraxelsOfInterest() {};
  virtual ~TraxelsOfInterest() {};
  virtual const std::string& name() const = 0;
  /**
  \brief Extract the traxels of interest of the hypotheses graph.

  \param[in] graph Hypotheses graph with its MAP-solution written into the
    property maps <TT>node_active</TT> and <TT>arc_active</TT>.
  \return vector of vectors of traxel references: ((t1, t2, t3), (t1, t3), ...).
    They are stored as std::vector<std::vector<const *pgmlink::Traxel> >
  */
  virtual const std::vector<ConstTraxelRefVector>& operator()(
    const HypothesesGraph& graph
  ) = 0;
};

////
//// class TraxelsFeatureExtractor
////
/**
\brief virtual class for extracting the interesting features of the traxels.

*/
class TraxelsFeatureExtractor {
 public:
  TraxelsFeatureExtractor() {};
  virtual ~TraxelsFeatureExtractor() {};
  virtual const std::string& name() const = 0;
  /**
  \brief extracts any features from the traxelvector given as the argument
  
  \param[in] traxelrefs Vector of references to traxels of which the features
    should be extracted.
  \param[out] feature_matrix two dimensional matrix that contains the extracted
    features
  */
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
/**
\brief virtual class for the calculation of features of any order.

The calculations of the higher order features are implemented as child classes
of this virtual class.
*/
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
/**
\brief extract the features of the traxels with the feature names given in the
  constructor, e.g. "size"

In the traxels there are the property maps stored which map a string to a 
vector. This class extracts with the extract method the vectors with the feature
names given in the constructor as a matrix. The columns are the extracted
feature vectors for the traxels.
*/
class TraxelsFeaturesIdentity : public TraxelsFeatureExtractor {
 public:
  /**
  \brief Takes the feature name that should be extracted in the extract method.

  The extract method the extracts this feature vector of the feature map.
  */
  TraxelsFeaturesIdentity(const std::string& feature_name);
  /**
  \brief Takes many feature names that should be extracted in the extract
    method.
  */
  TraxelsFeaturesIdentity(const std::vector<std::string>& feature_names);
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
/**
\brief gets all tracks in the MAP solution. On element in the return vector
  is a set of traxels that belong to the same track.

A track is a sequence of one to one connected nodes. A track begins if a node
has no incoming arc, more than one incoming arcs or a parent with two outgoing
arcs. A track is terminated if there is no outgoing arc, more than
one outgoing arc or a child with two incoming arcs.
The tracks in the following example are:
(n1, n2), (n3, n4), (n5, n6), (n7, n8)
\code
 n1 - n2 - n3 - n4
         \
 n5 - n6 - n7 - n8
\endcode
*/
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
/**
\brief identifies all cell divisions in the tracking an returns the references
  to the involved traxels.

A cell division is detected if one node has two outgoing arcs. One can specify
the depth to which the traxels are extracted. That is to say how many parent
and child nodes are returned.
Division of depth 2:
\code
          n2 - n3
        /
n1 - n0
        \
          n4 - n5
\endcode
The numbers indicate their position in the return vector as well. The dividing
cell is always in the 0th position.
*/
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
template<int N>
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
template<int N>
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
template<int N>
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
template<int N>
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
//// class StdCalculator
////
template<int N>
class StdCalculator : public TraxelsFeatureCalculator {
 public:
  StdCalculator() {};
  virtual ~StdCalculator() {};
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
  MeanCalculator<0> mean_calculator_;
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
  MeanCalculator<0> mean_calculator_;
  static const std::string name_;
};

} // namespace pgmlink

#endif // PGMLINK_HIGHER_ORDER_FEATURES_H
