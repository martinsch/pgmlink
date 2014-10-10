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

Minimal example:
\code
#include <pgmlink/higher_order_features.h>

using namespace pgmlink

// Get the hypotheses graph with the MAP solution written into the node_active
// and arc_active property map
HypothesesGraph graph;
...

// Extract all tracks in the tracking
TrackTraxels track_extractor;
std::vector<ConstTraxelRefVector> track_traxels = track_extractor(graph);

// Extract the positions of the traxels in the first track
TraxelsFeaturesIdentity position_extractor("com");
FeatureMatrix positions;
position_extractor.extract(track_traxels[0], positions);

// Calculate the diffusion coefficient of those positions
DiffusionCalculator diffusion_calculator;
FeatureMatrix return_value;
diffusion_calculator.calculate(positions, return_value);

// The return_value is now a 1x1 matrix with the diffusion coefficient in the
// (0,0) position
\endcode

*/

#ifndef PGMLINK_HIGHER_ORDER_FEATURES_H
#define PGMLINK_HIGHER_ORDER_FEATURES_H

// stl
#include <vector>

// pgmlink
#include "pgmlink/traxels.h" /* for traxels */
#include "pgmlink/hypotheses.h" /* for hypotheses graph */
#include "pgmlink/classifier_auxiliary.h" /* for class FeatureCalculator */

// boost
#include <boost/serialization/serialization.hpp> /* for serialization */
#include <boost/shared_ptr.hpp> /* for shared_ptr */

// vigra
#include <vigra/multi_array.hxx> /* for the feature extractors */

namespace pgmlink {
namespace features {

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
\brief virtual class to extract subsets of the graph that are of interest, e.g.
  tracks or division.

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
/**
\brief simplifies the workflow TraxelsOfInterest -> TraxelsFeatureExtractor ->
  TraxelsFeatureCalculator

\code
                          TraxelsOfInterest
                    Graph ----------------------> Vector of vectors of traxels

                          TraxelsFeatureExtractor
   i-th vector of Traxels ----------------------> matrix of traxel features

                          TraxelsFeatureCalculator
matrix of traxel features ----------------------> nxm matrix of traxel features
\endcode
The GraphFeatureCalculator uses the TraxelsFeatureExtractor and 
TraxelsFeatureCalculator given in the constructor to calculate the features 
for of every subset returned by the TraxelsOfInterest class in the constructor.
The return_matrix is a row vector of all calculated feature matrix entries.
*/
class GraphFeatureCalculator {
 public:
  /**
  \brief takes shared pointers to the TraxelsOfInterest, TraxelsFeatureExtractor
    and TraxelsFeatureCalculator that should be used in the calculation workflow
  
  Example usage:
  \code
  boost::shared_ptr<TraxelsOfInterest> track_extractor_ptr(
    new TrackExtractor
  );
  boost::shared_ptr<TraxelsFeatureExtractor> position_identity_ptr(
    new TraxelsFeaturesIdentity("com")
  );
  boost::shared_ptr<TraxelsFeatureCalculator> diffusion_calculator_ptr(
    new DiffusionCalculator
  );

  GraphFeatureCalculator diffusion_coefficients(
    track_extractor_ptr,
    position_identity_ptr,
    diffusion_calculator_ptr
  );
  
  FeatureMatrix return_matrix;
  diffusion_coefficients.calculate(graph, return_matrix);
  // the return matrix is now a row vector of all diffusion coefficients of all
  // tracks in the hypotheses graph.
  \endcode
  */
  GraphFeatureCalculator(
    boost::shared_ptr<TraxelsOfInterest> subsets_extractor_ptr,
    boost::shared_ptr<TraxelsFeatureExtractor> feature_extractor_ptr,
    boost::shared_ptr<TraxelsFeatureCalculator> feature_calculator_ptr
  ) : 
    subsets_extractor_ptr_(subsets_extractor_ptr),
    feature_extractor_ptr_(feature_extractor_ptr),
    feature_calculator_ptr_(feature_calculator_ptr) {};
  virtual ~GraphFeatureCalculator() {}
  virtual void calculate(
    const HypothesesGraph& graph,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  boost::shared_ptr<TraxelsOfInterest> subsets_extractor_ptr_;
  boost::shared_ptr<TraxelsFeatureExtractor> feature_extractor_ptr_;
  boost::shared_ptr<TraxelsFeatureCalculator> feature_calculator_ptr_;
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
/**
\brief Composition of two TraxelsFeatureCalculator
*/
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
//// class TCompositionCalculator
//// template version
/**
\brief templated version of the CompositionCalculator
*/
template<typename FC1, typename FC2>
class TCompositionCalculator : public TraxelsFeatureCalculator {
 public:
  TCompositionCalculator() {};
  virtual ~TCompositionCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  static const std::string name_;
  FC1 first_calculator_;
  FC2 second_calculator_;
};

template<typename FC1, typename FC2>
const std::string TCompositionCalculator<FC1, FC2>::name_ = "CompositionCalculator";

template<typename FC1, typename FC2>
const std::string& TCompositionCalculator<FC1, FC2>::name() const {
  return name_;
}

template<typename FC1, typename FC2>
void TCompositionCalculator<FC1, FC2>::calculate(
  const FeatureMatrix& feature_matrix,
  FeatureMatrix& return_matrix
) const {
  FeatureMatrix temp;
  first_calculator_.calculate(feature_matrix, temp);
  second_calculator_.calculate(temp, return_matrix);
}

////
//// class TraxelsFCFromFC
////
/**
\brief helper class to use the feature calculators from the classifier auxiliary
  with the higher order features interface.

The pgmlink::TraxelsFCFromFC is constructed with a pgmlink::FeatureCalculator.
The calculate method then calcultes the result matrix with this feature
calculator. The order parameter given in the constructor determines how many
neighbouring columns are used to calculate one column in the return matrix.

Let \f$f(\vec{x_1}, \vec{x_2}, \vec{x_3})\f$ be a feature calculator of order 3.
Then TraxelsFCFromFC::calculate method does the following:

\f[
  (\vec{x_1}, \vec{x_2}, \dots, \vec{x_n})
  \mapsto
  (f(\vec{x_1}, \vec{x_2}, \vec{x_3}),
  f(\vec{x_2}, \vec{x_3}, \vec{x_4}),
  \dots,
  f(\vec{x_{n-2}}, \vec{x_{n-1}}, \vec{x_n}))
\f]

If there are fewer column vectors in the input matrix than the order of the
feature calculator, this calculator returns a 1x1 matrix with a 0 in the (0,0)
position.
*/
class TraxelsFCFromFC : public TraxelsFeatureCalculator {
 public:
  TraxelsFCFromFC(
    boost::shared_ptr<FeatureCalculator> feature_calculator_ptr,
    size_t order
  ) : feature_calculator_ptr_(feature_calculator_ptr), order_(order) {};
  virtual ~TraxelsFCFromFC() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_vector
  ) const;
 protected:
  boost::shared_ptr<FeatureCalculator> feature_calculator_ptr_;
  size_t order_;
  static const std::string name_;
};

////
//// class SumCalculator
////
/**
\brief Calculates the sum along the matrix axis specified in the template
  argument

\tparam N axis along which the sum is taken

- N=0: The sum is taken along the rows, returns a column vector.
- N=1: The sum is taken along the columns, returns a row vector.
- N=-1: The sum is taken over all elements of the matrix, returns a 1x1 matrix.
*/
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
//// class DeviationCalculator
////
/**
\brief calculates the deviation of the column vectors to the mean column vector
*/
class DeviationCalculator : public TraxelsFeatureCalculator {
 public:
  DeviationCalculator() {};
  virtual ~DeviationCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  MeanCalculator<0> mean_calculator_;
  static const std::string name_;
};

/**
\brief the square of each element
*/
class SquareCalculator : public TraxelsFeatureCalculator {
 public:
  SquareCalculator() {};
  virtual ~SquareCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  static const std::string name_;
};

/**
\brief the square root of each element
*/
class SquareRootCalculator : public TraxelsFeatureCalculator {
 public:
  SquareRootCalculator() {};
  virtual ~SquareRootCalculator() {};
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
/**
\brief calculates the squared norm of the row or column vectors

\tparam N axis along which the norm is taken

- N=0 Calculates the squared norm of each column vector
- N=1 Calculates the squared norm of each row vector
*/
template<int N>
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

/**
\brief calculates the euclidean norm for each column vector
*/
typedef TCompositionCalculator<
  SquaredNormCalculator<0>,
  SquareRootCalculator
> EuclideanNormCalculator;


/**
\brief calculates the variance for each row.
*/
typedef TCompositionCalculator<
  TCompositionCalculator<
    DeviationCalculator,
    SquareCalculator
  >,
  MeanCalculator<0>
> VarianceCalculator;

/**
\brief calculates the squared norm of the difference of neighbouring column
  vectors of the
*/
typedef TCompositionCalculator<
  DiffCalculator,
  SquaredNormCalculator<0>
> SquaredDiffCalculator;

////
//// class DiffusionCalculator
////
/**
\brief calculates the diffusion coefficient of the column vectors of the input
  matrix.

Returns \f$D = \langle \vec{v}^2 \rangle - {\langle \vec{v} \rangle}^2\f$, where
\f$\vec{v}\f$ is the difference of two neigbouring column vectors. It is
therefore variance of the mean move distance if the column vectors are the cell
positions.
*/
typedef TCompositionCalculator<
  SquaredDiffCalculator,
  MeanCalculator<0>
> DiffusionCalculator;

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
//// class AngleCosineCalculator
////
/**
\brief calculates the cosine of the angle of three neighbouring column vectors

AngleCosineCalculator returns the cosine of an angle. This angle is the change
of direction of three neighbouring column vectors. \f$ \vec{x_1}, \vec{x_2},
\vec{x_3} \f$. The cosine of alpha is calculated with:
\f[
  \cos(\alpha) =
  \frac{(\vec{x_2} - \vec{x_1}) \cdot (\vec{x_3} - \vec{x_2})}
       {\Vert\vec{x_2} - \vec{x_1}\Vert \cdot \Vert\vec{x_3} - \vec{x_2}\Vert}
\f]
The 2-norm is taken as the norm.
\code
 x_1
   \
    \
     x_2 ----- x_3
      \  )
       \
\endcode
*/
class AngleCosineCalculator : public TraxelsFeatureCalculator {
 public:
  AngleCosineCalculator() {};
  virtual ~AngleCosineCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  DiffCalculator diff_calculator_;
  EuclideanNormCalculator norm_calculator_;
  DotProductCalculator dot_product_calculator_;
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
//// class ChildDeceleration
////
/**
\brief calculates the deceleration of the childrens movement in a cell division

Calculates the ratios of the childrens squared move distance between t and t+1
where t is the time of the cell division. It returns a row vector with two
values. One for the first child, one for the second child.
\code
       t2--t3
      /
t1--t0
      \
       t4--t5
\endcode
*/
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
//// class CovarianceCalculator
////
/**
\brief calculates the covariance or inverse covariance matrix for the column
  vectors of the feature matrix

\tparam INV determines wheter the covariance matrix or the inverse covariance
  matrix is returned
*/
template<bool INV>
class CovarianceCalculator : public TraxelsFeatureCalculator {
 public:
  CovarianceCalculator() {};
  virtual ~CovarianceCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
 protected:
  static const std::string name_;
};

////
//// class SquaredMahalanobisCalculator
////
/**
\brief calculates the mahalanobis distance between each column vector and either
  the 0-vector or the mean vector (see description).

The distance matrix can be set in the calculate() method. If it is given, the
distance between the 0-vector and the column vectors is calculated. If no
distance matrix is given the inverse covariance matrix of the column vectors is
taken. In this case the mahalanobis distance to the mean vector is calculated.

With a specified distance matrix \f$A\f$:
\f[
  (\vec{x_1}, \dots, \vec{x_n})
  \mapsto
  (\vec{x_1}^T A \vec{x_1}, \dots, \vec{x_n}^T A \vec{x_n})
\f]
With no distance matrix given:
\f[
  (\vec{x_1}, \dots, \vec{x_n})
  \mapsto
  ((\vec{x_1}-\vec{\mu})^T \Sigma^{-1} (\vec{x_1}-\vec{\mu}), \dots, (\vec{x_n}-\vec{\mu})^T \Sigma^{-1} (\vec{x_n}-\vec{\mu}))
\f]

Where \f$\vec{\mu}\f$ denotes the mean vector of all column vectors.
\f$\Sigma^{-1}\f$ is the inverse of the covariance matrix of the row vectors.
*/
class SquaredMahalanobisCalculator : public TraxelsFeatureCalculator {
 public:
  SquaredMahalanobisCalculator() {};
  virtual ~SquaredMahalanobisCalculator() {};
  virtual const std::string& name() const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
  ) const;
  virtual void calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix,
    const FeatureMatrix& inv_covariance_matrix
  ) const;
 protected:
  DeviationCalculator deviation_calculator_;
  CovarianceCalculator<true> inv_covariance_calculator_;
  static const std::string name_;
};

////
//// class MVNOutlierCalculator
////
/**
\brief calculates count of outliers normalized to the track length

The outliers are calculated with the multivariant gaussian distribution. The
covariance matrix \f$\Sigma\f$ and the mean value \f$\vec{\mu}\f$ are calculated
with the column vectors. An outlier is defined as a column vector that differs
for more than three sigma from the mean value. This sigma threshold can also be
varied.

The count of outliers returned by MVNOutlierCalculator::calculate is normalized
to the count of column vectors in the matrix.

Definition of outlier:
\f[\sigma < (\vec{x}-\vec{\mu})^T \Sigma^{-1} (\vec{x}-\vec{\mu})\f]
*/
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
 protected:
  SquaredMahalanobisCalculator mahalanobis_calculator_;
  static const std::string name_;
};

} // end namespace features
} // end namespace pgmlink

#endif // PGMLINK_HIGHER_ORDER_FEATURES_H
