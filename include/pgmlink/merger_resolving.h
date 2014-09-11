#ifndef MERGER_RESOLVING_H
#define MERGER_RESOLVING_H

// stl headers
#include <vector>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <string>
#include <map>


// external headers
#include <lemon/maps.h>
#include <lemon/adaptors.h>
#include <armadillo>
#include <boost/shared_ptr.hpp>
#include <vigra/multi_iterator_coupled.hxx>
#include <vigra/tinyvector.hxx>


// pgmlink headers
#include "pgmlink/hypotheses.h"
#include "pgmlink/event.h"
#include "pgmlink/traxels.h"
#include "pgmlink/reasoner.h"
#include "pgmlink/merger_resolving_grammar.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/feature.h"
#include "pgmlink/pgmlink_export.h"

/**
 * @brief Implementation of ideas for merger resolution in the HypothesesGraph environment.
 * @file
 *
 *
 * This header contains all the tools neccessary for resolve mergers on a HypothesesGraph. It provides specifications
 * of base classes for the use with the conservsation tracking. It is possible to derive classes from the base classes to
 * use the resolver for your own specific problem
 */

namespace pgmlink {

// typedef std::map<unsigned, arma::mat> IdCoordinateMap;

typedef std::map<std::pair<int, unsigned>, arma::mat> TimestepIdCoordinateMap;

typedef boost::shared_ptr<TimestepIdCoordinateMap > TimestepIdCoordinateMapPtr;

////
//// ClusteringMlpackBase
////

class ClusteringMlpackBase 
{
 protected:
  PGMLINK_EXPORT void copy_centers_to_feature_array(const arma::mat& centers, feature_array& c);
 public:
  PGMLINK_EXPORT virtual ~ClusteringMlpackBase() {}
  PGMLINK_EXPORT virtual feature_array operator()() = 0;
  PGMLINK_EXPORT virtual double score() const {return 0.0;}
};


////
//// KMeans
////
/**
 * @class KMeans 
 * @brief compatibility class for kMeans as an interface between feature_array and the mlpack library used for
 * kMeans
 *
 *
 * The library mlpack provides several clustering algorithms, one of which is kMeans. Instead of doing repetetive work
 * as rewriting kMeans would be, we are using the kMeans implementation provided in mlpack. To do so we need to convert
 * the data given in the form of a feature_array into an appropriate armadillo matrix (arma::mat), that can be used by
 * mlpack
 */
class KMeans 
: public ClusteringMlpackBase 
{
 private:
  KMeans();
  int k_;
  const feature_array& data_;
  // void copy_centers_to_feature_array(const arma::mat& centers, feature_array& c);
 public:
  // tested
  /**
   * @brief Constructor
   * @param [in] k number of clusters
   * @param [in] data feature_array storing data
   */
  PGMLINK_EXPORT KMeans(int k, const feature_array& data) 
  : k_(k), data_(data)
  {}

  // tested
  /**
   * @brief compute cluster centers and labels for datapoints
   * @returns feature_array that contains the coordinates of k clusters
   */
  PGMLINK_EXPORT virtual feature_array operator()();
};


class GMM 
: public ClusteringMlpackBase 
{
 private:
  GMM();
  int k_;
  int n_;
  const feature_array& data_;
  double score_;
  int n_trials_;
 public:
  // constructor needs to specify number of dimensions
  // for 2D data, ilastik provides coordinates with 3rd dimension 0
  // which will cause singular covariance matrix
  // therefore add option for dimensionality
  PGMLINK_EXPORT GMM(int k, int n, const feature_array& data, int n_trials=1) 
  : k_(k), n_(n), data_(data), score_(0.0), n_trials_(n_trials)
  {}

  PGMLINK_EXPORT virtual feature_array operator()();
  PGMLINK_EXPORT double score() const;
    
};


class GMMWithInitialized 
: public ClusteringMlpackBase 
{
 private:
  GMMWithInitialized();
  int k_;
  int n_;
  const feature_array& data_;
  double score_;
  int n_trials_;
  const std::vector<arma::vec>& means_;
  const std::vector<arma::mat>& covs_;
  const arma::vec& weights_;
 public:
  // constructor needs to specify number of dimensions
  // for 2D data, ilastik provides coordinates with 3rd dimension 0
  // which will cause singular covariance matrix
  // therefore add option for dimensionality
  PGMLINK_EXPORT GMMWithInitialized(int k, int n, const feature_array& data, int n_trials,
                                    const std::vector<arma::vec>& means, 
                                    const std::vector<arma::mat>& covs, const arma::vec& weights)
  : k_(k), n_(n), data_(data), score_(0.0), n_trials_(n_trials), means_(means), covs_(covs), weights_(weights) 
  {}

  PGMLINK_EXPORT virtual feature_array operator()();
  PGMLINK_EXPORT double score() const;
    
};


class GMMInitializeArma 
: public ClusteringMlpackBase 
{

 public:
  
  PGMLINK_EXPORT GMMInitializeArma(int k, const arma::mat& data, int n_trials=1, int n_iterations=30, double threshold=0.000001);
  PGMLINK_EXPORT ~GMMInitializeArma();

  PGMLINK_EXPORT virtual feature_array operator()();
  PGMLINK_EXPORT double score() const;

 private:
  GMMInitializeArma();
  int k_;
  const arma::mat& data_;
  double score_;
  int n_trials_;
  int n_iterations_;
  double threshold_;

};

    

////
//// helper functions
////
template <typename T, typename U>
// tested
/**
 * @brief Helper function to convert feature_array to arma::Mat.
 * @param [in] in original data; specifying T=float will make in a feature_array
 * @param [in,out] out arma::Mat<U> that holds the converted data. For the use in
 * KMeans specify U=double
 */
void feature_array_to_arma_mat(const std::vector<T>& in, arma::Mat<U>& out);


template <typename T, typename U>
void feature_array_to_arma_mat_skip_last_dimension(const std::vector<T>& in, arma::Mat<U>& out, unsigned int last_dimension);

  
template <typename T>
// tested
/**
 * @brief Helper function to calculate center coordinates from data assignments.
 * @param [in] data data points (coordinates)
 * @param [in] labels assignments after running kMeans
 * @param [in,out] centers arma::Mat to hold the coordinates of the cluster centers
 * @param [in] k number of clusters used for kMeans
 *
 * The mlpack kMeans implementation does not return the coordinates of the cluster centers.
 * The centers can be computed using the original data and the assignments.
 */
void get_centers(const arma::Mat<T>& data, const arma::Col<size_t> labels, arma::Mat<T>& centers, int k);


////
//// FeatureExtractorBase
////
/*
 * @brief Base class for feature extraction used when resolving merger nodes
 * @class FeatureExtractorBase
 *
 * 
 */
class FeatureExtractorBase {
 public:
  virtual std::vector<Traxel> operator()(Traxel& trax, size_t nMergers, unsigned int max_id) = 0;
 protected:
};


////
//// FeatureExtractorMCOMsFromPCOMs
////
class FeatureExtractorMCOMsFromPCOMs 
: public FeatureExtractorBase 
{
 public:
  PGMLINK_EXPORT virtual std::vector<Traxel> operator()(Traxel& trax, size_t nMergers, unsigned int max_id);
};

  
////
//// FeatureExtractorMCOMsFromMCOMs
////
class FeatureExtractorMCOMsFromMCOMs 
: public FeatureExtractorBase
{
 public:
  PGMLINK_EXPORT virtual std::vector<Traxel> operator()(Traxel& trax, size_t nMergers, unsigned int max_id);
};
    

////
//// FeatureExtractorMCOMsFromKMeans
////
class FeatureExtractorMCOMsFromKMeans 
: public FeatureExtractorBase 
{
 public:
  PGMLINK_EXPORT virtual std::vector<Traxel> operator()(Traxel& trax, size_t nMergers, unsigned int max_id);
};


////
//// FeatureExtractorMCOMsFromGMM
////
class FeatureExtractorMCOMsFromGMM 
: public FeatureExtractorBase
{
 private:
  int n_dim_;
 public:
  PGMLINK_EXPORT FeatureExtractorMCOMsFromGMM(int n_dim) 
  : n_dim_(n_dim) 
  {}
  
  PGMLINK_EXPORT virtual std::vector<Traxel> operator()(Traxel& trax, size_t nMergers, unsigned int max_id);
};

////
//// FeatureExtractorArmadillo
////
class FeatureExtractorArmadillo 
: public FeatureExtractorBase 
{
 public:
  PGMLINK_EXPORT FeatureExtractorArmadillo(TimestepIdCoordinateMapPtr coordinates);
  PGMLINK_EXPORT virtual std::vector<Traxel> operator()(Traxel& trax, size_t nMergers, unsigned int max_id);
 private:
  FeatureExtractorArmadillo();
  TimestepIdCoordinateMapPtr coordinates_;
};
  


////
//// DistanceBase
////
class DistanceBase {
 public:
  virtual double operator()(const HypothesesGraph& g, HypothesesGraph::Node from, HypothesesGraph::Node to) = 0;
  virtual double operator()(Traxel from,  Traxel to) = 0;
};


////
//// DistanceFromCOMs
////
class DistanceFromCOMs 
: public DistanceBase 
{
 public:
  PGMLINK_EXPORT virtual double operator()(const HypothesesGraph& g, HypothesesGraph::Node from, HypothesesGraph::Node to);
  PGMLINK_EXPORT virtual double operator()(Traxel from,  Traxel to);
};

  
////
//// FeatureHandlerBase
////
class FeatureHandlerBase 
{
 public:
  PGMLINK_EXPORT 
  void add_arcs_for_replacement_node(HypothesesGraph& g,
                                     HypothesesGraph::Node n,
                                     const std::vector<HypothesesGraph::base_graph::Arc>& sources,
                                     const std::vector<HypothesesGraph::base_graph::Arc>& targets,
                                     DistanceBase& distance);
 public:
  PGMLINK_EXPORT 
  virtual void operator()(HypothesesGraph& g,
                          HypothesesGraph::Node n,
                          std::size_t n_merger,
                          unsigned int max_id,
                          int timestep,
                          const std::vector<HypothesesGraph::base_graph::Arc>& sources,
                          const std::vector<HypothesesGraph::base_graph::Arc>& targets,
                          std::vector<unsigned int>& new_ids
                          ) = 0;
};

  
////
//// FeatureHandlerFromTraxels
////
class FeatureHandlerFromTraxels 
: public FeatureHandlerBase 
{
 private:
  FeatureExtractorBase& extractor_;
  DistanceBase& base_;
  // FeatureHandlerFromTraxelsMCOMsFromPCOMs() {};
 public:
  PGMLINK_EXPORT 
  FeatureHandlerFromTraxels(FeatureExtractorBase& extractor,
                            DistanceBase& base)
  : extractor_(extractor), base_(base)
  {}

  PGMLINK_EXPORT 
  virtual void operator()(HypothesesGraph& g,
                          HypothesesGraph::Node n,
                          std::size_t n_merger,
                          unsigned int max_id,
                          int timestep,
                          const std::vector<HypothesesGraph::base_graph::Arc>& sources,
                          const std::vector<HypothesesGraph::base_graph::Arc>& targets,
                          std::vector<unsigned int>& new_ids
                          );
};


////
//// ResolveAmbiguousArcsBase
////
class ResolveAmbiguousArcsBase {
 public:
  virtual HypothesesGraph& operator()(HypothesesGraph* g) = 0;
};
  

////
//// ResolveAmbiguousArcsGreedy
////
class ResolveAmbiguousArcsGreedy 
: public ResolveAmbiguousArcsBase 
{
 public:
  PGMLINK_EXPORT virtual HypothesesGraph& operator()(HypothesesGraph* g);
};


////
//// ReasonerMaxOneArc
////
class ReasonerMaxOneArc : public Reasoner {
};
  

////
//// ResolveAmbiguousArcsPgm
////
class ResolveAmbiguousArcsPgm : public ReasonerMaxOneArc, private ResolveAmbiguousArcsBase {
};


////
//// MergerResolver
////
/**
 * @brief Resolve mergers on a HypothesesGraph
 *
 * Using HypothesesGraph and property_map it is possible to build an algorithm that is capable of merger detection.
 * However, to fully solve the merger problem, the mergers need to be resolved into new objects and the tracking
 * has to be fed with the additional information from those new objects. This class gives an implementation that
 * is as general as possible to allow for application in various settings.
 * The model must provide a HypothesesGraph and the properties node_active2, arc_active, arc_distance.
 * The classes for application in the conservation tracking environment are provided. For the use in other settings, 
 * the appropriate classes have to be specified accordingly.
 */
   
class MergerResolver 
{
 private:
  HypothesesGraph* g_;
    
  // default constructor should be private (no object without specified graph allowed)
  MergerResolver();
    
  // collect arcs from ArcIterator and store them in vector
  // tested
  template <typename ArcIterator>
  void collect_arcs(ArcIterator,
                    std::vector<HypothesesGraph::base_graph::Arc>&);

  // template <typename ClusteringAlg>
  // do a more general way!
  // void calculate_centers(HypothesesGraph::Node,
  // int nMergers);

  // Add arcs to nodes created to replace merger node.
  void add_arcs_for_replacement_node(HypothesesGraph::Node node,
                                     Traxel& trax,
                                     std::vector<HypothesesGraph::base_graph::Arc> src,
                                     std::vector<HypothesesGraph::base_graph::Arc> dest,
                                     DistanceBase& distance);

  // Deactivate arcs of merger node.
  void deactivate_arcs(std::vector<HypothesesGraph::base_graph::Arc> arcs);

  // Deactivate all resolved merger nodes
  void deactivate_nodes(std::vector<HypothesesGraph::Node> nodes);

  // Get maximum id for given timestep
  unsigned int get_max_id(int ts);

  // Split merger node into appropiately many new nodes.
  void refine_node(HypothesesGraph::Node,
                   std::size_t,
                   FeatureHandlerBase& handler);

 public:
  PGMLINK_EXPORT MergerResolver(HypothesesGraph* g) 
  : g_(g)
  {
    if (!g_)
      throw std::runtime_error("HypotesesGraph* g_ is a null pointer!");
    if (!g_->has_property(merger_resolved_to()))
      g_->add(merger_resolved_to());
    if (!g_->has_property(node_active2()))
      throw std::runtime_error("HypothesesGraph* g_ does not have property node_active2!");
    if (!g_->has_property(arc_active()))
      throw std::runtime_error("HypothesesGraph* g_ does not have property arc_active!");
    if (!g_->has_property(arc_distance()))
      throw std::runtime_error("HypothesesGraph* g_ does not have property arc_distance!");
    if (!g_->has_property(node_originated_from()))
      g_->add(node_originated_from());
    if (!g_->has_property(node_resolution_candidate()))
      g_->add(node_resolution_candidate());
    if (!g_->has_property(arc_resolution_candidate()))
      g_->add(arc_resolution_candidate());
  }

  PGMLINK_EXPORT HypothesesGraph* resolve_mergers(FeatureHandlerBase& handler);
};


////
//// given a graph, do retracking
////
PGMLINK_EXPORT void resolve_graph(HypothesesGraph& src, HypothesesGraph& dest, boost::function<double(const double)> transition, double ep_gap, bool with_tracklets,
                   const double transition_parameter=5, const bool with_constraints=true);
// void resolve_graph(HypothesesGraph& src, HypothesesGraph& dest);

  
////
//// transfer graph to graph containing only subset of nodes based on tags
////
template <typename NodePropertyTag, typename ArcPropertyTag>
void copy_hypotheses_graph_subset(const HypothesesGraph& src,
                                  HypothesesGraph& dest,
                                  std::map<HypothesesGraph::Node, HypothesesGraph::Node>& nr,
                                  std::map<HypothesesGraph::Arc, HypothesesGraph::Arc>& ar,
                                  std::map<HypothesesGraph::Node, HypothesesGraph::Node>& ncr,
                                  std::map<HypothesesGraph::Arc, HypothesesGraph::Arc>& acr
                                  );


template <typename PropertyTag, typename KeyType>
void translate_property_value_map(const HypothesesGraph& src,
                                  const HypothesesGraph& dest,
                                  std::map<KeyType, KeyType> dict
                                  );


template <typename PropertyTag, typename KeyType>
void translate_property_bool_map(const HypothesesGraph& src,
                                 const HypothesesGraph& dest,
                                 std::map<KeyType,KeyType> dict
                                 );


template <typename NodePropertyTag, typename ArcPropertyTag>
void get_subset(const HypothesesGraph& src,
                HypothesesGraph& dest,
                HypothesesGraph::NodeMap<HypothesesGraph::Node>& nr,
                HypothesesGraph::ArcMap<HypothesesGraph::Arc>& ar,
                HypothesesGraph::NodeMap<HypothesesGraph::Node>& ncr,
                HypothesesGraph::ArcMap<HypothesesGraph::Arc>& acr
                );

  
PGMLINK_EXPORT double calculate_BIC(int k, int n, double weight);


PGMLINK_EXPORT void gmm_priors_and_centers(const feature_array& data, feature_array& priors, feature_array& centers, int k_max, int n, double weight);

PGMLINK_EXPORT void gmm_priors_and_centers_arma(const arma::mat& data, feature_array& priors, feature_array& centers, int k_max, int ndim, double regularization_weight);


////
//// duplicate division nodes in subset graph
////
PGMLINK_EXPORT void duplicate_division_nodes(HypothesesGraph& graph,
                              std::map<HypothesesGraph::Node, HypothesesGraph::Node>& division_splits,
                              std::map<HypothesesGraph::Arc, HypothesesGraph::Arc>& arc_cross_reference);


////
//// merge previously split divisions after inference
////
PGMLINK_EXPORT void merge_split_divisions(const HypothesesGraph& graph,
                           std::map<HypothesesGraph::Node, HypothesesGraph::Node>& division_splits,
                           std::map<HypothesesGraph::Arc, HypothesesGraph::Arc>& arc_cross_reference);


PGMLINK_EXPORT void calculate_gmm_beforehand(HypothesesGraph& g, int n_trials, int n_dimensions);


// extract coordinates in arma::mat
template<int N, typename T>
void extract_coordinates(TimestepIdCoordinateMapPtr coordinates,
                         const vigra::MultiArrayView<N, T>& image,
                         const vigra::TinyVector<long int, N>& offsets,
                         const Traxel& trax);

template<int N, typename T>
void extract_coord_by_timestep_id(TimestepIdCoordinateMapPtr coordinates,
                                  const vigra::MultiArrayView<N, T>& image,
                                  const vigra::TinyVector<long int, N>& offsets,
                                  const size_t timestep,
                                  const size_t traxel_id,
                                  const size_t traxel_size);


////
//// IMPLEMENTATIONS ////
////



template <typename T, typename U>
void feature_array_to_arma_mat(const std::vector<T>& in, arma::Mat<U>& out) {
  int stepSize = out.n_rows;
  int n = out.n_cols;
  if (stepSize*n != (int)in.size()) {
    throw std::range_error("Source vector dimension and matrix dimensions do not agree!");
  }
  int count = 0;
  typename std::vector<T>::const_iterator srcIt = in.begin();
  while (count < n) {
    arma::Col<U> col(stepSize);
    std::copy(srcIt, srcIt+stepSize, col.begin());
    out.col(count) = col;
    ++count;
    srcIt += stepSize;
  }
}

  
template <typename T, typename U>
void feature_array_to_arma_mat_skip_last_dimension(const std::vector<T>& in, arma::Mat<U>& out, unsigned int last_dimension) {
  unsigned int stepSize = out.n_rows;
  unsigned int n = out.n_cols;
  unsigned int count = 0;
  assert(last_dimension*n == in.size());
  assert(stepSize == last_dimension-1);
  typename std::vector<T>::const_iterator srcIt = in.begin();
  while (count < n) {
    arma::Col<U> col(stepSize);
    std::copy(srcIt, srcIt+stepSize, col.begin());
    out.col(count) = col;
    ++count;
    srcIt += last_dimension;
  }
}

  
template <typename T>
void get_centers(const arma::Mat<T>& data, const arma::Col<size_t> labels, arma::Mat<T>& centers, int k) {
  arma::Col<size_t>::const_iterator labelIt = labels.begin();
  std::vector<int> clusterSize(k, 0);
  centers.zeros();
  for (unsigned int n = 0; n < data.n_cols; ++n, ++labelIt) {
    ++clusterSize[*labelIt];
    centers.col(*labelIt) = centers.col(*labelIt) + data.col(n);
  }
  for (int i = 0; i < k; ++i) {
    centers.col(i) /= clusterSize[i];
  }
}
  
  
template <typename ArcIterator>
void MergerResolver::collect_arcs(ArcIterator arcIt,
                                  std::vector<HypothesesGraph::base_graph::Arc>& res) {
  assert(res.size() == 0);
  // check if arc is active
  property_map<arc_active, HypothesesGraph::base_graph>::type& arc_active_map = g_->get(arc_active());
  for (; arcIt != lemon::INVALID; ++arcIt) {
    if (arc_active_map[arcIt]) {
      res.push_back(arcIt);
    }
  }
}

  
template <typename NodePropertyTag, typename ArcPropertyTag>
void copy_hypotheses_graph_subset(const HypothesesGraph& src,
                                  HypothesesGraph& dest,
                                  std::map<HypothesesGraph::Node, HypothesesGraph::Node>& nr,
                                  std::map<HypothesesGraph::Arc, HypothesesGraph::Arc>& ar,
                                  std::map<HypothesesGraph::Node, HypothesesGraph::Node>& ncr,
                                  std::map<HypothesesGraph::Arc, HypothesesGraph::Arc>& acr
                                  ) {
  LOG(logDEBUG) << "copy_hypotheses_graph_subset(): entered";
  property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = src.get(node_traxel());

  dest.add(node_originated_from());

  typedef typename property_map<NodePropertyTag, HypothesesGraph::base_graph>::type NodeFilter;
  typedef typename property_map<ArcPropertyTag, HypothesesGraph::base_graph>::type ArcFilter;
  typedef property_map<node_originated_from, HypothesesGraph::base_graph>::type OriginMap;
    
  property_map<node_timestep, HypothesesGraph::base_graph>::type& time_map = src.get(node_timestep());
  NodeFilter& node_filter_map = src.get(NodePropertyTag());
  ArcFilter& arc_filter_map = src.get(ArcPropertyTag());
  OriginMap& origin_map_src = src.get(node_originated_from());
  OriginMap& origin_map_dest = dest.get(node_originated_from());


  for (typename NodeFilter::TrueIt nodeIt(node_filter_map); nodeIt != lemon::INVALID; ++nodeIt) {
    HypothesesGraph::Node node = dest.add_node(time_map[nodeIt]);
    nr[nodeIt] = node;
    ncr[node] = nodeIt;
    origin_map_dest.set(node, origin_map_src[nodeIt]);
  }

  for (typename ArcFilter::TrueIt arcIt(arc_filter_map); arcIt != lemon::INVALID; ++arcIt) {
    HypothesesGraph::Node from = src.source(arcIt);
    HypothesesGraph::Node to = src.target(arcIt);
    if (nr.count(from) == 0) {
      HypothesesGraph::Node node = dest.add_node(time_map[from]);
      nr[from] = node;
      ncr[node] = from;
      origin_map_dest.set(node, origin_map_src[from]);
      LOG(logDEBUG3) << "copy_hypotheses_graph_subset(): copied node: " << traxel_map[from];
    }
    if (nr.count(to) == 0) {
      HypothesesGraph::Node node = dest.add_node(time_map[to]);
      nr[to] = node;
      ncr[node] = to;
      origin_map_dest.set(node, origin_map_src[to]);
      LOG(logDEBUG3) << "copy_hypotheses_graph_subset(): copied node: " << traxel_map[to];
    }
    HypothesesGraph::Arc arc = dest.addArc(nr[from], nr[to]);
    ar[arcIt] = arc;
    acr[arc] = arcIt;
    LOG(logDEBUG3) << "copy_hypotheses_graph_subset(): copied arc " << src.id(arcIt) << ": ("
                   << traxel_map[src.source(arcIt)] << ','
                   << traxel_map[src.target(arcIt)] << ')';
  }
  LOG(logDEBUG) << "copy_hypotheses_graph_subset(): done";
}

  
/* template <typename PropertyTag, typename KeyType, typename MapType>
   void translate_property_map(const HypothesesGraph& src, const HypothesesGraph& dest, const std::map<KeyType, KeyType> dict); */

  
template <typename PropertyTag, typename KeyType>
void translate_property_value_map(const HypothesesGraph& src,
                                  const HypothesesGraph& dest,
                                  std::map<KeyType, KeyType> dict
                                  ) {
  typedef typename property_map<PropertyTag, HypothesesGraph::base_graph>::type IterableMap;
  IterableMap& src_map = src.get(PropertyTag());
  IterableMap& dest_map = dest.get(PropertyTag());
  typename IterableMap::ValueIt v_it;
  for (v_it = src_map.beginValue(); v_it != src_map.endValue(); ++ v_it) {
    typename IterableMap::ItemIt i_it(src_map, *v_it);
    for (; i_it != lemon::INVALID; ++i_it) {
      if (dict.count(i_it) > 0) {
        dest_map.set(dict[i_it], *v_it);
      }
    }
  }
}


template <typename PropertyTag, typename KeyType>
void translate_property_bool_map(const HypothesesGraph& src,
                                 const HypothesesGraph& dest,
                                 std::map<KeyType,KeyType> dict
                                 ) {
  // use c++11 and change for loop to:
  // for(bool b : {false, true});

  // C++11 !!!

  LOG(logDEBUG) << "translate_property_bool_map(): entering";
    
  bool const bools[] = {false, true};
  typedef typename property_map<PropertyTag, HypothesesGraph::base_graph>::type IterableMap;
  IterableMap& src_map = src.get(PropertyTag());
  IterableMap& dest_map = dest.get(PropertyTag());
  for (bool const* v_it(bools); v_it != bools + 2; ++v_it) {
    typename IterableMap::ItemIt i_it(src_map, *v_it);
    for(; i_it != lemon::INVALID; ++i_it) {
      if (dict.count(i_it) > 0) {
        dest_map.set(dict[i_it], *v_it);
      }
    }
  }
}

template <typename NodePropertyTag, typename ArcPropertyTag>
void get_subset(const HypothesesGraph& src,
                HypothesesGraph& dest,
                HypothesesGraph::NodeMap<HypothesesGraph::Node>& nr,
                HypothesesGraph::ArcMap<HypothesesGraph::Arc>& ar,
                HypothesesGraph::NodeMap<HypothesesGraph::Node>& ncr,
                HypothesesGraph::ArcMap<HypothesesGraph::Arc>& acr
                ) {
  typedef typename property_map<NodePropertyTag, HypothesesGraph::base_graph>::type NodeFilter;
  typedef typename property_map<ArcPropertyTag, HypothesesGraph::base_graph>::type ArcFilter;

  typedef lemon::SubDigraph<HypothesesGraph::base_graph, NodeFilter, ArcFilter> CopyGraph;
  //typedef HypothesesGraph::base_graph GRAPH;
    
  NodeFilter& node_filter_map = src.get(NodePropertyTag());
  ArcFilter& arc_filter_map = src.get(ArcPropertyTag());
  CopyGraph sub(src, node_filter_map, arc_filter_map);

  lemon::digraphCopy<CopyGraph, HypothesesGraph::base_graph>(sub,dest).nodeRef(nr).nodeCrossRef(ncr).arcRef(ar).arcCrossRef(acr).run();
}


template<int N, typename T>
void extract_coordinates(TimestepIdCoordinateMapPtr coordinates,
                         const vigra::MultiArrayView<N, T>& image,
                         const vigra::TinyVector<long int, N>& offsets,
                         const Traxel& trax) {
  LOG(logDEBUG3) << "extract_coordinates -- entered for " << trax;
  typedef typename vigra::CoupledIteratorType<N, T>::type Iterator;
  Iterator start = createCoupledIterator(image);
  Iterator end = start.getEndIterator();
  arma::mat& coord = (*coordinates)[std::make_pair(trax.Timestep, trax.Id)];
  coord = arma::mat(N, trax.features.find("count")->second[0]);
  // coord stores coordinates: each row is a spatial dimension and each column is a pixel
  LOG(logDEBUG4) << "extract_coordinates -- coordinate matrix has "
                 << coord.n_rows << " rows and "
                 << coord.n_cols << " cols. The traxel size is "
                 << trax.features.find("count")->second[0] << ".";
  {
    size_t index = 0;
    for (; start != end; ++start) {
      if (start.template get<1>() == trax.Id) {
        const vigra::TinyVector<long int, N>& position = start.template get<0>();
        for (int i = 0; i < N; ++i) {
          coord(i, index) = position[i] + offsets[i];
        }        
        ++index;
      } else {
        continue;
      }
    }
    assert(index == coord.n_cols);
  }
  LOG(logDEBUG3) << "extract_coordinates -- done";
}

template<int N, typename T>
void extract_coord_by_timestep_id(TimestepIdCoordinateMapPtr coordinates,
                                  const vigra::MultiArrayView<N, T>& image,
                                  const vigra::TinyVector<long int, N>& offsets,
                                  const size_t timestep,
                                  const size_t traxel_id,
                                  const size_t traxel_size) {
  LOG(logDEBUG3) << "extract_coordinates -- entered for " << traxel_id;
  typedef typename vigra::CoupledIteratorType<N, T>::type Iterator;
  Iterator start = createCoupledIterator(image);
  Iterator end = start.getEndIterator();
  arma::mat& coord = (*coordinates)[std::make_pair(timestep, traxel_id)];
  coord = arma::mat(N, traxel_size);
  // coord stores coordinates: each row is a spatial dimension and each column is a pixel
  LOG(logDEBUG4) << "extract_coordinates -- coordinate matrix has "
                 << coord.n_rows << " rows and "
                 << coord.n_cols << " cols. The traxel size is "
                 << traxel_size << ".";
  {
    size_t index = 0;
    for (; start != end; ++start) {
      if (start.template get<1>() == traxel_id) {
        const vigra::TinyVector<long int, N>& position = start.template get<0>();
        for (int i = 0; i < N; ++i) {
          coord(i, index) = position[i] + offsets[i];
        }        
        ++index;
      } else {
        continue;
      }
    }
    LOG(logDEBUG4) << "matrix has " << index << " columns while there should be "
                  << coord.n_cols << " columns";
    assert(index == coord.n_cols);
  }
  LOG(logDEBUG3) << "extract_coordinates -- done";
}

/* template <typename ClusteringAlg>
   void MergerResolver::calculate_centers(HypothesesGraph::Node node,					 
   int nMergers) {
   // get traxel map from graph to access traxel
   property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g_->get(node_traxel());
   Traxel trax = traxel_map[node];
   feature_array mergerCOMs;

   // assert mergerCOMs does not exist
   assert(trax.features.find("mergerCOMs") == trax.features.end());
    
   // check for feature possibleCOMs. If present, read appropriate coordinates. Otherwise calculate mergerCOMs from coordinate list
   if (trax.features.find("possibleCOMs") != trax.features.end()) {
   int index1 = 3*nMergers*(nMergers-1)/2;
   int index2 = 3*nMergers*(nMergers+1)/2;
   mergerCOMs.assign(trax.features["possibleCOMs"].begin() + index1, trax.features["possibleCOMs"].begin() + index2);
   } else {
   // throw exception if list of coordinates is not stored int traxel
   if (trax.features.find("Coord<ValueList>") == trax.features.end()) {
   throw std::runtime_error("List of coordinates not stored in traxel!");
   }
   // calculate merger centers using clustering algorithm calg
   ClusteringAlg calg(nMergers, trax.features["Coord<ValueList>"]);
   mergerCOMs = calg();
   }
   trax.features["mergerCOMs"] = mergerCOMs;
   assert((int)mergerCOMs.size() == 3*nMergers);
   traxel_map.set(node, trax);
   } */

  
  
}


#endif /* MERGER_RESOLVING_H */
