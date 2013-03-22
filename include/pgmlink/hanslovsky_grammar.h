#ifndef HANSLOVSKY_GRAMMAR_H
#define HANSLOVSKY_GRAMMAR_H

// stl headers
#include <vector>


// external headers
#include <armadillo>


// pgmlink headers

namespace pgmlink {

  typedef std::vector<std::string> feature_list;

  template <typename NodeFilter, typename ArcFilter>
  struct SubHypothesesGraph {
    typedef typename lemon::SubDigraph<HypothesesGraph::base_graph,
                                       typename property_map<NodeFilter, HypothesesGraph::base_graph>::type,
                                       typename property_map<ArcFilter, HypothesesGraph::base_graph>::type > type;
  };

  struct MapTypeValue {};
  struct MapTypeBool {};

  typedef SubHypothesesGraph<node_resolution_candidate, arc_resolution_candidate>::type SubResolver;

  class KMeans;

  template <typename T, typename U>
  void feature_array_to_arma_mat(const std::vector<T>& in,
				 arma::Mat<U>& out);

  template <typename T>
  void get_centers(const arma::Mat<T>& data,
		   const arma::Col<size_t> labels,
		   arma::Mat<T>& centers,
		   int k);

  class FeatureExtractorBase;

  class FeatureExtractorMCOMsFromPCOMs;

  class FeatureExtractorMCOMsFromMCOMs;

  class FeatureHandlerBase;

  class FeatureHandlerFromTraxels;

  class DistanceBase;

  class DistanceFromCOMs;

  class ResolveAmbiguousArcsBase;

  class ResolveAmbiguousArcsGreedy;

  class ResolveAmbiguousArcsPgm;

  class ReasonerMaxOneArc;

  class MergerResolver;

}

#endif /* HANSLOVSKY_GRAMMAR_H */
