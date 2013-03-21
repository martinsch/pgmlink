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
  typedef lemon::SubDigraph<HypothesesGraph::base_graph, NodeFilter, ArcFilter> SubHypothesesGraph

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
