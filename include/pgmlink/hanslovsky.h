#ifndef HANSLOVSKY_H
#define HANSLOVSKY_H

// stl headers
#include <vector>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <string>
#include <map>


// external headers
#include <lemon/maps.h>
#include <armadillo>


// pgmlink headers
#include "pgmlink/hypotheses.h"
#include "pgmlink/event.h"
#include "pgmlink/traxels.h"


namespace pgmlink {

  typedef std::vector<std::string> feature_list;
  
  ////
  //// KMeans
  ////
  class KMeans {
  private:
    int k_;
    const feature_array& data_;
    
    void copy_centers_to_feature_array(const arma::mat& centers, feature_array& c);
  public:
    // tested
    KMeans(int k, const feature_array& data) :
      k_(k), data_(data) {}

    // tested
    feature_array operator()();
  };

  ////
  //// helper functions
  ////
  template <typename T, typename U>
  // tested
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

  
  template <typename T>
  // tested
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


  ////
  //// FeatureExtractorBase
  ////
  class FeatureExtractorBase {
  public:
    virtual std::vector<Traxel> operator()(Traxel trax, size_t nMergers, unsigned int max_id) = 0;
  };


  ////
  //// FeatureExtractorMCOMsFromPCOMs
  ////
  class FeatureExtractorMCOMsFromPCOMs : public FeatureExtractorBase {
  public:
    virtual std::vector<Traxel> operator()(Traxel trax, size_t nMergers, unsigned int max_id);
  };

  
  ////
  //// FeatureExtractorMCOMsFromMCOMs
  ////
  class FeatureExtractorMCOMsFromMCOMs : public FeatureExtractorBase {
  public:
    virtual std::vector<Traxel> operator()(Traxel trax, size_t nMergers, unsigned int max_id);
  };
    

  ////
  //// FeatureExtractorMCOMsFromKMeans
  ////
  class FeatureExtractorMCOMsFromKMeans : public FeatureExtractorBase {
  public:
    virtual std::vector<Traxel> operator()(Traxel trax, size_t nMergers, unsigned int max_id);
  };


  ////
  //// DistanceBase
  ////
  class DistanceBase {
  public:
    virtual double operator()(Traxel from, Traxel to) = 0;
  };


  ////
  //// DistanceFromCOMs
  ////
  class DistanceFromCOMs : public DistanceBase{
  public:
    virtual double operator()(Traxel from, Traxel to);
  };


  ////
  //// ReasonerMaxOneArc
  ////
  class ReasonerMaxOneArc : public Reasoner {
  private:
    
  public:
    ReasonerMaxOneArc();
    virtual void formulate(const HypothesesGraph& g);
    virtual void infer();
    virtual void conclude(HypothesesGraph& g);
  };


  ////
  //// MergerResolver
  ////
  class MergerResolver {
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
    // tested
    void add_arcs_for_replacement_node(HypothesesGraph::Node node,
				       Traxel trax,
				       std::vector<HypothesesGraph::base_graph::Arc> src,
				       std::vector<HypothesesGraph::base_graph::Arc> dest,
				       DistanceBase& distance);

    // Deactivate arcs of merger node.
    // tested
    void deactivate_arcs(std::vector<HypothesesGraph::base_graph::Arc> arcs);

    // Deactivate all resolved merger nodes
    // tested
    void deactivate_nodes(std::vector<HypothesesGraph::Node> nodes);

    // Get maximum id for given timestep
    // tested
    unsigned int get_max_id(int ts);

    // Split merger node into appropiately many new nodes.
    // tested
    void refine_node(HypothesesGraph::Node,
		     std::size_t,
		     FeatureExtractorBase& extractor,
		     DistanceBase& distance);
    
  public:
    MergerResolver(HypothesesGraph* g) : g_(g)
    {
      if (!g_)
	throw std::runtime_error("HypotesesGraph g_ is a null pointer!");
      if (!g_->has_property(merger_resolved_to()))
	g_->add(merger_resolved_to());
      if (!g_->has_property(node_active2()))
	throw std::runtime_error("HypothesesGraph g_ does not have property node_active2!");
      if (!g_->has_property(arc_active()))
	throw std::runtime_error("HypothesesGraph g_ does not have property arc_active!");
      if (!g_->has_property(arc_distance()))
	throw std::runtime_error("HypothesesGraph g_ does not have property arc_distance!");
      if (!g_->has_property(node_originated_from()))
	g_->add(node_originated_from());
    }
    HypothesesGraph* resolve_mergers(FeatureExtractorBase& extractor,
				     DistanceBase& distance);
  };
  
  template <typename ArcIterator>
  void MergerResolver::collect_arcs(ArcIterator arcIt,
				    std::vector<HypothesesGraph::base_graph::Arc>& res) {
    assert(res.size() == 0);
    for (; arcIt != lemon::INVALID; ++arcIt) {
      res.push_back(arcIt);
    }
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


#endif /* HANSLOVSKY_H */
