#ifndef HANSLOVSKY_H
#define HANSLOVSKY_H

// stl headers
#include <vector>
#include <stdexcept>
#include <cassert>
#include <algorithm>


// external headers
#include <lemon/maps.h>
#include <armadillo>


// pgmlink headers
#include "pgmlink/hypotheses.h"
#include "pgmlink/event.h"
#include "pgmlink/traxels.h"


namespace pgmlink {
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

    // Calculate cluster centers for merged cells if feature is not yet present at traxel.
    // tested
    template <typename ClusteringAlg>
    void calculate_centers(HypothesesGraph::Node,
			   ClusteringAlg calg,
			   int nMergers);

    // Add arcs to nodes created to replace merger node.
    // tested
    void add_arcs_for_replacement_node(HypothesesGraph::Node node,
				       Traxel trax,
				       std::vector<HypothesesGraph::base_graph::Arc> src,
				       std::vector<HypothesesGraph::base_graph::Arc> dest);

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
    void refine_node(HypothesesGraph::Node,
		     std::size_t);
    
  public:
    MergerResolver(HypothesesGraph* g) : g_(g) {}
    template <typename ClusteringAlg>
    HypothesesGraph* resolve_mergers(ClusteringAlg);
  };
  
  template <typename ArcIterator>
  void MergerResolver::collect_arcs(ArcIterator arcIt,
				    std::vector<HypothesesGraph::base_graph::Arc>& res) {
    assert(res.size() == 0);
    for (; arcIt != lemon::INVALID; ++arcIt) {
      res.push_back(arcIt);
    }
  }
  
  template <typename ClusteringAlg>
  void MergerResolver::calculate_centers(HypothesesGraph::Node node,
					 ClusteringAlg calg,
					 int nMergers) {
    // get traxel map from graph to access traxel
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g_->get(node_traxel());
    Traxel trax = traxel_map[node];
    feature_array mergerCOMs(nMergers*3);

    // assert mergerCOMs does not exist
    assert(trax.features.find("mergerCOMs") == trax.features.end());
    
    // check for feature possibleCOMs. If present, read appropriate coordinates. Otherwise calculate mergerCOMs from coordinate list
    if (trax.features.find("possibleCOMs") != trax.features.end()) {
      int index1 = nMergers*(nMergers-1)/2;
      int index2 = nMergers*(nMergers+1)/2;
      mergerCOMs.assign(trax.features["possibleCOMs"].begin() + index1, trax.features["possibleCOMs"].begin() + index2);
    } else {
      // throw exception if list of coordinates is not stored int traxel
      if (trax.features.find("Coord<ValueList>") == trax.features.end()) {
	throw std::runtime_error("List of coordinates not stored in traxel!");
      }
      // calculate merger centers using clustering algorithm calg
      mergerCOMs = calg(trax.features["Coord<ValueList>"]);
    }
    trax.features["mergerCOMs"] = mergerCOMs;
    traxel_map.set(node, trax);
  }
  
  template <typename ClusteringAlg>
  HypothesesGraph* MergerResolver::resolve_mergers(ClusteringAlg calg) {
    // extract property maps and iterators from graph
    property_map<node_active2, HypothesesGraph::base_graph>::type& active_map = g_->get(node_active2());
    property_map<node_active2, HypothesesGraph::base_graph>::type::ValueIt active_valueIt = active_map.beginValue();

    
    // iterate over mergers and replace merger nodes
    // keep track of merger nodes to deactivate them later
    std::vector<HypothesesGraph::Node> nodes_to_deactivate;
    for (; active_valueIt != active_map.endValue(); ++active_valueIt) {
      if (*active_valueIt > 1) {
	property_map<node_active2, HypothesesGraph::base_graph>::type::ItemIt active_itemIt(active_map, *active_valueIt);
	
	for (; active_itemIt != lemon::INVALID; ++active_itemIt) {
	  calculate_centers(active_itemIt, calg, *active_valueIt);
	  // for each object create new node and set arcs to old merger node inactive (neccessary for pruning)
	  refine_node(active_itemIt, *active_valueIt);
	  nodes_to_deactivate.push_back(active_itemIt);
	}
      }
    }
    deactivate_nodes(nodes_to_deactivate);
    prune_inactive(*g_);
    return g_;
  }
  
}


#endif /* HANSLOVSKY_H */
