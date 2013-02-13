// stl headers
#include <vector>
#include <stdexcept>
#include <cassert>
#include <algorithm>


// external headers
#include <lemon/maps.h>
#include <armadillo>
#include <mlpack/core.hpp>
#include <mlpack/methods/kmeans/kmeans.hpp>


// pgmlink headers
#include "pgmlink/hanslovsky.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/event.h"
#include "pgmlink/traxels.h"


namespace pgmlink {
  ////
  //// KMeans
  ////
  feature_array KMeans::operator()() {
    mlpack::kmeans::KMeans<> kMeans;
    int n = data_.size()/3;
    arma::mat data(3, n);
    arma::mat centers(3, k_);
    arma::Col<size_t> labels;
    feature_array_to_arma_mat(data_, data);
    kMeans.Cluster(data, k_, labels);
    get_centers(data, labels, centers, k_);
    feature_array fa_centers(3*k_);
    copy_centers_to_feature_array(centers, fa_centers);
    return fa_centers;
  }
  

  void KMeans::copy_centers_to_feature_array(const arma::mat& centers, feature_array& c) {
    int n = centers.n_cols;
    int stepSize = centers.n_rows;
    
    if (stepSize*n != (int)c.size()) {
      throw std::range_error("Source matrix dimensions and vector dimension do not agree!");
    }

    feature_array::iterator it = c.begin();
    for (int i = 0; i < n; ++i, it += stepSize) {
      arma::vec c = centers.col(i);
      std::copy(c.begin(), c.end(), it);
    }
  }
 


  ////
  //// FeatureExtractorMCOMsFromPCOMs
  ////
  std::vector<Traxel> FeatureExtractorMCOMsFromPCOMs::operator()(
								 Traxel trax,
								 size_t nMerger,
								 unsigned int start_id
								 ) {
    std::map<std::string, feature_array>::iterator it = trax.features.find("possibleCOMs");
    assert(it != trax.features.end());
    std::vector<Traxel> res;
    unsigned int index1 = 3*(nMerger*(nMerger-1))/2;
    unsigned int index2 = 3*(nMerger*(nMerger+1))/2;
    feature_array range(it->second.begin()+index1, it->second.begin()+index2);
    for (unsigned int n = 0; n < nMerger; ++n, ++start_id) {
      trax.Id = start_id;
      trax.features["com"] = feature_array(range.begin()+(3*n), range.begin()+(3*(n+1)));
      res.push_back(trax);
    }
    return res;
  }


  ////
  //// FeatureExtractorMCOMsFromMCOMs
  ////
  std::vector<Traxel> FeatureExtractorMCOMsFromMCOMs::operator()(
										  Traxel trax,
										  size_t nMerger,
										  unsigned int start_id
										  ) {
    std::map<std::string, feature_array>::iterator it = trax.features.find("mergerCOMs");
    assert(it != trax.features.end());
    std::vector<Traxel> res;
    for (unsigned int n = 0; n < nMerger; ++n, ++start_id) {
      trax.Id = start_id;
      trax.features["com"] = feature_array(it->second.begin()+(3*n), it->second.begin()+(3*(n+1)));
      res.push_back(trax);
    }
    return res;
  }

  
  ////
  //// DistanceFromCOMs
  ////
  double DistanceFromCOMs::operator()(Traxel from, Traxel to) {
    return from.distance_to(to);
  }

  ////
  //// ReasonerMaxOneArc
  ////
  void ReasonerMaxOneArc::formulate(const HypothesesGraph& g) {
  }
  void ReasonerMaxOneArc::infer() {
  }
  void ReasonerMaxOneArc::conclude(HypothesesGraph& g) {
  }


  ////
  //// MergerResolver
  ////
  void MergerResolver::add_arcs_for_replacement_node(HypothesesGraph::Node node,
						     Traxel trax,
						     std::vector<HypothesesGraph::base_graph::Arc> src,
						     std::vector<HypothesesGraph::base_graph::Arc> dest,
						     DistanceBase& distance) {
    // add Arcs for new node that is a replacement node for merger node
    // get the property_maps needed for adding Arcs
    std::vector<HypothesesGraph::base_graph::Arc>::iterator it;
    property_map<arc_distance, HypothesesGraph::base_graph>::type& arc_distances = g_->get(arc_distance());
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g_->get(node_traxel());
    property_map<arc_active, HypothesesGraph::base_graph>::type& arc_active_map = g_->get(arc_active());

    // add incoming arcs
    for(it = src.begin(); it != src.end(); ++it) {
      HypothesesGraph::Node from = g_->source(*it);
      Traxel from_tr = traxel_map[from];
      double dist = distance(from_tr, trax);
      HypothesesGraph::Arc arc = g_->addArc(from, node);
      // add distance and activate arcs
      arc_distances.set(arc, dist);
      arc_active_map.set(arc, true);
    }

    // add outgoing arcs
    for(it = dest.begin(); it != dest.end(); ++it) {
      HypothesesGraph::Node to = g_->target(*it);
      Traxel to_tr = traxel_map[to];
      double dist = distance(trax, to_tr);
      HypothesesGraph::Arc arc = g_->addArc(node, to);
      // add distance and activate arcs
      arc_distances.set(arc, dist);
      arc_active_map.set(arc, true);
    }
  }

  void MergerResolver::deactivate_arcs(std::vector<HypothesesGraph::base_graph::Arc> arcs) {
    // deactivate Arcs provided by arcs
    // useful to deactivate Arcs of merger Node
    std::vector<HypothesesGraph::base_graph::Arc>::iterator it = arcs.begin();
    property_map<arc_active, HypothesesGraph::base_graph>::type& arc_active_map = g_->get(arc_active());
    for (; it != arcs.end(); ++it) {
      arc_active_map.set(*it, false);
    }
  }

  void MergerResolver::deactivate_nodes(std::vector<HypothesesGraph::Node> nodes) {
    // deactivate Nodes provided by nodes
    // needed to set all resolved merger nodes inactive
    std::vector<HypothesesGraph::Node>::iterator it = nodes.begin();
    property_map<node_active2, HypothesesGraph::base_graph>::type& node_active_map = g_->get(node_active2());
    for (; it != nodes.end(); ++it) {
      LOG(logDEBUG1) << "MergerResolver::deactivate_nodes(): setting Node " << g_->id(*it) << " property node_active2 to 0";
      node_active_map.set(*it, 0);
    }
  }

  unsigned int MergerResolver::get_max_id(int ts) {
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g_->get(node_traxel());
    property_map<node_timestep, HypothesesGraph::base_graph>::type& time_map = g_->get(node_timestep());
    property_map<node_timestep, HypothesesGraph::base_graph>::type::ItemIt timeIt(time_map, ts);
    unsigned int max_id = 0;
    for (; timeIt != lemon::INVALID; ++timeIt) {
      if (traxel_map[timeIt].Id > max_id) {
	max_id = traxel_map[timeIt].Id;
      }
    }
    return max_id;
  }

  void MergerResolver::refine_node(HypothesesGraph::Node node,
				   std::size_t nMerger,
				   FeatureExtractorBase& extractor,
				   DistanceBase& distance) {
    property_map<node_active2, HypothesesGraph::base_graph>::type& active_map = g_->get(node_active2());
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g_->get(node_traxel());
    property_map<node_timestep, HypothesesGraph::base_graph>::type& time_map = g_->get(node_timestep());
    property_map<merger_resolved_to, HypothesesGraph::base_graph>::type& resolved_map = g_->get(merger_resolved_to());
    property_map<node_originated_from, HypothesesGraph::base_graph>::type& origin_map = g_->get(node_originated_from());
    
    
    Traxel trax = traxel_map[node];
    int timestep = time_map[node];
    unsigned int max_id = get_max_id(timestep)+1;

    // get incoming and outgoing arcs for reorganizing arcs
    std::vector<HypothesesGraph::base_graph::Arc> sources;
    std::vector<HypothesesGraph::base_graph::Arc> targets;
    collect_arcs(HypothesesGraph::base_graph::InArcIt(*g_, node), sources);
    collect_arcs(HypothesesGraph::base_graph::OutArcIt(*g_, node), targets);


    // create new node for each of the objects merged into node
    std::vector<unsigned int> new_ids;
    std::vector<Traxel> ft = extractor(trax, nMerger, max_id);
    for (std::vector<Traxel>::iterator it = ft.begin(); it != ft.end(); ++it) {
      // set traxel features, most of which can be copied from the merger node
      // set new center of mass as calculated from GMM
      // add node to graph and activate it
      HypothesesGraph::Node newNode = g_->add_node(timestep);
      traxel_map.set(newNode, *it);
      active_map.set(newNode, 1);
      time_map.set(newNode, timestep);
      // add arc candidates for new nodes (todo: need to somehow choose which ones are active)
      add_arcs_for_replacement_node(newNode, *it, sources, targets, distance);
      // save new id from merger node to new_ids;
      new_ids.push_back(it->Id);
      // store parent (merger) node. this is used for creating the resolved_to event later
      origin_map.set(newNode, std::vector<unsigned int>(1, traxel_map[node].Id));
      
    }
    // deactivate incoming and outgoing arcs of merger node
    // merger node will be deactivated after pruning
    deactivate_arcs(sources);
    deactivate_arcs(targets);
    // save information on new ids in property map
    resolved_map.set(node, new_ids);
    LOG(logINFO) << "refine_node(): Node resolved to " << resolved_map[node].size() << " merger(s)";
    // g_->get(merger_resolved_to())[node] = new_ids;
  }


  HypothesesGraph* MergerResolver::resolve_mergers(FeatureExtractorBase& extractor,
						   DistanceBase& distance) {
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
	  // calculate_centers<ClusteringAlg>(active_itemIt, *active_valueIt);
	  // for each object create new node and set arcs to old merger node inactive (neccessary for pruning)
	  refine_node(active_itemIt, *active_valueIt, extractor, distance);
	  nodes_to_deactivate.push_back(active_itemIt);
	}
      }
    }
    // maybe keep merger nodes active for event extraction
    deactivate_nodes(nodes_to_deactivate);
    prune_inactive(*g_);
    return g_;
  }

}
