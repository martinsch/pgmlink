// stl headers
#include <vector>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <iterator>


// external headers
#include <lemon/maps.h>
#include <armadillo>
#include <mlpack/core.hpp>
#include <mlpack/methods/kmeans/kmeans.hpp>
#include <mlpack/methods/gmm/gmm.hpp>


// pgmlink headers
#include "pgmlink/merger_resolving.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/event.h"
#include "pgmlink/traxels.h"


namespace pgmlink {
  ////
  //// ClusteringMlpackBase
  ////
  void ClusteringMlpackBase::copy_centers_to_feature_array(const arma::mat& centers, feature_array& c) {
    int n = centers.n_cols;
    int stepSize = centers.n_rows;
    
    if (stepSize*n != (int)c.size()) {
      throw std::range_error("Source matrix dimensions and vector dimension do not agree!");
    }

    feature_array::iterator it = c.begin();
    for (int i = 0; i < n; ++i, it += stepSize) {
      arma::vec col = centers.col(i);
      std::copy(col.begin(), col.end(), it);
    }
  }

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


  ////
  //// GMM
  ////
  feature_array GMM::operator()() {
    mlpack::gmm::GMM<> gmm(k_, n_);
    int n_samples = data_.size()/n_;
    arma::mat data(n_,n_samples);
    arma::Col<size_t> labels;
    feature_array_to_arma_mat(data_, data);
    score_ = gmm.Estimate(data, n_trials_);
    std::vector<arma::vec> centers = gmm.Means();
    feature_array fa_centers;
    for (std::vector<arma::vec>::iterator it = centers.begin(); it != centers.end(); ++it) {
      std::copy(it->begin(), it->end(), std::back_insert_iterator<feature_array >(fa_centers));
    }
    return fa_centers;
  }


  double GMM::score() const {
    return score_;
  }
  

  /*void KMeans::copy_centers_to_feature_array(const arma::mat& centers, feature_array& c) {
    int n = centers.n_cols;
    int stepSize = centers.n_rows;
    
    if (stepSize*n != (int)c.size()) {
      throw std::range_error("Source matrix dimensions and vector dimension do not agree!");
    }

    feature_array::iterator it = c.begin();
    for (int i = 0; i < n; ++i, it += stepSize) {
      arma::vec col = centers.col(i);
      std::copy(col.begin(), col.end(), it);
    }
  }*/
 


  ////
  //// FeatureExtractorMCOMsFromKMeans
  ////
  std::vector<Traxel> FeatureExtractorMCOMsFromKMeans::operator()(Traxel trax,
                                                                  size_t nMergers,
                                                                  unsigned int max_id
                                                                  ) {
    std::map<std::string, feature_array>::iterator it = trax.features.find("coordinates");
    assert(it != trax.features.end());
    std::vector<Traxel> res;
    KMeans kmeans(nMergers, it->second);
    trax.features["mergerCOMs"] = kmeans();
    FeatureExtractorMCOMsFromMCOMs extractor;
    return extractor(trax, nMergers, max_id);
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
      LOG(logDEBUG3) << "FeatureExtractorMCOMsFromPCOMs::operator()(): Appended traxel with com (" << trax.features["com"][0] << "," << trax.features["com"][1] << "," << trax.features["com"][2] << ") and id " << trax.Id;
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
  //// FeatureHandlerBase
  ////
  void FeatureHandlerBase::add_arcs_for_replacement_node(HypothesesGraph& g,
                                                         HypothesesGraph::Node n,
                                                         const std::vector<HypothesesGraph::base_graph::Arc>& sources,
                                                         const std::vector<HypothesesGraph::base_graph::Arc>& targets,
                                                         DistanceBase& distance) {
    // add Arcs for new node that is a replacement node for merger node
    // get the property_maps needed for adding Arcs
    std::vector<HypothesesGraph::base_graph::Arc>::const_iterator it;
    property_map<arc_distance, HypothesesGraph::base_graph>::type& arc_distances = g.get(arc_distance());
    // property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    property_map<arc_active, HypothesesGraph::base_graph>::type& arc_active_map = g.get(arc_active());
    property_map<arc_resolution_candidate, HypothesesGraph::base_graph>::type& arc_resolution_map = g.get(arc_resolution_candidate());


    // add incoming arcs
    for (it = sources.begin(); it != sources.end(); ++it) {
      HypothesesGraph::Node from = g.source(*it);
      double dist = distance(g, from, n);
      HypothesesGraph::Arc arc = g.addArc(from, n);
      arc_distances.set(arc, dist);
      arc_active_map.set(arc, true);
      arc_resolution_map.set(arc, true);
      // LOG?
    }

    // add outgoing arcs
    for (it = targets.begin(); it != targets.end(); ++it) {
      HypothesesGraph::Node to = g.target(*it);
      double dist = distance(g, n, to);
      HypothesesGraph::Arc arc = g.addArc(n, to);
      arc_distances.set(arc, dist);
      arc_active_map.set(arc, true);
      arc_resolution_map.set(arc, true);
      // LOG?
    }
  }
    

  ////
  //// FeatureHandlerFromTraxels
  ////
  void FeatureHandlerFromTraxels::operator()(
                                             HypothesesGraph& g,
                                             HypothesesGraph::Node n,
                                             std::size_t n_merger,
                                             unsigned int max_id,
                                             int timestep,
                                             const std::vector<HypothesesGraph::base_graph::Arc>& sources,
                                             const std::vector<HypothesesGraph::base_graph::Arc>& targets,
                                             std::vector<unsigned int>& new_ids
                                             ) {
    // property maps
    property_map<node_active2, HypothesesGraph::base_graph>::type& active_map = g.get(node_active2());
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    property_map<node_timestep, HypothesesGraph::base_graph>::type& time_map = g.get(node_timestep());
    property_map<node_originated_from, HypothesesGraph::base_graph>::type& origin_map = g.get(node_originated_from());
    property_map<node_resolution_candidate, HypothesesGraph::base_graph>::type& node_resolution_map = g.get(node_resolution_candidate());

    // traxel and vector of replacement traxels
    Traxel trax = traxel_map[n];
    std::vector<Traxel> ft = extractor_(trax, n_merger, max_id);
    // MAYBE LOG
    for (std::vector<Traxel>::iterator it = ft.begin(); it != ft.end(); ++it) {
      // set traxel features, most of which can be copied from the merger node
      // set new center of mass as calculated from GMM
      // add node to graph and activate it
      HypothesesGraph::Node new_node = g.add_node(timestep);
      // MAYBE LOG
      traxel_map.set(new_node, *it);
      active_map.set(new_node, 1);
      time_map.set(new_node, timestep);

      // add arc candidates for new nodes (todo: need to somehow choose which ones are active)
      this->add_arcs_for_replacement_node(g, new_node, sources, targets, base_);
      // save new id from merger node to new_ids;
      new_ids.push_back(it->Id);
      // store parent (merger) node. this is used for creating the resolved_to event later
      origin_map.set(new_node, std::vector<unsigned int>(1, trax.Id));
      LOG(logINFO) << "FeatureHandlerFromTraxels::operator(): added " << trax.Id << " to origin_map[" << g.id(new_node) << "]";
      node_resolution_map.set(new_node, true);
 
    }
  }

  
  ////
  //// DistanceFromCOMs
  ////
  double DistanceFromCOMs::operator()(const HypothesesGraph& g, HypothesesGraph::Node from, HypothesesGraph::Node to) {
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    return traxel_map[from].distance_to(traxel_map[to]);
  }

  double DistanceFromCOMs::operator()(Traxel from, Traxel to) {
    return from.distance_to(to);
  }

  ////
  //// ResolveAmbiguousArcsGreedy
  ////
  HypothesesGraph& ResolveAmbiguousArcsGreedy::operator()(HypothesesGraph* g) {
    
    return *g;
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
      LOG(logDEBUG3) << "MergerResolver::add_arcs_for_replacement_node(): added incoming arc: (" << from_tr.Id << "," << trax.Id << "), dist=" << dist << ", state=" << arc_active_map[arc];
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
      LOG(logDEBUG3) << "MergerResolver::add_arcs_for_replacement_node(): added outgoing arc " << g_->id(arc)  << " (" << trax.Id << "," << to_tr.Id << "), dist=" << dist << ", state=" << arc_active_map[arc];
    }
  }

  void MergerResolver::deactivate_arcs(std::vector<HypothesesGraph::base_graph::Arc> arcs) {
    // deactivate Arcs provided by arcs
    // useful to deactivate Arcs of merger Node
    std::vector<HypothesesGraph::base_graph::Arc>::iterator it = arcs.begin();
    property_map<arc_active, HypothesesGraph::base_graph>::type& arc_active_map = g_->get(arc_active());
    for (; it != arcs.end(); ++it) {
      //      if (!g_->valid(g_->source(*it)) || !g_->valid(g_->target(*it))) {
      LOG(logDEBUG3) << "MergerResolver::deactivate_arcs(): setting arc " << g_->id(*it)  << " (" << g_->get(node_traxel())[g_->source((*it))].Id << "," << g_->get(node_traxel())[g_->target((*it))].Id << ") property arc_active to false";
	arc_active_map.set(*it, false);
	//      }
	g_->erase(*it);
    }
  }

  void MergerResolver::deactivate_nodes(std::vector<HypothesesGraph::Node> nodes) {
    // deactivate Nodes provided by nodes
    // needed to set all resolved merger nodes inactive
    std::vector<HypothesesGraph::Node>::iterator it = nodes.begin();
    property_map<node_active2, HypothesesGraph::base_graph>::type& node_active_map = g_->get(node_active2());
    for (; it != nodes.end(); ++it) {
      LOG(logDEBUG3) << "MergerResolver::deactivate_nodes(): setting Node " << g_->id(*it) << " property node_active2 to 0";
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
				   FeatureHandlerBase& handler) {
    // property_map<node_active2, HypothesesGraph::base_graph>::type& active_map = g_->get(node_active2());
    // property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g_->get(node_traxel());
    property_map<node_timestep, HypothesesGraph::base_graph>::type& time_map = g_->get(node_timestep());
    property_map<merger_resolved_to, HypothesesGraph::base_graph>::type& resolved_map = g_->get(merger_resolved_to());
    // property_map<node_originated_from, HypothesesGraph::base_graph>::type& origin_map = g_->get(node_originated_from());
    

    
    
    // Traxel trax = traxel_map[node];
    int timestep = time_map[node];
    unsigned int max_id = get_max_id(timestep)+1;

    // get incoming and outgoing arcs for reorganizing arcs
    std::vector<HypothesesGraph::base_graph::Arc> sources;
    std::vector<HypothesesGraph::base_graph::Arc> targets;
    collect_arcs(HypothesesGraph::base_graph::InArcIt(*g_, node), sources);
    collect_arcs(HypothesesGraph::base_graph::OutArcIt(*g_, node), targets);

    // instead of calling everything with traxels:
    // write functor for extracting and setting features:
    // FeatureHandlerBase ft_handler_(*g, node, nMerger, max_id, sources, targets, vector<unsigned int>& new_ids);


    // create new node for each of the objects merged into node
    std::vector<unsigned int> new_ids;
    handler(*g_, node, nMerger, max_id, timestep, sources, targets, new_ids);
    // std::vector<Traxel> ft = extractor(trax, nMerger, max_id);
    // LOG(logDEBUG) << "MergerResolver::refine_node(): Resolving node " << g_->id(node) << " (" << trax << ") to " << active_map[node] << " merger(s)";
    /* for (std::vector<Traxel>::iterator it = ft.begin(); it != ft.end(); ++it) {
      // set traxel features, most of which can be copied from the merger node
      // set new center of mass as calculated from GMM
      // add node to graph and activate it
      HypothesesGraph::Node newNode = g_->add_node(timestep);
      LOG(logDEBUG3) << "MergerResolver::refine_node(): new node " << g_->id(newNode) << " (" << *it << ") at (" << it->features["com"][0] << "," << it->features["com"][1] << "," << it->features["com"][2] << ")";
      // ft_setter(newNode, 
      traxel_map.set(newNode, *it);
      active_map.set(newNode, 1);
      time_map.set(newNode, timestep);
      // add arc candidates for new nodes (todo: need to somehow choose which ones are active)
      add_arcs_for_replacement_node(newNode, *it, sources, targets, distance);
      // save new id from merger node to new_ids;
      new_ids.push_back(it->Id);
      // store parent (merger) node. this is used for creating the resolved_to event later
      origin_map.set(newNode, std::vector<unsigned int>(1, traxel_map[node].Id));
    } */
    // deactivate incoming and outgoing arcs of merger node
    // merger node will be deactivated after pruning
    deactivate_arcs(sources);
    deactivate_arcs(targets);
    // save information on new ids in property map
    resolved_map.set(node, new_ids);
  }


  HypothesesGraph* MergerResolver::resolve_mergers(FeatureHandlerBase& handler) {
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
	  refine_node(active_itemIt, *active_valueIt, handler);
	  nodes_to_deactivate.push_back(active_itemIt);
	}
      }
    }
    // maybe keep merger nodes active for event extraction
    deactivate_nodes(nodes_to_deactivate);
    prune_inactive(*g_);
    return g_;
  }


  void resolve_graph(HypothesesGraph& src, HypothesesGraph& dest) {
    if (!dest.has_property(node_traxel())) {
      dest.add(node_traxel());
    }
    if (!dest.has_property(node_tracklet())) {
      dest.add(node_tracklet());
    }
    if (!dest.has_property(tracklet_intern_dist())) {
      dest.add(tracklet_intern_dist());
    }
    if (!dest.has_property(arc_distance())) {
      dest.add(arc_distance());
    }
    if (!dest.has_property(node_originated_from())) {
      dest.add(node_originated_from());
    }

    /* HypothesesGraph::NodeMap<HypothesesGraph::Node> nr(dummy_sub);
       HypothesesGraph::ArcMap<HypothesesGraph::Arc> ar(dummy_sub);
       HypothesesGraph::NodeMap<HypothesesGraph::Node> ncr(dummy_sub);
       HypothesesGraph::ArcMap<HypothesesGraph::Arc> acr(dummy_sub); */ 
    /* SubResolver::NodeMap<HypothesesGraph::base_graph::Node> nr(dummy_sub);
       SubResolver::ArcMap<HypothesesGraph::base_graph::Arc> ar(dummy_sub);
       HypothesesGraph::base_graph::NodeMap<SubResolver::Node> ncr(dummy_sub);
       HypothesesGraph::base_graph::ArcMap<SubResolver::Arc> acr(dummy_sub); */
    std::map<HypothesesGraph::Node, HypothesesGraph::Node> nr;
    std::map<HypothesesGraph::Arc, HypothesesGraph::Arc> ar;
    std::map<HypothesesGraph::Node, HypothesesGraph::Node> ncr;
    std::map<HypothesesGraph::Arc, HypothesesGraph::Arc> acr;
    copy_hypotheses_graph_subset<node_resolution_candidate, arc_resolution_candidate>(src, dest, nr, ar, ncr, acr);
    std::vector<double> prob;
    prob.push_back(0.0);
    prob.push_back(1.0);
    boost::function<double(const Traxel&, const size_t)> division = NegLnDivision(1); // weight 1
    boost::function<double(const double)> transition = NegLnTransition(1); // weight 1
    boost::function<double(const Traxel&, const size_t)> detection = boost::bind<double>(NegLnConstant(1,prob), _2);
    translate_property_value_map<node_traxel, HypothesesGraph::Node>(src, dest, nr);
    translate_property_value_map<arc_distance, HypothesesGraph::Arc>(src, dest, ar);
    translate_property_value_map<node_originated_from, HypothesesGraph::Node>(src, dest, nr);

    ConservationTracking pgm(
                             1, //max_number_objects_,
                             detection, //detection,
                             division, // division
                             transition, // transition
                             0, // forbidden_cost_,
                             0.05, // ep_gap_,
                             true, // with_tracklets_,
                             false, // with_divisions_,
                             0, // disappearance_cost_
                             0, // appearance_cost
                             false, // with_misdetections_allowed
                             false, // with appearance
                             false // with disappearance
                             );
                             
    /*
      false, //with_appearance_,
      false, //with_disappearance_,
      false //with_tracklets_
      );*/
    pgm.formulate(dest);
    pgm.infer();
    pgm.conclude(dest);

    translate_property_bool_map<arc_active, HypothesesGraph::Arc>(dest, src, acr);
  }


  double calculate_BIC(int k, int n_samples, double weight, const GMM& gmm) {
    return gmm.score()/n_samples - weight*k; 
  }


  void gmm_priors_and_centers(const feature_array& data, feature_array& priors, feature_array& centers, int k_max, int n, double weight, double n_zero_prob=0.1) {
    assert(priors.size() == 0);
    assert(centers.size() == 0);
    assert(n_zero_prob > 0.0);
    assert(n_zero_prob < 1.0);
    priors.push_back(n_zero_prob);
    double n_not_zero_prob = 1 - n_zero_prob;
    std::back_insert_iterator<feature_array > push_back_iterator(centers);
    int n_samples = data.size()/n;
    double min_bic = 0.0;
    double partition_function = 0.0;
    for (int k = 1; k <= k_max; ++k) {
      GMM gmm(k, n, data);
      feature_array means = gmm();
      double curr_bic = calculate_BIC(k, n_samples, weight, gmm);
      if (curr_bic < min_bic) {
        min_bic = curr_bic;
      }
      priors.push_back(curr_bic);
      partition_function += curr_bic;
      std::copy(means.begin(), means.end(), push_back_iterator);
    }
    partition_function -= k_max*min_bic;
    for (feature_array::iterator it = priors.begin(); it != priors.end(); ++it) {
      *it = n_not_zero_prob*(*it - min_bic)/partition_function;
    }
  }

}



