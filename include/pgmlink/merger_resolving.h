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


// pgmlink headers
#include "pgmlink/hypotheses.h"
#include "pgmlink/event.h"
#include "pgmlink/traxels.h"
#include "pgmlink/reasoner.h"
#include "pgmlink/merger_resolving_grammar.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/feature.h"
#include "pgmlink/clustering.h"

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
  //// FeatureExtractorMCOMsFromGMM
  ////
  class FeatureExtractorMCOMsFromGMM : public FeatureExtractorBase {
  private:
    int n_dim_;
  public:
    FeatureExtractorMCOMsFromGMM(int n_dim) : n_dim_(n_dim) {}
    virtual std::vector<Traxel> operator()(Traxel trax, size_t nMergers, unsigned int max_id);
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
  class DistanceFromCOMs : public DistanceBase {
  public:
    virtual double operator()(const HypothesesGraph& g, HypothesesGraph::Node from, HypothesesGraph::Node to);
    virtual double operator()(Traxel from,  Traxel to);
  };

  
  ////
  //// FeatureHandlerBase
  ////
  class FeatureHandlerBase {
  public:
    void add_arcs_for_replacement_node(HypothesesGraph& g,
                                       HypothesesGraph::Node n,
                                       const std::vector<HypothesesGraph::base_graph::Arc>& sources,
                                       const std::vector<HypothesesGraph::base_graph::Arc>& targets,
                                       DistanceBase& distance);
  public:
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
  class FeatureHandlerFromTraxels : public FeatureHandlerBase {
  private:
    FeatureExtractorBase& extractor_;
    DistanceBase& base_;
    // FeatureHandlerFromTraxelsMCOMsFromPCOMs() {};
  public:
    FeatureHandlerFromTraxels(FeatureExtractorBase& extractor,
                                            DistanceBase& base) :
      extractor_(extractor), base_(base) {}

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
  class ResolveAmbiguousArcsGreedy : public ResolveAmbiguousArcsBase {
  public:
    virtual HypothesesGraph& operator()(HypothesesGraph* g);
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
  //// ResolveAmbiguousArcsPgm
  ////
  class ResolveAmbiguousArcsPgm : public ReasonerMaxOneArc, private ResolveAmbiguousArcsBase {
    virtual HypothesesGraph& operator()(HypothesesGraph* g);
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
    void add_arcs_for_replacement_node(HypothesesGraph::Node node,
				       Traxel trax,
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
    MergerResolver(HypothesesGraph* g) : g_(g)
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
    HypothesesGraph* resolve_mergers(FeatureHandlerBase& handler);
  };


  ////
  //// given a graph, do retracking
  ////
  void resolve_graph(HypothesesGraph& src, HypothesesGraph& dest, boost::function<double(const double)> transition, double ep_gap, bool with_tracklets,
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

  
  double calculate_BIC(int k, int n, double weight);


  void gmm_priors_and_centers(const feature_array& data, feature_array& priors, feature_array& centers, int k_max, int n, double weight);

  void gmm_priors_and_centers_arma(const arma::mat& data, feature_array& priors, feature_array& centers, int k_max, int ndim, double regularization_weight);


  ////
  //// duplicate division nodes in subset graph
  ////
  void duplicate_division_nodes(HypothesesGraph& graph,
                                std::map<HypothesesGraph::Node, HypothesesGraph::Node>& division_splits,
                                std::map<HypothesesGraph::Arc, HypothesesGraph::Arc>& arc_cross_reference);


  ////
  //// merge previously split divisions after inference
  ////
  void merge_split_divisions(const HypothesesGraph& graph,
                             std::map<HypothesesGraph::Node, HypothesesGraph::Node>& division_splits,
                             std::map<HypothesesGraph::Arc, HypothesesGraph::Arc>& arc_cross_reference);



  ////
  //// IMPLEMENTATIONS ////
  ////

  
  
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
    
    NodeFilter& node_filter_map = src.get(NodePropertyTag());
    ArcFilter& arc_filter_map = src.get(ArcPropertyTag());
    CopyGraph sub(src, node_filter_map, arc_filter_map);

    lemon::digraphCopy<CopyGraph, HypothesesGraph::base_graph>(sub,dest).nodeRef(nr).nodeCrossRef(ncr).arcRef(ar).arcCrossRef(acr).run();
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
