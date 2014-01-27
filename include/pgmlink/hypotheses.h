/**
   @file
   @ingroup matching
   @brief Hypotheses Graph
*/

#ifndef HYPOTHESES_H
#define HYPOTHESES_H

#include <vector>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <boost/serialization/set.hpp>
#include <boost/shared_ptr.hpp>
#include <lemon/list_graph.h>
#include <lemon/maps.h>

#include "pgmlink/event.h"
#include "pgmlink/graph.h"
#include "pgmlink/log.h"
#include "pgmlink/pgmlink_export.h"
#include "pgmlink/traxels.h"

namespace pgmlink {

////
//// IterableEditableValueMap
////
template <typename Graph, typename Key, typename Value>
class IterableEditableValueMap : public lemon::IterableValueMap<Graph, Key, Value> {
 private:
 public:
  Value& get_value(const Key& key);
  explicit IterableEditableValueMap(const Graph& graph,
                                    const Value& value = Value());
};

template <typename Graph, typename Key, typename Value>
IterableEditableValueMap<Graph, Key, Value>::IterableEditableValueMap(const Graph& graph,
                                                                      const Value& value) :
    lemon::IterableValueMap<Graph,Key, Value>(graph, value) {
    
}

template <typename Graph, typename Key, typename Value>
Value& IterableEditableValueMap<Graph, Key, Value>::get_value(const Key& key) {
  return lemon::IterableValueMap<Graph, Key, Value>::Parent::operator[](key).value;
}


  ////
  //// HypothesesGraph
  ////

  // Properties of a HypothesesGraph

  // node_timestep
  struct node_timestep {};
  template <typename Graph>
    struct property_map<node_timestep, Graph> {
    typedef lemon::IterableValueMap< Graph, typename Graph::Node, int > type;
    static const std::string name;
  };
  template <typename Graph>
    const std::string property_map<node_timestep,Graph>::name = "node_timestep";

  // node_traxel
  struct node_traxel {};
  class Traxel;
  template <typename Graph>
    struct property_map<node_traxel, Graph> {
    typedef IterableEditableValueMap< Graph, typename Graph::Node, Traxel > type;
    static const std::string name;
  };
  template <typename Graph>
    const std::string property_map<node_traxel,Graph>::name = "node_traxel";

  // node_traxel
	struct node_tracklet {};
	template <typename Graph>
	  struct property_map<node_tracklet, Graph> {
	  typedef lemon::IterableValueMap< Graph, typename Graph::Node, std::vector<Traxel> > type;
	  static const std::string name;
	};
	template <typename Graph>
	  const std::string property_map<node_tracklet,Graph>::name = "node_tracklet";


	// tracklet_arcs
	struct tracklet_intern_dist {};
	template <typename Graph>
	  struct property_map<tracklet_intern_dist, Graph> {
	  typedef lemon::IterableValueMap< Graph, typename Graph::Node, std::vector<double> > type;
	  static const std::string name;
	};
	template <typename Graph>
	  const std::string property_map<tracklet_intern_dist,Graph>::name = "tracklet_intern_dist";

	// tracklet_arcs
	struct tracklet_intern_arc_ids {};
	template <typename Graph>
	  struct property_map<tracklet_intern_arc_ids, Graph> {
	  typedef lemon::IterableValueMap< Graph, typename Graph::Node, std::vector<int> > type;
	  static const std::string name;
	};
	template <typename Graph>
	  const std::string property_map<tracklet_intern_arc_ids,Graph>::name = "tracklet_intern_arc_ids";

  // node_active
  struct node_active {};
  template <typename Graph>
    struct property_map<node_active, Graph> {
    typedef lemon::IterableBoolMap< Graph, typename Graph::Node> type;
    static const std::string name;
  };
  template <typename Graph>
    const std::string property_map<node_active,Graph>::name = "node_active";

  // node_active2
    struct node_active2 {};
    template <typename Graph>
      struct property_map<node_active2, Graph> {
      typedef lemon::IterableValueMap< Graph, typename Graph::Node, std::size_t> type;
      static const std::string name;
    };
    template <typename Graph>
      const std::string property_map<node_active2,Graph>::name = "node_active2";

  // node_offered
  struct node_offered {};
  template <typename Graph>
    struct property_map<node_offered, Graph> {
    typedef lemon::IterableBoolMap< Graph, typename Graph::Node> type;
    static const std::string name;
  };
  template <typename Graph>
    const std::string property_map<node_offered,Graph>::name = "node_offered";

  // arc_distance
  struct arc_distance {};
  template <typename Graph>
    struct property_map<arc_distance, Graph> {
    typedef lemon::IterableValueMap< Graph, typename Graph::Arc, double> type;
    static const std::string name;
  };
  template <typename Graph>
    const std::string property_map<arc_distance,Graph>::name = "arc_distance";

  // traxel_arc_id
  struct traxel_arc_id {};
    template <typename Graph>
      struct property_map<traxel_arc_id, Graph> {
      typedef lemon::IterableValueMap< Graph, typename Graph::Arc, int> type;
      static const std::string name;
    };
    template <typename Graph>
      const std::string property_map<traxel_arc_id,Graph>::name = "traxel_arc_id";

  struct arc_vol_ratio {};
    template <typename Graph>
      struct property_map<arc_vol_ratio, Graph> {
      typedef lemon::IterableValueMap< Graph, typename Graph::Arc, double> type;
      static const std::string name;
    };
    template <typename Graph>
      const std::string property_map<arc_vol_ratio,Graph>::name = "arc_vol_ratio";

  // split_into
  struct split_from {};
  template <typename Graph>
    struct property_map<split_from, Graph> {
    typedef lemon::IterableValueMap< Graph, typename Graph::Node, int> type;
    static const std::string name;
  };
  template <typename Graph>
    const std::string property_map<split_from,Graph>::name = "split_from";

  // arc_from_timestep
  struct arc_from_timestep {};
  template <typename Graph>
    struct property_map<arc_from_timestep, Graph> {
    typedef typename Graph::template ArcMap<int> type;
    static const std::string name;
  };
  template <typename Graph>
    const std::string property_map<arc_from_timestep,Graph>::name = "arc_from_timestep";

  // arc_to_timestep
  struct arc_to_timestep {};
  template <typename Graph>
    struct property_map<arc_to_timestep, Graph> {
    typedef typename Graph::template ArcMap<int> type;
    static const std::string name;
  };
  template <typename Graph>
    const std::string property_map<arc_to_timestep,Graph>::name = "arc_to_timestep";

  // arc_active
  struct arc_active {};
  template <typename Graph>
    struct property_map<arc_active, Graph> {
    typedef lemon::IterableBoolMap< Graph, typename Graph::Arc> type;
    static const std::string name;
  };
  template <typename Graph>
    const std::string property_map<arc_active,Graph>::name = "arc_active";

  // division_active
    struct division_active {};
    template <typename Graph>
      struct property_map<division_active, Graph> {
      typedef lemon::IterableBoolMap< Graph, typename Graph::Node> type;
      static const std::string name;
    };
    template <typename Graph>
      const std::string property_map<division_active,Graph>::name = "division_active";

  // merger_resolved_to
  struct merger_resolved_to {};
  template <typename Graph>
  struct property_map<merger_resolved_to, Graph> {
    // typedef std::map<typename Graph::Node, std::vector<unsigned int> > type;
    typedef lemon::IterableValueMap<Graph, typename Graph::Node, std::vector<unsigned int> > type;
    static const std::string name;
  };
  template <typename Graph>
  const std::string property_map<merger_resolved_to, Graph>::name = "merger_resolved_to";

  // node_originated_from
  struct node_originated_from {};
  template <typename Graph>
  struct property_map<node_originated_from, Graph> {
    typedef lemon::IterableValueMap<Graph, typename Graph::Node, std::vector<unsigned int> > type;
    static const std::string name;
  };
  template <typename Graph>
  const std::string property_map<node_originated_from, Graph>::name = "node_originated_from";

  // node_resolution_candidate
  struct node_resolution_candidate {};
  template <typename Graph>
  struct property_map<node_resolution_candidate, Graph> {
    typedef lemon::IterableBoolMap<Graph, typename Graph::Node> type;
    static const std::string name;
  };
  template <typename Graph>
  const std::string property_map<node_resolution_candidate, Graph>::name = "node_resolution_candidate";

  // arc_resolution_candidate
  struct arc_resolution_candidate {};
  template <typename Graph>
  struct property_map<arc_resolution_candidate, Graph> {
    typedef lemon::IterableBoolMap<Graph, typename Graph::Arc> type;
    static const std::string name;
  };
  template <typename Graph>
  const std::string property_map<arc_resolution_candidate, Graph>::name = "arc_resolution_candidate";



  class PGMLINK_EXPORT HypothesesGraph : public PropertyGraph<lemon::ListDigraph> {
  public:
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map;

    HypothesesGraph() {
      // Properties attached to every HypothesesGraph
      add(node_timestep());
      add(arc_from_timestep());
      add(arc_to_timestep());
    };

    // use this instead of calling the parent graph directly
    HypothesesGraph::Node add_node(node_timestep_map::Value timestep);
    // call this function to add a multi-temporal node (e.g. for tracklets)
    HypothesesGraph::Node add_node(std::vector<node_timestep_map::Value> timesteps);

    const std::set<HypothesesGraph::node_timestep_map::Value>& timesteps() const;
    node_timestep_map::Value earliest_timestep() const;
    node_timestep_map::Value latest_timestep() const;
    
  protected:
    std::set<node_timestep_map::Value> timesteps_;
  private:
    // boost serialize
    friend class boost::serialization::access;
    template< typename Archive >
      void save( Archive&, const unsigned int /*version*/ ) const;
    template< typename Archive >
      void load( Archive&, const unsigned int /*version*/ );
    BOOST_SERIALIZATION_SPLIT_MEMBER()
  };

  void generateTrackletGraph(const HypothesesGraph& traxel_graph, HypothesesGraph& tracklet_graph);
  std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > generateTrackletGraph2(
		  const HypothesesGraph& traxel_graph, HypothesesGraph& tracklet_graph);
  PGMLINK_EXPORT HypothesesGraph& prune_inactive(HypothesesGraph&);
  PGMLINK_EXPORT boost::shared_ptr<std::vector< std::vector<Event> > > events(const HypothesesGraph&);
  PGMLINK_EXPORT boost::shared_ptr<std::vector< std::vector<Event> > > multi_frame_move_events(const HypothesesGraph& g);
  PGMLINK_EXPORT boost::shared_ptr<std::vector< std::vector<Event> > > merge_event_vectors(const std::vector<std::vector<Event> >& ev1, const std::vector<std::vector<Event> >& ev2);
  PGMLINK_EXPORT boost::shared_ptr<std::vector< std::map<unsigned int, bool> > > state_of_nodes(const HypothesesGraph&);

  // lemon graph format (lgf) serialization
  PGMLINK_EXPORT void write_lgf( const HypothesesGraph&, std::ostream& os=std::cout,
		  bool with_n_traxel=false );
  PGMLINK_EXPORT void read_lgf( HypothesesGraph&, std::istream& is=std::cin,
		 bool with_n_traxel=false);



  ////
  //// HypothesesBuilder
  ////
  class PGMLINK_EXPORT HypothesesBuilder {
  public:
    virtual HypothesesGraph* build() const;

  protected:
    // template methods
    virtual HypothesesGraph* construct() const = 0;
    virtual HypothesesGraph* add_nodes(HypothesesGraph*) const = 0;
    virtual HypothesesGraph* add_edges(HypothesesGraph*) const = 0;
  };



  ////
  //// SingleTimestepTraxel_HypothesesBuilder
  ////
  class PGMLINK_EXPORT SingleTimestepTraxel_HypothesesBuilder : public HypothesesBuilder {
  public:
    struct Options {
	Options(unsigned int mnn = 6, double dt = 50,
			bool forward_backward=false, bool consider_divisions=false,
			double division_threshold = 0.5) :
  		max_nearest_neighbors(mnn), distance_threshold(dt), forward_backward(forward_backward),
  		consider_divisions(consider_divisions),
  		division_threshold(division_threshold){};
  	unsigned int max_nearest_neighbors;
  	double distance_threshold;
  	bool forward_backward, consider_divisions;
  	double division_threshold;
    };

  SingleTimestepTraxel_HypothesesBuilder(const TraxelStore* ts, const Options& o = Options()) : ts_(ts), options_(o) {};

  protected:
    // builder method implementations
    virtual HypothesesGraph* construct() const;
    virtual HypothesesGraph* add_nodes(HypothesesGraph*) const;
    virtual HypothesesGraph* add_edges(HypothesesGraph*) const;

    const TraxelStore* ts_;
    Options options_;
  private:
    HypothesesGraph* add_edges_at(HypothesesGraph*, int timestep, bool reverse=false) const;
  };




  /**/
  /* implementation */
  /**/
  template< typename Archive >
    void HypothesesGraph::save( Archive& ar, const unsigned int /*version*/ ) const {
    ar & timesteps_;
    
    bool with_n_traxel = false;
    try {
      this->get(node_traxel());
      with_n_traxel = true;
    } catch( std::runtime_error& e ) {}
    ar & with_n_traxel;

    std::string lgf;
    {
      std::stringstream ss;
      write_lgf(*this, ss, with_n_traxel);
      lgf = ss.str();
    }
    ar & lgf;
  }

  template< typename Archive >
    void HypothesesGraph::load( Archive& ar, const unsigned int /*version*/ ) {
    ar & timesteps_;

    bool with_n_traxel;
    ar & with_n_traxel;

    std::string lgf;
    ar & lgf;
    {
      std::stringstream ss(lgf);
      read_lgf(*this, ss, with_n_traxel);
    }
   }

}
#endif /* HYPOTHESES_H */
