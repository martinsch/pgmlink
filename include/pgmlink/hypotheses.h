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
    typedef lemon::IterableValueMap< Graph, typename Graph::Node, Traxel > type;
    static const std::string name;
  };
  template <typename Graph>
    const std::string property_map<node_traxel,Graph>::name = "node_traxel";

  // node_active
  struct node_active {};
  template <typename Graph>
    struct property_map<node_active, Graph> {
    typedef lemon::IterableBoolMap< Graph, typename Graph::Node> type;
    static const std::string name;
  };
  template <typename Graph>
    const std::string property_map<node_active,Graph>::name = "node_active";

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

    const std::set<HypothesesGraph::node_timestep_map::Value>& timesteps() const;
    node_timestep_map::Value earliest_timestep() const;
    node_timestep_map::Value latest_timestep() const;
    
  private:
    // boost serialize
    friend class boost::serialization::access;
    template< typename Archive >
      void save( Archive&, const unsigned int /*version*/ ) const;
    template< typename Archive >
      void load( Archive&, const unsigned int /*version*/ );
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    std::set<node_timestep_map::Value> timesteps_;      
  };

  PGMLINK_EXPORT HypothesesGraph& prune_inactive(HypothesesGraph&);
  PGMLINK_EXPORT boost::shared_ptr<std::vector< std::vector<Event> > > events(const HypothesesGraph&);
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
	Options(unsigned int mnn = 6, double dt = 50) : max_nearest_neighbors(mnn),
							distance_threshold(dt) {};
	unsigned int max_nearest_neighbors;
	double distance_threshold;
    };

  SingleTimestepTraxel_HypothesesBuilder(const TraxelStore* ts, const Options& o = Options()) : ts_(ts), options_(o) {};

  protected:
    // builder method implementations
    virtual HypothesesGraph* construct() const;
    virtual HypothesesGraph* add_nodes(HypothesesGraph*) const;
    virtual HypothesesGraph* add_edges(HypothesesGraph*) const;

  private:
    HypothesesGraph* add_edges_at(HypothesesGraph*, int timestep) const;

    const TraxelStore* ts_;
    Options options_;
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
