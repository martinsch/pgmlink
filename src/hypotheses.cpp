#include <cassert>
#include <utility>
#include <vector>
#include <boost/tuple/tuple.hpp>
#include "hypotheses.h"
#include "log.h"
#include "nearest_neighbors.h"

using namespace std;

namespace Tracking {
  ////
  //// class HypothesesGraph
  ////
  HypothesesGraph::Node HypothesesGraph::add_node(node_timestep_map::Value timestep) {
      node_timestep_map& timestep_m = get(node_timestep());

      HypothesesGraph::Node node = addNode();
      timestep_m.set(node, timestep);
      timesteps_.insert( timestep );
      return node;
  }

  const std::set<HypothesesGraph::node_timestep_map::Value>& HypothesesGraph::timesteps() const {
    return timesteps_;
  }

  HypothesesGraph::node_timestep_map::Value HypothesesGraph::earliest_timestep() const {
    return *(timesteps_.begin());
  }

  HypothesesGraph::node_timestep_map::Value HypothesesGraph::latest_timestep() const {
    return *(timesteps_.rbegin());
  }


    
  HypothesesGraph& prune_inactive(HypothesesGraph& g) {
      LOG(logDEBUG) << "prune_inactive(): entered";
      property_map<node_active, HypothesesGraph::base_graph>::type& active_nodes = g.get(node_active());
      property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = g.get(arc_active());
    
      // prune inactive arcs
      LOG(logDEBUG) << "prune_inactive: prune inactive arcs";
      typedef property_map<arc_active, HypothesesGraph::base_graph>::type::FalseIt inactive_arc_it;

      // we first collect and then erase (iterator will be invalid after erase; therefore a two-step
      // procedure)

      // collect inactive arcs
      vector<HypothesesGraph::Arc> arcs_to_prune;

      for(inactive_arc_it it(active_arcs); it!=lemon::INVALID; ++it) {
	arcs_to_prune.push_back(it);
      } 

      // prune inactive arcs
      for(vector<HypothesesGraph::Arc>::const_iterator it = arcs_to_prune.begin(); it!= arcs_to_prune.end(); ++it) {
	LOG(logDEBUG3) << "prune_inactive: pruned arc: " << g.id(*it);
	g.erase(*it);
      }

      // prune inactive nodes 
      LOG(logDEBUG) << "prune_inactive: prune inactive nodes";
      typedef property_map<node_active, HypothesesGraph::base_graph>::type::FalseIt inactive_node_it;
      // collect inactive nodes
      vector<HypothesesGraph::Node> nodes_to_prune;

      for(inactive_node_it it(active_nodes); it!=lemon::INVALID; ++it) {
	nodes_to_prune.push_back(it);
      } 

      // prune inactive nodes
      for(vector<HypothesesGraph::Node>::const_iterator it = nodes_to_prune.begin(); it!= nodes_to_prune.end(); ++it) {
	LOG(logDEBUG3) << "prune_inactive: pruned node: " << g.id(*it);
	g.erase(*it);
      } 

      return g;
  }


  boost::shared_ptr<std::vector< std::vector<Event> > > events(const HypothesesGraph& g) {
    LOG(logDEBUG) << "events(): entered";
    shared_ptr<std::vector< std::vector<Event> > > ret(new vector< vector<Event> >);
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    node_timestep_map_t& node_timestep_map = g.get(node_timestep());
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map_t;
    node_traxel_map_t& node_traxel_map = g.get(node_traxel());



    // for every timestep
    LOG(logDEBUG1) << "events(): earliest_timestep: " << g.earliest_timestep();
    LOG(logDEBUG1) << "events(): latest_timestep: " << g.latest_timestep();
    for(int t = g.earliest_timestep(); t < g.latest_timestep(); ++t) {
        LOG(logDEBUG2) << "events(): processing timestep: " << t;
	ret->push_back(vector<Event>());

	// for every node: destiny
	LOG(logDEBUG2) << "events(): for every node: destiny";
	for(node_timestep_map_t::ItemIt node_at(node_timestep_map, t); node_at!=lemon::INVALID; ++node_at) {
	    assert(node_traxel_map[node_at].Timestep == t);
	    // count ougoing arcs
	    int count = 0;
	    for(HypothesesGraph::base_graph::OutArcIt a(g, node_at); a!=lemon::INVALID; ++a) ++count;
	    LOG(logDEBUG3) << "events(): counted outgoing arcs: " << count;
	    // construct suitable Event object
	    switch(count) {
		// Disappearance
		case 0: {
		    Event e;
		    e.type = Event::Disappearance;
		    e.traxel_ids.push_back(node_traxel_map[node_at].Id);
		    (*ret)[t-g.earliest_timestep()].push_back(e);
		    LOG(logDEBUG3) << e;
		    break;
		    }
		// Move
		case 1: {
		    Event e;
		    e.type = Event::Move;
		    e.traxel_ids.push_back(node_traxel_map[node_at].Id);
		    HypothesesGraph::base_graph::OutArcIt a(g, node_at);
		    e.traxel_ids.push_back(node_traxel_map[g.target(a)].Id);
		    (*ret)[t-g.earliest_timestep()].push_back(e);
		    LOG(logDEBUG3) << e;
		    break;
		    }
		// Division
		case 2: {
		    Event e;
		    e.type = Event::Division;
		    e.traxel_ids.push_back(node_traxel_map[node_at].Id);
		    HypothesesGraph::base_graph::OutArcIt a(g, node_at);
		    e.traxel_ids.push_back(node_traxel_map[g.target(a)].Id);
		    ++a;
		    e.traxel_ids.push_back(node_traxel_map[g.target(a)].Id);
		    (*ret)[t-g.earliest_timestep()].push_back(e);
		    LOG(logDEBUG3) << e;
		break;
	        }
		default:
		    throw runtime_error("events(): encountered node dividing in three or more nodes in graph");
		break;
	    }
	}

	// appearances in next timestep
	LOG(logDEBUG2) << "events(): appearances in next timestep";
	for(node_timestep_map_t::ItemIt node_at(node_timestep_map, t+1); node_at!=lemon::INVALID; ++node_at) {
	    // count incoming arcs
	    int count = 0;
	    for(HypothesesGraph::base_graph::InArcIt a(g, node_at); a!=lemon::INVALID; ++a) ++count;
    	    LOG(logDEBUG3) << "events(): counted incoming arcs in next timestep: " << count;
	   
	    // no incoming arcs => appearance
	    if(count == 0) {
		    Event e;
		    e.type = Event::Appearance;
		    e.traxel_ids.push_back(node_traxel_map[node_at].Id);
		    (*ret)[t-g.earliest_timestep()].push_back(e);
		    LOG(logDEBUG3) << e;	      
	    }
	}
    }

    return ret;
  }



  //
  // state_of_nodes()
  //
  boost::shared_ptr<std::vector< std::map<unsigned int, bool> > > state_of_nodes(const HypothesesGraph& g) {
    LOG(logDEBUG) << "detections(): entered";
    shared_ptr<vector< map<unsigned int, bool> > > ret(new vector< map<unsigned int, bool> >);

    // required node properties: timestep, traxel, active
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    node_timestep_map_t& node_timestep_map = g.get(node_timestep());
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map_t;
    node_traxel_map_t& node_traxel_map = g.get(node_traxel());
    typedef property_map<node_active, HypothesesGraph::base_graph>::type node_active_map_t;
    node_active_map_t& node_active_map = g.get(node_active());


    // for every timestep
    for(int t = g.earliest_timestep(); t <= g.latest_timestep(); ++t) {
      ret->push_back(map<unsigned int, bool>());
      for(node_timestep_map_t::ItemIt node_at(node_timestep_map, t); node_at!=lemon::INVALID; ++node_at) {
	assert(node_traxel_map[node_at].Timestep == t);
	unsigned int id = node_traxel_map[node_at].Id;
	bool active = node_active_map[node_at];
	(*ret)[t-g.earliest_timestep()][id] = active;
      }
    }

    return ret;
  }
    


  ////
  //// class HypothesesBuilder
  ////
  HypothesesGraph* HypothesesBuilder::build() const {
    // construct an empty HypothesesGraph with all desired additional 
    // properties added
    HypothesesGraph* graph = construct();
    // add object nodes and set node properties
    graph = add_nodes(graph);
    // connect object nodes and set edge properties
    graph = add_edges(graph);

    return graph;
  }



  ////
  //// class SingleTimestepTraxel_HypothesesBuilder
  ////
  HypothesesGraph* SingleTimestepTraxel_HypothesesBuilder::construct() const {
    HypothesesGraph* graph = new HypothesesGraph();
    // store traxels inside the graph data structure
    graph->add(node_traxel());
    return graph;
  }

  HypothesesGraph* SingleTimestepTraxel_HypothesesBuilder::add_nodes(HypothesesGraph* graph) const {
    LOG(logDEBUG) << "SingleTimestepTraxel_HypothesesBuilder::add_nodes(): entered";
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_m = graph->get(node_traxel());

    for(TraxelStoreByTimestep::const_iterator it = ts_->begin(); it!= ts_->end(); ++it) {
      HypothesesGraph::Node node = graph->add_node(it->Timestep);
      traxel_m.set(node, *it);
    }

    return graph;
  }

  HypothesesGraph* SingleTimestepTraxel_HypothesesBuilder::add_edges(HypothesesGraph* graph) const {
    typedef HypothesesGraph::node_timestep_map::Value timestep_t;
    const set<timestep_t>& timesteps = graph->timesteps();
    // iterate over all timesteps except the last
    for(set<timestep_t>::const_iterator t = timesteps.begin(); t != (--timesteps.end()); ++t) {
	add_edges_at(graph, *t);
    }

    return graph;
  }

  HypothesesGraph* SingleTimestepTraxel_HypothesesBuilder::add_edges_at(HypothesesGraph* graph, int timestep) const {
      const HypothesesGraph::node_timestep_map& timemap = graph->get(node_timestep());
      typedef property_map<node_traxel, HypothesesGraph::base_graph>::type traxelmap_t;
      const traxelmap_t& traxelmap = graph->get(node_traxel());
      const TraxelStoreByTimeid& traxels_by_timeid = ts_->get<by_timeid>();
      const TraxelStoreByTimestep& traxels_by_timestep = ts_->get<by_timestep>();

      // establish transition edges between a current node and appropriate nodes in next timestep
      for(HypothesesGraph::node_timestep_map::ItemIt curr_node(timemap, timestep); curr_node!=lemon::INVALID; ++curr_node) {
	  assert(timemap[curr_node] == timestep);
	  assert(traxelmap[curr_node].Timestep == timestep);

	  //// find k nearest neighbors in next timestep
	  // init nearest neighbor search
	  pair<TraxelStoreByTimestep::const_iterator,TraxelStoreByTimestep::const_iterator> traxels_at =
	    traxels_by_timestep.equal_range(timestep+1);

	  for(TraxelStoreByTimestep::const_iterator it = traxels_at.first; it != traxels_at.second; ++it) {
	    assert(it->Timestep == (timestep+1));
	  }

	  NearestNeighborSearch nns(traxels_at.first, traxels_at.second);

	  // search
	  map<unsigned int, double> nearest_neighbors =
	    nns.knn_in_range(traxelmap[curr_node], options_.distance_threshold, options_.max_nearest_neighbors);
	  
	  //// connect current node with k nearest neighbor nodes
	  for(map<unsigned int, double>::const_iterator neighbor = nearest_neighbors.begin();
	    neighbor != nearest_neighbors.end(); 
            ++neighbor) {
		// connect with one of the neighbor nodes
		TraxelStoreByTimeid::iterator neighbor_traxel = traxels_by_timeid.find(boost::make_tuple(timestep+1, neighbor->first));
		assert(neighbor_traxel->Timestep == (timestep + 1));
		assert(neighbor_traxel->Timestep != traxelmap[curr_node].Timestep);
		traxelmap_t::ItemIt neighbor_node(traxelmap, *neighbor_traxel);
		assert(curr_node != neighbor_node);
		graph->addArc(curr_node, neighbor_node);
	  }
      }
	
      return graph;
  }

} /* namespace Tracking */
