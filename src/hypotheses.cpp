#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/tuple/tuple.hpp>
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/nearest_neighbors.h"
#include "pgmlink/traxels.h"

using namespace std;

namespace pgmlink {
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

  HypothesesGraph::Node HypothesesGraph::add_node(std::vector<node_timestep_map::Value> timesteps) {
	  node_timestep_map& timestep_m = get(node_timestep());

	  HypothesesGraph::Node node = addNode();
	  for (std::vector<node_timestep_map::Value>::const_iterator it = timesteps.begin(); it!=timesteps.end(); ++it) {
			timestep_m.set(node, *it);
			timesteps_.insert( *it );
	  }
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
      property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = g.get(arc_active());

      property_map<node_active, HypothesesGraph::base_graph>::type* active_nodes;
	  property_map<node_active2, HypothesesGraph::base_graph>::type* active2_nodes;
	  bool active2_used = false;
	  if (g.getProperties().count("node_active") > 0) {
		active_nodes = &g.get(node_active());
	  } else if (g.getProperties().count("node_active2") > 0) {
		active2_nodes = &g.get(node_active2());
		active2_used = true;
	  }
    
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
      // collect inactive nodes
      vector<HypothesesGraph::Node> nodes_to_prune;

      if (active2_used) {
  		for (HypothesesGraph::NodeIt it(g); it != lemon::INVALID; ++it) {
  			if (!(*active2_nodes)[it]) {
  				nodes_to_prune.push_back(it);
  			}
  		}
      } else {
    	typedef property_map<node_active, HypothesesGraph::base_graph>::type::FalseIt inactive_node_it;
		for (inactive_node_it it(*active_nodes); it != lemon::INVALID; ++it) {
			nodes_to_prune.push_back(it);
		}
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
    property_map<division_active, HypothesesGraph::base_graph>::type* division_node_map;
    bool with_division_detection = false;
    if (g.getProperties().count("division_active") > 0) {
    	division_node_map = &g.get(division_active());
    	with_division_detection = true;
    }
    property_map<node_active2, HypothesesGraph::base_graph>::type* node_number_of_objects;
    bool with_mergers = false;
    if (g.getProperties().count("node_active2") > 0) {
    	node_number_of_objects = &g.get(node_active2());
    	with_mergers = true;
    }

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

	    if(with_mergers && (*node_number_of_objects)[node_at] > 1) {
			Event e;
			e.type = Event::Merger;
			e.traxel_ids.push_back(node_traxel_map[node_at].Id);
			e.traxel_ids.push_back((*node_number_of_objects)[node_at]);
			(*ret)[t-g.earliest_timestep()].push_back(e);
			LOG(logDEBUG3) << e;
	    }

	    // count outgoing arcs
	    size_t count = 0;
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
		// Division or Splitting
		default: {
			Event e;
		    if (with_division_detection) {
		    	if (count == 2 && (*division_node_map)[node_at]) {
		    		e.type = Event::Division;
					e.traxel_ids.push_back(node_traxel_map[node_at].Id);
					HypothesesGraph::base_graph::OutArcIt a(g, node_at);
					e.traxel_ids.push_back(node_traxel_map[g.target(a)].Id);
					++a;
					e.traxel_ids.push_back(node_traxel_map[g.target(a)].Id);
					(*ret)[t-g.earliest_timestep()].push_back(e);
					LOG(logDEBUG3) << e;
		    	} else {
					for(HypothesesGraph::base_graph::OutArcIt a(g, node_at); a != lemon::INVALID; ++a) {
						e.type = Event::Move;
						e.traxel_ids.clear();
						e.traxel_ids.push_back(node_traxel_map[node_at].Id);
						e.traxel_ids.push_back(node_traxel_map[g.target(a)].Id);
						(*ret)[t-g.earliest_timestep()].push_back(e);
						LOG(logDEBUG3) << e;
					}
		    	}
		    } else { // for backward compatibility
		    	if (count != 2) {
				    throw runtime_error("events(): encountered node dividing in three or more nodes in graph");
		    	}
				e.type = Event::Division;
				e.traxel_ids.push_back(node_traxel_map[node_at].Id);
				HypothesesGraph::base_graph::OutArcIt a(g, node_at);
				e.traxel_ids.push_back(node_traxel_map[g.target(a)].Id);
				++a;
				e.traxel_ids.push_back(node_traxel_map[g.target(a)].Id);
				(*ret)[t-g.earliest_timestep()].push_back(e);
				LOG(logDEBUG3) << e;
		    }
		break;
	        }
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
//    // mergers in last timestep
//	if(with_mergers) {
//		LOG(logDEBUG2) << "events(): mergers in last timestep";
//		int t = g.latest_timestep();
//		for(node_timestep_map_t::ItemIt node_at(node_timestep_map, t); node_at!=lemon::INVALID; ++node_at) {
//			if((*node_number_of_objects)[node_at] > 1) {
//				Event e;
//				e.type = Event::Merger;
//				e.traxel_ids.push_back(node_traxel_map[node_at].Id);
//				e.traxel_ids.push_back((*node_number_of_objects)[node_at]);
//				(*ret)[t-g.earliest_timestep()].push_back(e);
//				LOG(logDEBUG3) << e;
//			}
//		}
//	}

    return ret;
  }

  //
  // generateTrackletGraph
  //
  // a traxel graph (containing nodes of traxels) with some active arcs is converted into a tracklet graph
  // where all nodes connected by active paths are summarized in a single tracklet node
  void generateTrackletGraph(const HypothesesGraph& traxel_graph, HypothesesGraph& tracklet_graph) {

  	// go through the traxels graph, add each node which doesn't have an active incoming arc, and
  	// follow the active outgoing path to add those nodes to the tracklet
  	property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = traxel_graph.get(arc_active());
	property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = traxel_graph.get(node_traxel());

	property_map<node_tracklet, HypothesesGraph::base_graph>::type* traxels_tracklet_map;
	bool traxel_nodes_are_tracklets = false;
	if (traxel_graph.getProperties().count("node_tracklet") > 0) {
		traxel_nodes_are_tracklets = true;
		traxels_tracklet_map = &traxel_graph.get(node_tracklet());
	}

  	typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
	node_timestep_map_t& node_timestep_map = traxel_graph.get(node_timestep());

	// add empty traxel_map to the tracklet graph in order to make the tracklet graph equivalent to traxelgraphs
	tracklet_graph.add(node_traxel());

  	tracklet_graph.add(node_tracklet());
  	property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = tracklet_graph.get(node_tracklet());

  	// this map stores the traxel_graph nodes at t+1 as keys which will be linked by the
  	// list of tracklet nodes at t stored as map-values
  	std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > node_link_map;

  	// add nodes
  	for(int t = traxel_graph.earliest_timestep(); t <= traxel_graph.latest_timestep(); ++t) {
		for(node_timestep_map_t::ItemIt traxel_node(node_timestep_map, t); traxel_node!=lemon::INVALID; ++traxel_node) {
			bool has_active_incoming = false;
			std::vector<HypothesesGraph::Arc> tracklet_incoming_arcs;
			for(HypothesesGraph::InArcIt a(traxel_graph, traxel_node); a != lemon::INVALID; ++a) {
				tracklet_incoming_arcs.push_back(a);
				if (active_arcs[a]) {
					has_active_incoming = true;
					break;
				}
			}

			if (has_active_incoming) {
				// if the traxel node has an active incoming arc, it has already been added to some tracklet
				continue;
			}
			std::vector<int> timesteps;
			std::vector<Traxel> tracklet;
			std::vector<Traxel> traxels;
			if (!traxel_nodes_are_tracklets) {
				traxels.push_back(traxel_map[traxel_node]);
			} else {
				traxels.insert(traxels.begin(), (*traxels_tracklet_map)[traxel_node].begin(),
						(*traxels_tracklet_map)[traxel_node].end());
			}

			for (std::vector<Traxel>::const_iterator tr = traxels.begin(); tr != traxels.end(); ++tr) {
				timesteps.push_back(tr->Timestep);
				tracklet.push_back(*tr);
			}
			HypothesesGraph::Node tn = traxel_node;
			HypothesesGraph::Node tn_next;
			std::vector<HypothesesGraph::Arc> tn_outarcs;

			// follow the active outgoing path to add those nodes to the tracklet
			while (true) {
				bool active_outgoing = false;
				tn_outarcs.clear();

				for(HypothesesGraph::OutArcIt a(traxel_graph, tn); a != lemon::INVALID; ++a) {
					tn_outarcs.push_back(a);

					if (active_arcs[a]) {
						assert(!active_outgoing); // "found more than one active outgoing arc"
						tn_next = traxel_graph.target(a);

						traxels.clear();
						if (!traxel_nodes_are_tracklets) {
							traxels.push_back(traxel_map[tn_next]);
						} else {
							traxels.insert(traxels.begin(), (*traxels_tracklet_map)[tn_next].begin(),
									(*traxels_tracklet_map)[tn_next].end());
						}

						for (std::vector<Traxel>::const_iterator tr = traxels.begin(); tr != traxels.end(); ++tr) {
							timesteps.push_back(tr->Timestep);
							tracklet.push_back(*tr);
						}
						active_outgoing = true;
					}
				}
				if (active_outgoing) {
					tn = tn_next;
				} else { // no active outgoing arc found -- tracklet ends
					break;
				}
			}

			HypothesesGraph::Node curr_tracklet_node = tracklet_graph.add_node(timesteps);
			tracklet_map.set(curr_tracklet_node, tracklet);

			// store the outgoing arcs of the traxel node
			if (tn_outarcs.size() > 0) {
				for(std::vector<HypothesesGraph::Arc>::const_iterator arc_it = tn_outarcs.begin(); arc_it!=tn_outarcs.end(); ++arc_it) {
					HypothesesGraph::Node traxel_node_to = traxel_graph.target(*arc_it);
					if (node_link_map.find(traxel_node_to) == node_link_map.end()) {
						node_link_map.insert(
								std::pair<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > (
										traxel_node_to, std::vector<HypothesesGraph::Node>() ) );
					}
					node_link_map[traxel_node_to].push_back(curr_tracklet_node);
				}
			}

			// set the incoming arcs of the tracklet (traxel_node is the first node in the tracklet)
			if (!traxel_nodes_are_tracklets) {
				assert(traxel_map[traxel_node].Id==tracklet[0].Id);
				assert(traxel_map[traxel_node].Timestep==tracklet[0].Timestep);
			} else {
				assert((*traxels_tracklet_map)[traxel_node][0].Id==tracklet[0].Id);
				assert((*traxels_tracklet_map)[traxel_node][0].Timestep==tracklet[0].Timestep);
			}
			std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> >::const_iterator traxel_to_it = node_link_map.find(traxel_node);

			if (traxel_to_it == node_link_map.end()) {
				LOG(logDEBUG) << "There is no incoming arc for traxel_node " << traxel_graph.id(traxel_node);
			} else {
				for(std::vector<HypothesesGraph::Node>::const_iterator tracklet_from_it = (*traxel_to_it).second.begin();
						tracklet_from_it != (*traxel_to_it).second.end(); ++tracklet_from_it) {
					tracklet_graph.addArc(*tracklet_from_it,curr_tracklet_node);
				}
			}
		}
  	}
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
    property_map<node_active, HypothesesGraph::base_graph>::type* node_active_map;
    property_map<node_active2, HypothesesGraph::base_graph>::type* node_active2_map;
    bool active2_used = false;
    if (g.getProperties().count("node_active") > 0) {
    	node_active_map = &g.get(node_active());
    } else if (g.getProperties().count("node_active2") > 0) {
    	node_active2_map = &g.get(node_active2());
    	active2_used = true;
    }

    // for every timestep
    for(int t = g.earliest_timestep(); t <= g.latest_timestep(); ++t) {
      ret->push_back(map<unsigned int, bool>());
      for(node_timestep_map_t::ItemIt node_at(node_timestep_map, t); node_at!=lemon::INVALID; ++node_at) {
	assert(node_traxel_map[node_at].Timestep == t);
	unsigned int id = node_traxel_map[node_at].Id;
	bool active = false;
	if (active2_used) {
		active = ((*node_active2_map)[node_at] > 0);
	} else {
		active = (*node_active_map)[node_at];
	}
	(*ret)[t-g.earliest_timestep()][id] = active;
      }
    }

    return ret;
  }


  //
  // write_lgf()
  //
  namespace {
    struct TraxelToStrConverter {
      std::string operator()(const Traxel& t) {
	stringstream ss;
	boost::archive::text_oarchive oa(ss);
	oa & t;
	return ss.str();
      }
    };
  }

  void write_lgf( const HypothesesGraph& g, std::ostream& os, bool with_n_traxel ) {
    lemon::DigraphWriter<HypothesesGraph> writer( g, os );
    writer.
      nodeMap("timestep", g.get(node_timestep())).
      arcMap("from_timestep", g.get(arc_from_timestep())).
      arcMap("to_timestep", g.get(arc_to_timestep()));
    if(with_n_traxel) {
      writer.nodeMap("traxel", g.get(node_traxel()), TraxelToStrConverter());
    }
    writer.run();
  }



  //
  // read_lgf()
  //
  namespace {
    struct StrToTraxelConverter {
      Traxel operator()(const std::string& s) {
	stringstream ss(s);
	boost::archive::text_iarchive ia(ss);
	Traxel t;
	ia & t;
	return t;
      }
    };
  }

void read_lgf( HypothesesGraph& g, std::istream& is, bool with_n_traxel ) {
    lemon::DigraphReader<HypothesesGraph> reader( g, is );
    reader.
      nodeMap("timestep", g.get(node_timestep())).
      arcMap("from_timestep", g.get(arc_from_timestep())).
      arcMap("to_timestep", g.get(arc_to_timestep()));
    if( with_n_traxel ) {
      g.add(node_traxel());
      reader.nodeMap("traxel", g.get(node_traxel()), StrToTraxelConverter());
    }
    reader.run();
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


namespace {
  double getDivisionProbability(const Traxel& tr) {
  	FeatureMap::const_iterator it = tr.features.find("divProb");
  	if (it == tr.features.end()) {
  		throw runtime_error("getDivisionProbability(): divProb feature not in traxel");
  	}
  	return it->second[0];
  }
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

	HypothesesGraph* SingleTimestepTraxel_HypothesesBuilder::add_edges(
			HypothesesGraph* graph) const {
		LOG(logDEBUG) << "SingleTimestepTraxel_HypothesesBuilder::add_edges(): entered";
		typedef HypothesesGraph::node_timestep_map::Value timestep_t;
		const set<timestep_t>& timesteps = graph->timesteps();
		// iterate over all timesteps except the last
		for (set<timestep_t>::const_iterator t = timesteps.begin();
				t != (--timesteps.end()); ++t) {
			add_edges_at(graph, *t);
		}

		// if the forward_backward option is enabled, go again through the graph, this time
		// reversely and add the nearest neighbors if not already present
		if (options_.forward_backward) {
			// reversely iterate over all timesteps except the first
			for (set<timestep_t>::const_reverse_iterator t = timesteps.rbegin();
					t != (--timesteps.rend()); ++t) {
				add_edges_at(graph, *t, true);
			}
		}

		return graph;
	}

  HypothesesGraph* SingleTimestepTraxel_HypothesesBuilder::add_edges_at(HypothesesGraph* graph,
		  int timestep, bool reverse) const {
	const HypothesesGraph::node_timestep_map& timemap = graph->get(
			node_timestep());
	typedef property_map<node_traxel, HypothesesGraph::base_graph>::type traxelmap_t;
	const traxelmap_t& traxelmap = graph->get(node_traxel());
	const TraxelStoreByTimeid& traxels_by_timeid = ts_->get<by_timeid>();
	const TraxelStoreByTimestep& traxels_by_timestep = ts_->get<by_timestep>();

	int to_timestep = timestep + 1;
	if (reverse) {
		// iterating through the graph backward in time
		to_timestep = timestep - 1;
	}

	//// find k nearest neighbors in next timestep
	// init nearest neighbor search
	pair<TraxelStoreByTimestep::const_iterator,
			TraxelStoreByTimestep::const_iterator> traxels_at =
			traxels_by_timestep.equal_range(to_timestep);

	for (TraxelStoreByTimestep::const_iterator it = traxels_at.first;
			it != traxels_at.second; ++it) {
		assert(it->Timestep == to_timestep);
	}

	NearestNeighborSearch nns(traxels_at.first, traxels_at.second);


	// establish transition edges between a current node and appropriate nodes in next timestep
	for (HypothesesGraph::node_timestep_map::ItemIt curr_node(timemap,
			timestep); curr_node != lemon::INVALID; ++curr_node) {
		assert(timemap[curr_node] == timestep);
		assert(traxelmap[curr_node].Timestep == timestep);

		// if we want to consider divisions already in the Hypotheses graph
		// make sure that each potentially dividing cell has 2 nearest neighbors
		// (but only if we go through the graph forward in time)
		unsigned int max_nn = options_.max_nearest_neighbors;
		if (options_.consider_divisions && !reverse && max_nn < 2) {
			double div_prob = getDivisionProbability(traxelmap[curr_node]);
			if (div_prob > options_.division_threshold) {
				max_nn = 2;
			}
		}

		// search
		map<unsigned int, double> nearest_neighbors = nns.knn_in_range(
				traxelmap[curr_node], options_.distance_threshold,
				max_nn);

		//// connect current node with k nearest neighbor nodes
		for (map<unsigned int, double>::const_iterator neighbor =
				nearest_neighbors.begin(); neighbor != nearest_neighbors.end();
				++neighbor) {
			// connect with one of the neighbor nodes
			TraxelStoreByTimeid::iterator neighbor_traxel =
					traxels_by_timeid.find(
							boost::make_tuple(to_timestep, neighbor->first));
			assert(neighbor_traxel->Timestep == to_timestep);
			assert(neighbor_traxel->Timestep != traxelmap[curr_node].Timestep);
			traxelmap_t::ItemIt neighbor_node(traxelmap, *neighbor_traxel);
			assert(curr_node != neighbor_node);
			if (!reverse) {
				// if we go through the graph forward in time, add an arc from curr_node to neighbor_node
				graph->addArc(curr_node, neighbor_node);
			} else {
				// if we go through the graph backward in time, add an arc from neighbor_node to curr_node
				// if not already present
				if (lemon::findArc(*graph,neighbor_node,curr_node) == lemon::INVALID) {
					graph->addArc(neighbor_node, curr_node);
				}
			}
		}
	}

	return graph;
  }

} /* namespace pgmlink */
