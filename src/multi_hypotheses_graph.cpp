// stl
#include <algorithm>
#include <set>
#include <vector>

// boost
#include <boost/assert.hpp>

// pgmlink
#include <pgmlink/multi_hypotheses_graph.h>
#include <pgmlink/traxels.h>
#include <pgmlink/nearest_neighbors.h>
#include <pgmlink/hypotheses.h> // for tags
#include <pgmlink/region_graph.h> // more tags and RegionGraph


namespace pgmlink {
  

////
//// EventNode
////
/*EventNode::EventNode(EventType type, unsigned from, unsigned to) :
  type_(type), from_(from), to_(to), distances_(feature_array(0)) {
  // nothing to be done here
  }
  
  
  EventNode::EventType EventNode::get_type() {
  return type_;
  }

  
  unsigned EventNode::from() {
  return from_;
  }


  unsigned EventNode::to() {
  return to_;
  }*/
  

////
//// MultiHypothesesGraph
////
MultiHypothesesGraph::MultiHypothesesGraph() {
  add(node_traxel());
  add(node_regions_in_component());
}


////
//// MultiHypothesesGraphBuilder
////
MultiHypothesesGraphBuilder::MultiHypothesesGraphBuilder(const Options& options) :
    options_(options),
    reference_map_(),
    cross_reference_map_() {
  // nothing to be done here
}


MultiHypothesesGraphPtr
MultiHypothesesGraphBuilder::build(RegionGraphVectorPtr graphs) {
  MultiHypothesesGraphPtr graph(new MultiHypothesesGraph);
  // for adjacent timesteps do:
  // RegionGraphVector::iterator time_iterator = graphs->begin();
  // RegionGraphVector::iterator time_plus_one_iterator = ++graphs->begin();
  // TraxelVectorPtr traxels_at_t, traxels_at_t_plus_one;
  // traxels_at_t_plus_one = extract_traxels(*time_iterator, 0u);
  add_nodes(graphs, graph);
  add_edges(graph);
  /* for (int timestep = 0;
       time_plus_one_iterator != graphs->end();
       ++time_iterator, ++time_plus_one_iterator, ++timestep) {
    add_nodes(*time_plus_one_iterator, graph, 0u, timestep);
    // find connected components in range (kNN)
    traxels_at_t = traxels_at_t_plus_one;
    traxels_at_t_plus_one = extract_traxels(*time_plus_one_iterator, 0u);
    NearestNeighborSearch nearest_neighbor_search(traxels_at_t_plus_one->begin(),
                                                  traxels_at_t_plus_one->end()
                                                  );
    for (TraxelVector::iterator traxel_it = traxels_at_t->begin();
         traxel_it != traxels_at_t->end();
         ++traxel_it) {
      std::map<unsigned, double> nearest_neighbors =
          nearest_neighbor_search.knn_in_range(*traxel_it,
                                               options_.distance_threshold,
                                               options_.max_nearest_neighbors,
                                               options_.forward_backward
                                               );
      TraxelVectorPtr nearest_neighbors_traxels = reduce_to_nearest_neighbors(traxels_at_t_plus_one,
                                                                              nearest_neighbors);

      // Create arcs from connected components in the nearest neighbor range.
      // Store information about the connected components in property maps
      
      
      create_events_for_component(*traxel_it,
                                  reduce_to_nearest_neighbors(traxels_at_t_plus_one,
                                                              nearest_neighbors),
                                  *time_iterator,
                                  graph);
    }
    // connect regions from those ccs using
    // appropriate event nodes and connection arcs
    // create conflict arcs for event nodes
  } */
  
  return graph;
}


MultiHypothesesGraphPtr
MultiHypothesesGraphBuilder::build(const MultiHypothesesTraxelStore& ts) {
  MultiHypothesesGraphPtr graph(new MultiHypothesesGraph);
  add_nodes(ts, graph);
  add_edges(graph);

  return graph;
}

void MultiHypothesesGraphBuilder::add_nodes(RegionGraphVectorPtr source_graphs,
                                            MultiHypothesesGraphPtr dest_graph) {
  LOG(logDEBUG) << "MultiHypothesesGraphBuilder::add_nodes -- entered";
  RegionGraphVector::iterator time_iterator = source_graphs->begin();
  for (unsigned timestep = 0; time_iterator != source_graphs->end(); ++time_iterator, ++timestep) {
    add_nodes_at(*time_iterator, dest_graph, timestep);
  }
}

void MultiHypothesesGraphBuilder::add_nodes(const MultiHypothesesTraxelStore& ts,
                                            MultiHypothesesGraphPtr dest_graph) {
  LOG(logDEBUG) << "MultiHypothesesGraphBuilder::add_nodes -- entered";
  MultiHypothesesGraph::ContainedRegionsMap& regions = dest_graph->get(node_regions_in_component());
  for (TimestepRegionMap::const_iterator timestep = ts.map.begin();
       timestep != ts.map.end();
       ++timestep) {
    for (std::map<unsigned, std::pair<Traxel, std::vector<Traxel> > >::const_iterator component
             = timestep->second.begin();
         component != timestep->second.end();
         ++component) {
      const MultiHypothesesGraph::Node& node = dest_graph->add_node(timestep->first);
      std::vector<Traxel>& traxels = regions.get_value(node);
      traxels.push_back(component->second.first);
      traxels.insert(traxels.end(), component->second.second.begin(), component->second.second.end());
    }
  }
}


void MultiHypothesesGraphBuilder::add_nodes_at(RegionGraphPtr source_graph,
                                               MultiHypothesesGraphPtr dest_graph,
                                               unsigned timestep) {
  LOG(logDEBUG) << "MultiHypothesesGraphBuilder::add_nodes_at -- entered";
  RegionGraph::ConnectedComponentMap& component_map = source_graph->get(node_connected_component());
  for (RegionGraph::ConnectedComponentMap::ItemIt source_it(component_map, 0u);
       source_it != lemon::INVALID;
       ++source_it) {
    add_node(source_graph, dest_graph, source_it, timestep);
  }
}



void MultiHypothesesGraphBuilder::add_node(RegionGraphPtr source_graph,
                                           MultiHypothesesGraphPtr dest_graph,
                                           const RegionGraph::Node& source_node,
                                           int timestep) {
  LOG(logDEBUG) << "MultiHypothesesGraphBuilder::add_node -- entered";
  RegionGraph::ConnectedComponentMap& component_map = source_graph->get(node_connected_component());
  // RegionGraph::LabelMap& label_map = source_graph->get(node_label());
  RegionGraph::TraxelMap& traxel_map = source_graph->get(node_traxel());
  RegionGraph::ConflictMap& conflict_map = source_graph->get(node_conflicts());
  RegionGraph::LevelMap& level_map = source_graph->get(node_level());
  MultiHypothesesGraph::ContainedRegionsMap& contained_regions_map = dest_graph->get(node_regions_in_component());

  MultiHypothesesGraph::TraxelMap& multi_traxel_map = dest_graph->get(node_traxel());
  
  const MultiHypothesesGraph::Node& node = dest_graph->add_node(timestep);
  
  
  Traxel trax = traxel_map[source_node];
  const std::set<RegionGraph::Node>& cc_conflicts = conflict_map[source_node];
  boost::shared_ptr<std::vector<label_type> > conflict_labels =
      source_graph->convert_nodes_to_property<node_label>(cc_conflicts.begin(), cc_conflicts.end());
  trax.features["conflicts"].assign(conflict_labels->begin(), conflict_labels->end());
  trax.features["root"].assign(1, 1.);
  trax.features["level"].assign(1, 0.);
  std::vector<Traxel>& contained_regions_traxel = contained_regions_map.get_value(node);
  contained_regions_traxel.push_back(trax);
  multi_traxel_map.set(node, trax);
  float max_level = 0;
  for (RegionGraph::ConnectedComponentMap::ItemIt region(component_map, trax.Id);
       region != lemon::INVALID;
       ++region) {
    trax = traxel_map[region];
    const std::set<RegionGraph::Node>& conflicts = conflict_map[region];
    conflict_labels = source_graph->convert_nodes_to_property<node_label>(conflicts.begin(), conflicts.end());
    trax.features["conflicts"].assign(conflict_labels->begin(), conflict_labels->end());
    trax.features["root"].assign(1, 0.);
    float level = level_map[region];
    max_level = std::max(max_level, level);
    trax.features["level"].assign(1, level);
    contained_regions_traxel.push_back(trax);
  }

  // set levels properly
  for (std::vector<Traxel>::iterator trax_it = contained_regions_traxel.begin()+1;
       trax_it != contained_regions_traxel.end();
       ++trax_it) {
    float& level = trax_it->features["level"][0];
    level = max_level - level + 1.; 
  }
}





void MultiHypothesesGraphBuilder::add_edges(MultiHypothesesGraphPtr graph) {
  LOG(logDEBUG) << "MultiHypothesesGraphBuilder::add_edges -- entered";
  int max_timestep = graph->latest_timestep();
  int min_timestep = graph->earliest_timestep();
  for (int timestep = min_timestep;
       timestep < max_timestep;
       timestep += FORWARD) {
    add_edges_at(graph, timestep, FORWARD);
  }
  if (options_.forward_backward) {
    for (int timestep = max_timestep;
         timestep > min_timestep;
         timestep += BACKWARD) {
      add_edges_at(graph, timestep, BACKWARD);
    }
  }
}


void MultiHypothesesGraphBuilder::add_edges_at(MultiHypothesesGraphPtr graph,
                                               int timestep,
                                               DIRECTION direction) {
  int next_timestep = timestep + direction;
  LOG(logDEBUG) << "MultiHypothesesGraphBuilder::add_edges_at -- entered, timestep="
                << timestep << " next_timestep=" << next_timestep;
  MultiHypothesesGraph::node_timestep_map& timestep_map = graph->get(node_timestep());
  MultiHypothesesGraph::TraxelMap& traxel_map = graph->get(node_traxel());
  TraxelVectorPtr traxels_at_next = graph->get_properties_at<node_traxel>(next_timestep);
  NearestNeighborSearch nearest_neighbor_search(traxels_at_next->begin(),
                                                traxels_at_next->end());
  for (MultiHypothesesGraph::node_timestep_map::ItemIt iterator_at(timestep_map, timestep);
       iterator_at != lemon::INVALID;
       ++iterator_at) {
    std::map<unsigned, double> nearest_neighbors =
        nearest_neighbor_search.knn_in_range(traxel_map[iterator_at],
                                             options_.distance_threshold,
                                             options_.max_nearest_neighbors,
                                             direction == BACKWARD
                                             );
    add_edges_for_node(graph, iterator_at, nearest_neighbors, timestep, direction);
  }
}


void MultiHypothesesGraphBuilder::add_edges_for_node(MultiHypothesesGraphPtr graph,
                                                     const MultiHypothesesGraph::Node& node,
                                                     std::map<unsigned, double>& neighbors,
                                                     int timestep,
                                                     DIRECTION direction) {
  LOG(logDEBUG) << "MultiHypothesesGraphBuilder::add_edges_for_node -- entered";
  MultiHypothesesGraph::TraxelMap& traxel_map = graph->get(node_traxel());
  for (std::map<unsigned, double>::iterator neighbor = neighbors.begin();
       neighbor != neighbors.end();
       ++neighbor) {
    MultiHypothesesGraph::TraxelMap::ItemIt neighbor_node(traxel_map, Traxel(neighbor->first, timestep + direction));
    assert(neighbor_node != lemon::INVALID);
    assert(traxel_map[neighbor_node].Timestep == timestep + direction);
    assert(neighbor_node != node);
    if (direction == FORWARD) {
      graph->addArc(node, neighbor_node);
    } else {
      graph->addArc(neighbor_node, node);
    }
  }
}


TraxelVectorPtr MultiHypothesesGraphBuilder::extract_traxels(RegionGraphPtr graph,
                                                             unsigned cc_label) {
  TraxelVectorPtr traxels(new TraxelVector);
  RegionGraph::NodeVectorPtr nodes = graph->get_nodes_in_component(cc_label);
  RegionGraph::TraxelMap& traxel_map = graph->get(node_traxel());
  for (RegionGraph::NodeVector::iterator node = nodes->begin();
       node != nodes->end();
       ++node) {
    traxels->push_back(traxel_map[*node]);
  }
  return traxels;
}


TraxelVectorPtr MultiHypothesesGraphBuilder::
reduce_to_nearest_neighbors(TraxelVectorPtr traxels,
                             std::map<unsigned, double>& neighbors) {
  TraxelVectorPtr nearest_traxels(new TraxelVector);
  for (std::map<unsigned, double>::iterator neighbor = neighbors.begin();
       neighbor != neighbors.end();
       ++neighbor) {
    TraxelVector::iterator find_res =
        std::find(traxels->begin(), traxels->end(), Traxel(neighbor->first));
    if (find_res != traxels->end()) {
      nearest_traxels->push_back(*find_res);
    }
  }
  return nearest_traxels;
}
                                                                          


void MultiHypothesesGraphBuilder::create_events_for_component(const Traxel& trax,
                                                              TraxelVectorPtr nearest_neighbors,
                                                              RegionGraphPtr single_timestep_graph,
                                                              MultiHypothesesGraphPtr graph) {
  TraxelVector single_component_all_regions;
  TraxelVector regions_in_component =
      *extract_traxels(single_timestep_graph, trax.Id);
  single_component_all_regions.push_back(trax);
  single_component_all_regions.insert(single_component_all_regions.end(),
                                      regions_in_component.begin(),
                                      regions_in_component.end()
                                      );
  NearestNeighborSearch nearest_neighbor_search(nearest_neighbors->begin(),
                                                nearest_neighbors->end()
                                                );
  for (TraxelVector::iterator region = single_component_all_regions.begin();
       region != single_component_all_regions.end();
       ++region) {
    std::map<unsigned, double> nearest_neighbors_map =
        nearest_neighbor_search.knn_in_range(*region,
                                             options_.distance_threshold,
                                             options_.max_nearest_neighbors,
                                             options_.forward_backward
                                             );
    add_move_events(*region, nearest_neighbors_map, graph);
    add_division_events(*region, nearest_neighbors_map, graph);
    add_appearance_events(*region, graph);
    add_disappearance_events(*region, graph);
  }
}


void MultiHypothesesGraphBuilder::add_move_events(const Traxel& trax,
                                                  std::map<unsigned, double>& neighbors,
                                                  MultiHypothesesGraphPtr graph) {
  
}


void MultiHypothesesGraphBuilder::add_division_events(const Traxel& trax,
                                                      std::map<unsigned, double>& neighbors,
                                                      MultiHypothesesGraphPtr graph) {

}


void MultiHypothesesGraphBuilder::add_appearance_events(const Traxel& trax,
                                                        MultiHypothesesGraphPtr graph) {

}


void MultiHypothesesGraphBuilder::add_disappearance_events(const Traxel& trax,
                                                           MultiHypothesesGraphPtr graph) {

}
  
}

