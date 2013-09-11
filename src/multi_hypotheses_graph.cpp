// stl
#include <algorithm>
#include <set>
#include <vector>

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


unsigned MultiHypothesesGraph::maximum_timestep() {
  return maximum_timestep_;
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
  RegionGraphVector::iterator time_iterator = graphs->begin();
  RegionGraphVector::iterator time_plus_one_iterator = ++graphs->begin();
  TraxelVectorPtr traxels_at_t, traxels_at_t_plus_one;
  traxels_at_t_plus_one = extract_traxels(*time_iterator, 0u);
  add_nodes(graphs, graph);
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

void MultiHypothesesGraphBuilder::add_nodes(RegionGraphVectorPtr source_graphs,
                                            MultiHypothesesGraphPtr dest_graph) {
  RegionGraphVector::iterator time_iterator = source_graphs->begin();
  for (unsigned timestep = 0; time_iterator != source_graphs->end(); ++time_iterator, ++timestep) {
    add_nodes_at(*time_iterator, dest_graph, timestep);
  }
}


void MultiHypothesesGraphBuilder::add_nodes_at(RegionGraphPtr source_graph,
                                               MultiHypothesesGraphPtr dest_graph,
                                               unsigned timestep) {
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
  RegionGraph::ConnectedComponentMap& component_map = source_graph->get(node_connected_component());
  RegionGraph::LabelMap& label_map = source_graph->get(node_label());
  RegionGraph::TraxelMap& traxel_map = source_graph->get(node_traxel());
  RegionGraph::ConflictMap& conflict_map = source_graph->get(node_conflicts());
  RegionGraph::LevelMap& level_map = source_graph->get(node_level());
  MultiHypothesesGraph::ContainedRegionsMap& contained_regions_map = dest_graph->get(node_regions_in_component());

  const MultiHypothesesGraph::Node& node = dest_graph->add_node(timestep);
  /* reference_map_[source_node] = node;
  cross_reference_map_[node] = source_node; */
  
  
  Traxel trax = traxel_map[source_node];
  const std::set<RegionGraph::Node>& cc_conflicts = conflict_map[source_node];
  boost::shared_ptr<std::vector<label_type> > conflict_labels =
      source_graph->convert_nodes_to_property<node_label>(cc_conflicts.begin(), cc_conflicts.end());
  trax.features["conflicts"].assign(conflict_labels->begin(), conflict_labels->end());
  trax.features["root"].assign(1, 1.);
  trax.features["level"].assign(1, 0.);
  std::vector<Traxel>& contained_regions_traxel = contained_regions_map.get_value(node);
  contained_regions_traxel.push_back(trax);
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

  for (std::vector<Traxel>::iterator trax_it = contained_regions_traxel.begin()+1;
       trax_it != contained_regions_traxel.end();
       ++trax_it) {
    float& level = trax_it->features["level"][0];
    level = max_level - level + 1.; 
  }

  // add_conflict_graph_to_component(); to be done!
  // need to add a list of traxels (property map!) that contains information about:
  // -conflicts for each region
  // -level of each region
  // -features (already present in traxels)
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

