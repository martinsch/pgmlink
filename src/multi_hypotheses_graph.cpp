// stl
#include <algorithm>
#include <set>
#include <vector>
#include <cassert>

// boost
#include <boost/assert.hpp>

// omp
#include <omp.h>

// pgmlink
#include <pgmlink/multi_hypotheses_graph.h>
#include <pgmlink/traxels.h>
#include <pgmlink/nearest_neighbors.h>
#include <pgmlink/hypotheses.h> // for tags
#include <pgmlink/region_graph.h> // more tags and RegionGraph
#include <pgmlink/classifier_auxiliary.h>


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


void MultiHypothesesGraph::add_classifier_features(ClassifierStrategy* move,
                                                   ClassifierStrategy* division,
                                                   ClassifierStrategy* count,
                                                   ClassifierStrategy* detection) {
  LOG(logDEBUG) << "MultiHypothesesGraph::add_classifier_features: entered";
  add(node_move_features());
  add(node_division_features());
  add(node_count_features());

  DivisionFeatureMap& division_map = get(node_division_features());
  CountFeatureMap& count_map = get(node_count_features());
  MoveFeatureMap& move_map = get(node_move_features());
  ContainedRegionsMap& regions = get(node_regions_in_component());
  
  std::vector<NodeT> nodes = this->nodes;

  // allocate space for results before going parallel
  for (size_t i = 0; i < nodes.size(); ++i) { // need this for omp
    	Node n = this->nodeFromId(i);
    	assert(this->valid(n));
      //  for (NodeIt n(*this); n != lemon::INVALID; ++n) {
      LOG(logDEBUG1) << "MultiHypothesesGraph::add_classifier_features: classifying region "
                     << get(node_traxel())[n].Id << " at timestep "
                     << get(node_traxel())[n].Timestep;
      std::vector<Traxel>& sources = regions.get_value(n);
      std::vector<Traxel> targets;
      for (OutArcIt a(*this, n); a != lemon::INVALID; ++a) {
        const std::vector<Traxel>& target = regions[this->target(a)];
        LOG(logDEBUG4) << "MultiHypothesesGraph::add_classifier_features: added "
                       << target.size() << " regions to target";
        targets.insert(targets.end(), target.begin(), target.end());
      }
      LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: initializing moves";
      move->classify(sources, targets, move_map.get_value(n), false); // with_predict = false
      LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: initializing divisions";
      division->classify(sources, targets, division_map.get_value(n), false); // with_predict = false
      LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: initializing detections";
      detection->classify(sources, false); // with_predict = false
    }


  // doing the actual predictions in parallel
  # pragma omp parallel for
  for (size_t i = 0; i < nodes.size(); ++i) { // need this for omp
  	Node n = this->nodeFromId(i);
  	assert(this->valid(n));
    //  for (NodeIt n(*this); n != lemon::INVALID; ++n) {
    LOG(logDEBUG1) << "MultiHypothesesGraph::add_classifier_features: classifying region "
                   << get(node_traxel())[n].Id << " at timestep "
                   << get(node_traxel())[n].Timestep;
    std::vector<Traxel>& sources = regions.get_value(n);
    std::vector<Traxel> targets; 
    for (OutArcIt a(*this, n); a != lemon::INVALID; ++a) {
      const std::vector<Traxel>& target = regions[this->target(a)];
      LOG(logDEBUG4) << "MultiHypothesesGraph::add_classifier_features: added "
                     << target.size() << " regions to target";
      targets.insert(targets.end(), target.begin(), target.end());
    }
    LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: classifying moves";
    move->classify(sources, targets, move_map.get_value(n), true); // with_predict = true
    LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: classifying divisions";
    division->classify(sources, targets, division_map.get_value(n), true); // with_predict = true
    LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: classifying detections";
    detection->classify(sources, true); // with_predict = true
    LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: classifying count";
    count->classify(sources);
  }
  
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
  graph->add(node_conflict_sets());
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
  LOG(logDEBUG) << "MultiHypothesesGraphBuilder::add_nodes() -- entered";
  MultiHypothesesGraph::ContainedRegionsMap& regions = dest_graph->get(node_regions_in_component());
  MultiHypothesesGraph::TraxelMap& traxel_map = dest_graph->get(node_traxel());
  MultiHypothesesGraph::ConflictSetMap& conflict_sets = dest_graph->get(node_conflict_sets());
  for (TimestepRegionMap::const_iterator timestep = ts.map.begin();
       timestep != ts.map.end();
       ++timestep) {
    for (std::map<unsigned, std::vector<Traxel> >::const_iterator component = timestep->second.begin();
         component != timestep->second.end();
         ++component) {
      const MultiHypothesesGraph::Node& node = dest_graph->add_node(timestep->first);
      std::vector<Traxel>& traxels = regions.get_value(node);


      LOG(logDEBUG4) << "MultiHypothesesGraphBuilder::add_nodes() -- rearranging traxel vector to make "
                     << "connected component traxel is at the beginning of the vector";
      std::vector<Traxel>::const_iterator base =
          std::find(component->second.begin(), component->second.end(), Traxel(component->first, timestep->first));
      assert(base != component->second.end());
      traxels.push_back(*base);
      for (std::vector<Traxel>::const_iterator r = component->second.begin();
           r != component->second.end();
           ++r) {
        // only add traxel if it is not the connected component which is already present in the vector
        if (r != base) {
          traxels.push_back(*r);
        }
      }
      
      assert(traxels.size() > 0);
      assert(traxels.size() == component->second.size());
      traxel_map.set(node, traxels[0]);
      assert(ts.conflicts_by_timestep.find(timestep->first) != ts.conflicts_by_timestep.end());
      const ConflictSetMap& conflicts_source = ts.conflicts_by_timestep.find(timestep->first)->second;
      ConflictSetMap::const_iterator conflicts_source_at = conflicts_source.find(component->first);
      assert(conflicts_source_at != conflicts_source.end());
      
      MultiHypothesesGraph::ConflictSetMap::Value& conflicts_target_at = conflict_sets.get_value(node);
      assert(conflicts_target_at.size() == 0);
      conflicts_target_at.insert(conflicts_target_at.end(), conflicts_source_at->second.begin(), conflicts_source_at->second.end());
      LOG(logDEBUG4) << "MultiHypothesesGraphBuilder::add_nodes -- "
                     << "added new component " << traxels[0].Id
                     << " at time " << traxels[0].Timestep;
      LOG(logDEBUG4) << "MultiHypothesesGraphBuilder::add_nodes -- "
                     <<  "new component has " << conflict_sets[node].size()
                     << " conflict sets (" << traxels[0] << " had " << conflicts_source_at->second.size() << " conflicts)";
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
  LOG(logDEBUG3) << "MultiHypothesesGraphBuilder::add_edges_at -- initializing knn";
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
    } else if (lemon::findArc(*graph, neighbor_node, node) == lemon::INVALID) {
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


////
//// MultiHypothesesTraxelStore
////
void MultiHypothesesTraxelStore::add(const Traxel& trax, unsigned component_id) {
  map[trax.Timestep][component_id].push_back(trax);
}


void MultiHypothesesTraxelStore::add_conflict_map(int timestep, const ConflictSetMap& conflicts) {
  conflicts_by_timestep[timestep] = conflicts;
}


void MultiHypothesesTraxelStore::add_signed_conflict_map(int timestep, const SignedConflictSetMap& conflicts) {
  LOG(logDEBUG2) << "MultiHypothesesTraxelStore::add_signed_conflict_map()";
  for (SignedConflictSetMap::const_iterator region = conflicts.begin(); region != conflicts.end(); ++region) {
    for (std::vector<std::vector<int> >::const_iterator conflict = region->second.begin();
         conflict != region->second.end();
         ++conflict) {
      LOG(logDEBUG4) << "MultiHypothesesTraxelStore::add_signed_conflict_map() -- adding set containing "
                     << conflict->size() << " conflicts @ " << timestep << "," << region->first;
      conflicts_by_timestep[timestep][region->first].push_back(std::vector<unsigned>());
      std::vector<unsigned>& new_conflict = *conflicts_by_timestep[timestep][region->first].rbegin();
      new_conflict.insert(new_conflict.end(), conflict->begin(), conflict->end());
    }
  }
}


void MultiHypothesesTraxelStore::start_component(const Traxel& trax) {
  LOG(logDEBUG4) << "MultiHypothesesTraxelStore::start_component: " << trax;
  map[trax.Timestep][trax.Id] = std::vector<Traxel>();
  assert(map[trax.Timestep][trax.Id].size() == 0);
  map[trax.Timestep][trax.Id].push_back(trax);
  LOG(logDEBUG4) << "MultiHypothesesTraxelStore::start_component: "
                 << "new trax has com?"
                 << map[trax.Timestep][trax.Id].rbegin()->features.count("com");
}


void MultiHypothesesTraxelStore::make_sane() {
  for (TimestepRegionMap::const_iterator t = map.begin(); t != map.end(); ++t) {
    if (conflicts_by_timestep.find(t->first) == conflicts_by_timestep.end()) {
      conflicts_by_timestep[t->first] = ConflictSetMap();
    }
    for (std::map<unsigned, std::vector<Traxel> >::const_iterator r = t->second.begin(); r != t->second.end(); ++r) {
      if (conflicts_by_timestep[t->first].find(r->first) == conflicts_by_timestep[t->first].end()) {
        conflicts_by_timestep[t->first][r->first] = std::vector<std::vector<unsigned> >();
      }
    }
  }
}


std::string MultiHypothesesTraxelStore::print() {
  std::stringstream ss;
  for (TimestepRegionMap::const_iterator it = map.begin(); it != map.end(); ++it) {
    ss << "t=" << it->first << '\n';
    for (std::map<unsigned, std::vector<Traxel> >::const_iterator it2 = it->second.begin();
         it2 != it->second.end();
         ++it2) {
      ss << ".. connected component=" << it2->first << '\n';
      for (std::vector<Traxel>::const_iterator it3 = it2->second.begin();
           it3 != it2->second.end();
           ++it3) {
        ss << ".... region=" << *it3 << '\n';
      }
    }
    ss << "\n";
  }
  return ss.str();
}

}

