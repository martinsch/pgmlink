// stl
#include <algorithm>

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
  for (; time_plus_one_iterator != graphs->end();
       ++time_iterator, ++time_plus_one_iterator) {
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
      
      create_events_for_component(*traxel_it,
                                  reduce_to_nearest_neighbors(traxels_at_t_plus_one,
                                                              nearest_neighbors),
                                  *time_iterator,
                                  graph);
    }
    // connect regions from those ccs using
    // appropriate event nodes and connection arcs
    // create conflict arcs for event nodes
  }
  
  return graph;
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
    nearest_traxels->push_back(
        *std::find(traxels->begin(), traxels->end(), Traxel(neighbor->first)));
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

