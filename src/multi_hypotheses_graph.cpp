// pgmlink
#include <pgmlink/multi_hypotheses_graph.h>
#include <pgmlink/traxels.h>
#include <pgmlink/nearest_neighbors.h>


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
    // for now: nothing to be done here
    // will change
  }


  unsigned MultiHypothesesGraph::maximum_timestep() {
    return maximum_timestep_;
  }


  ////
  //// MultiHypothesesGraphBuilder
  ////
  boost::shared_ptr<MultiHypothesesGraph>
  MultiHypothesesGraphBuilder::build(RegionGraphVectorPtr graphs) {
    boost::shared_ptr<MultiHypothesesGraph> graph(new MultiHypothesesGraph);
    // for adjacent timesteps do:
    RegionGraphVector::iterator time_iterator = graphs->begin();
    RegionGraphVector::iterator time_plus_one_iterator = ++graphs->begin();
    TraxelVectorPtr traxels_at_t, traxels_at_t_plus_one;
    // traxels_at_t_plus_one = extract_traxels(*time_iterator);
    for (; time_plus_one_iterator != graphs->end();
         ++time_iterator, ++time_plus_one_iterator) {
      // find connected components in range (kNN)
      traxels_at_t = traxels_at_t_plus_one;
      // traxels_at_t_plus_one = extract_traxels(*time_plus_one_iterator);
      NearestNeighborSearch nearest_neighbor_search(traxels_at_t->begin(),
                                                    traxels_at_t->end()
                                                    );
      for (TraxelVector::iterator traxel_it = traxels_at_t_plus_one->begin();
           traxel_it != traxels_at_t_plus_one->end();
           ++traxel_it) {
        map<unsigned, double> nearest_neighbors =
          nearest_neighbor_search.knn_in_range(*traxel_it,
                                               options_.distance_threshold,
                                               options_.max_nearest_neighbors,
                                               options_.forward_backward
                                               );
      }
      // connect regions from those ccs using
      // appropriate event nodes and connection arcs
      // create conflict arcs for event nodes
    }

    return graph;
  }
  

    
}

