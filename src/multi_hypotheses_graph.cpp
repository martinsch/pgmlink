

// pgmlink
#include <pgmlink/multi_hypotheses_graph.h>


namespace pgmlink {
  

  ////
  //// EventNode
  ////
  EventNode::EventNode(EventType type, unsigned from, unsigned to) :
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
  }
  

  ////
  //// MultiHypothesesGraph
  ////
  unsigned MultiHypothesesGraph::maximum_timestep() {
    return maximum_timestep_;
  }


  ////
  //// MultiHypothesesGraphBuilder
  ////
  boost::shared_ptr<MultiHypothesesGraph>
  MultiHypothesesGraphBuilder::build(RegionGraphVector) {
    // for adjacent timesteps do:
    // find connected components in range (kNN)
    // connect regions from those ccs using
    // appropriate event nodes and connection arcs
    // create conflict arcs for event nodes
  }
  

    
}

