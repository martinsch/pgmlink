#include "pgmlink/higher_order_features.h"

namespace pgmlink {

////
//// class track
////
Track::Track() {
  parent_id_ = 0;
  child_ids_.resize(2,0);
}

Track::Track(
  const size_t id,
  const size_t time_start,
  const Traxelvector& traxels,
  const size_t parent_id,
  const size_t left_child_id,
  const size_t right_child_id
) : id_(id), time_start_(time_start), traxels_(traxels), parent_id_(parent_id) {
  Track::set_child_ids(left_child_id, right_child_id);
}

void Track::set_id(const size_t id) {
  id_ = id;
}

size_t Track::get_id() const {
  return Track::id_;
}

void Track::set_time_start(const size_t time_start) {
  time_start_ = time_start;
}

size_t Track::get_time_start() const {
  return Track::time_start_;
}

void Track::set_parent_id(const size_t parent_id) {
  parent_id_ = parent_id;
}

size_t Track::get_parent_id() const {
  return Track::parent_id_;
}

void Track::set_child_ids(const size_t left, const size_t right) {
  child_ids_.resize(2);
  child_ids_[0] = left;
  child_ids_[1] = right;
}

const std::vector<size_t>& Track::get_child_ids() const {
  return Track::child_ids_;
}

size_t Track::get_length() const {
  return Track::traxels_.size();
}

////
//// Some auxiliary for building a tracking out of the hypotheses graph
////
typedef typename
  property_map<node_timestep, HypothesesGraph::base_graph>::type
  node_timestep_type;
typedef typename
  property_map<node_active_count, HypothesesGraph::base_graph>::type
  node_active_map_type;
typedef typename
  property_map<arc_active_count, HypothesesGraph::base_graph>::type
  arc_active_map_type;
typedef typename
  property_map<division_active_count, HypothesesGraph::base_graph>::type
  div_active_map_type;
typedef typename
  property_map<node_traxel, HypothesesGraph::base_graph>::type
  node_traxel_type;
typedef typename
  property_map<node_tracklet, HypothesesGraph::base_graph>::type
  node_tracklet_type;
typedef typename HypothesesGraph::NodeIt NodeIt;
typedef typename HypothesesGraph::InArcIt InArcIt;
typedef typename HypothesesGraph::OutArcIt OutArcIt;

Track track_from_start_node(
  const HypothesesGraph& graph,
  const HypothesesGraph::Node& node,
  const size_t index
) {
  bool tracklet_graph = graph.has_property(node_tracklet());
  node_timestep_type& timestep_map = graph.get(node_timestep());
  arc_active_map_type& arc_active_map = graph.get(arc_active_count());

  Traxelvector traxels;
  bool terminate = false;
  HypothesesGraph::Node n = node;
  while(not terminate) {
    // append the current traxel/tracklet to the end of the traxel vector
    if(tracklet_graph) {
      node_tracklet_type& tracklet_map = graph.get(node_tracklet());
      traxels.insert(
        traxels.end(),
        tracklet_map[n].begin(),
        tracklet_map[n].end()
      );
    } else {
      node_traxel_type& traxel_map = graph.get(node_traxel());
      traxels.push_back(traxel_map[n]);
    } // end if(tracklet_graph)

    std::vector<HypothesesGraph::Node> to_nodes;
    // get all outgoing arcs
    for(OutArcIt oa_it(graph, n); oa_it != lemon::INVALID; ++oa_it) {
      if(arc_active_map[oa_it][index]){
        to_nodes.push_back(graph.target(oa_it));
      }
    }
    // check if we have to terminate the track
    if(to_nodes.size() != 1) {
      terminate = true;
    } else {
      n = to_nodes.back();
    }// end if(to_nodes.size() !=1)
  } // end while(not terminate)
  return Track(0, timestep_map[node], traxels);
}

////
//// Class Tracking
////
Tracking::Tracking() {
}

Tracking::Tracking(const HypothesesGraph& graph, const size_t index)
: index_(index) {
  // check the if the properties are available
  assert(graph.has_property(node_active_count()));
  assert(graph.has_property(arc_active_count()));
  if (not graph.has_property(node_tracklet())) {
    assert(graph.has_property(node_traxel));
  }
  // get the property maps
  node_active_map_type& node_active_map = graph.get(node_active_count());
  arc_active_map_type& arc_active_map = graph.get(arc_active_count());
  bool has_div_active_map = graph.has_property(division_active_count());

  // Check for every node if we have to start a new track. Condition is:
  // node is active
  // AND
  // node has no incoming arc OR exactly one incoming arc from cell division
  for(NodeIt n_it(graph); n_it != lemon::INVALID; ++n_it) {
    if(node_active_map[n_it][index_]) {
      std::vector<HypothesesGraph::Node> from_nodes;
      // check for active incoming arcs
      for(InArcIt ia_it(graph, n_it); ia_it != lemon::INVALID; ++ia_it) {
        if(arc_active_map[ia_it][index_]){
          from_nodes.push_back(graph.source(ia_it));
        }
      }
      size_t num_in = from_nodes.size();
      bool parent_division_active = false;
      if (num_in==1) {
        HypothesesGraph::Node parent = from_nodes.front();
        if (has_div_active_map) {
          // if the graph has the division active property read the
          // parent_division_active property from this map
          parent_division_active = graph.get(
            division_active_count()
          )[parent][index_];
        } else {
          // otherwise iterate over all outgoing arcs of the parent to check if
          // there is a cell division
          size_t active_out_arcs = 0;
          for(OutArcIt oa_it(graph, parent); oa_it != lemon::INVALID; ++oa_it) {
            if(arc_active_map[oa_it][index_]) active_out_arcs++;
          }
          if (active_out_arcs==2) parent_division_active=true;
        }
      } // end if (num_in==1)
      if ((num_in==0) or (parent_division_active)) {
        // push back the track that starts from this node
        tracks_.push_back(track_from_start_node(graph, n_it, index));
        tracks_.back().set_id(tracks_.size());
      }
    }
  }
}

////
//// class TrackValue
////
template<typename TrackFeatureExtractor_T, typename FeatureAggregator_T>
feature_type TrackValue<
  TrackFeatureExtractor_T,
  FeatureAggregator_T
>::operator()(const Track& track) const {
  return featureaggregator_(trackfeatureextractor_(track));
}

} // namespace pgmlink
