#include "pgmlink/higher_order_features.h"
#include <lemon/adaptors.h> // for SubDigraph
#include <boost/pending/disjoint_sets.hpp>

namespace pgmlink {

////
//// Some useful typedefinitions
////
typedef typename
  property_map<node_active, HypothesesGraph::base_graph>::type
  node_active_map_type;
typedef typename
  property_map<arc_active, HypothesesGraph::base_graph>::type
  arc_active_map_type;
typedef typename
  property_map<node_timestep, HypothesesGraph::base_graph>::type
  node_timestep_type;
typedef typename
  property_map<node_active_count, HypothesesGraph::base_graph>::type
  nodes_active_map_type;
typedef typename
  property_map<arc_active_count, HypothesesGraph::base_graph>::type
  arcs_active_map_type;
typedef typename
  property_map<division_active_count, HypothesesGraph::base_graph>::type
  divs_active_count_map_type;
typedef typename
  property_map<node_traxel, HypothesesGraph::base_graph>::type
  node_traxel_type;
typedef typename
  property_map<node_tracklet, HypothesesGraph::base_graph>::type
  node_tracklet_type;
typedef typename HypothesesGraph::NodeIt NodeIt;
typedef typename HypothesesGraph::ArcIt ArcIt;
typedef typename HypothesesGraph::InArcIt InArcIt;
typedef typename HypothesesGraph::OutArcIt OutArcIt;

typedef typename node_active_map_type::TrueIt NodeActiveIt;
typedef typename arc_active_map_type::TrueIt ArcActiveIt;


////
//// function set_solution
////
void set_solution(HypothesesGraph& graph, const size_t solution_index) {
  // check if the graph has the necessary property maps
  if (not graph.has_property(node_active_count())) {
    throw std::runtime_error(
      "Graph doesn't have a \"node_active_count\" property map"
    );
  }
  if (not graph.has_property(arc_active_count())) {
    throw std::runtime_error(
      "Graph doesn't have an \"arc_active_count\" property map"
    );
  }

  // Get the property maps
  nodes_active_map_type& nodes_active_map = graph.get(node_active_count());
  arcs_active_map_type& arcs_active_map = graph.get(arc_active_count());

  // check if the solution_index is legal
  if (nodes_active_map.beginValue()->size() <= solution_index) {
    throw std::runtime_error("Index of solution out of range");
  }

  // create the a node_active and arc_active map
  if (not graph.has_property(node_active())) {
    graph.add(node_active());
  }
  if (not graph.has_property(arc_active())) {
    graph.add(arc_active());
  }

  // Get the property maps (caution: "node_active_map" not "nodes_active_map")
  node_active_map_type& node_active_map = graph.get(node_active());
  arc_active_map_type& arc_active_map = graph.get(arc_active());

  // Now we can start the with writing the solution with the index
  // solution_index into the node_active_map
  for (NodeIt n_it(graph); n_it != lemon::INVALID; ++n_it) {
    node_active_map[n_it] = nodes_active_map[n_it][solution_index];
  }
  for (ArcIt a_it(graph); a_it != lemon::INVALID; ++a_it) {
    arc_active_map[a_it] = arcs_active_map[a_it][solution_index];
  }
}

////
//// class track
////
Track::Track() {
  parent_ = NULL;
  children_.clear();
}

Track::Track(
  const size_t id,
  const size_t time_start,
  const Traxelvector& traxels,
  Track* const parent,
  Track* const left_child,
  Track* const right_child
) : id_(id), time_start_(time_start), traxels_(traxels), parent_(parent) {
  Track::set_child(left_child);
  Track::set_child(right_child);
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

void Track::set_parent(Track* const parent) {
  parent_ = parent;
}

Track* Track::get_parent() const {
  return Track::parent_;
}

void Track::set_child(Track* const child) {
  if (child != NULL) {
    LOG(logDEBUG4) << "Set child in track " << id_;
    children_.push_back(child);
  }
}

void Track::set_children(Track* const left_child, Track* const right_child) {
  children_.resize(2);
  children_[0] = left_child;
  children_[1] = right_child;
}

const std::vector<Track*>& Track::get_children() const {
  return Track::children_;
}

size_t Track::get_length() const {
  return Track::traxels_.size();
}

Track track_from_start_node(
  const HypothesesGraph& graph,
  const HypothesesGraph::Node& node,
  const size_t index,
  HypothesesGraph::Node& last_node
) {
  bool tracklet_graph = graph.has_property(node_tracklet());
  node_timestep_type& timestep_map = graph.get(node_timestep());
  arcs_active_map_type& arcs_active_map = graph.get(arc_active_count());

  Traxelvector traxels;
  bool terminate = false;
  last_node = node;
  while(not terminate) {
    // append the current traxel/tracklet to the end of the traxel vector
    if(tracklet_graph) {
      node_tracklet_type& tracklet_map = graph.get(node_tracklet());
      traxels.insert(
        traxels.end(),
        tracklet_map[last_node].begin(),
        tracklet_map[last_node].end()
      );
    } else {
      node_traxel_type& traxel_map = graph.get(node_traxel());
      traxels.push_back(traxel_map[last_node]);
    } // end if(tracklet_graph)

    std::vector<HypothesesGraph::Node> to_nodes;
    // get all outgoing arcs
    for(OutArcIt oa_it(graph, last_node); oa_it != lemon::INVALID; ++oa_it) {
      if(arcs_active_map[oa_it][index]){
        to_nodes.push_back(graph.target(oa_it));
      }
    }
    // check if we have to terminate the track
    if(to_nodes.size() != 1) {
      terminate = true;
    } else {
      last_node = to_nodes.back();
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
    assert(graph.has_property(node_traxel()));
  }
  // stores reference to the track with the last node of the track as the key
  std::map<int, size_t> node_in_track;
  // stores the parent node of a track with the track as the key
  std::map<size_t, int> has_parent_node;

  // get the property maps
  nodes_active_map_type& nodes_active_map = graph.get(node_active_count());
  arcs_active_map_type& arcs_active_map = graph.get(arc_active_count());
  bool has_div_active_map = graph.has_property(division_active_count());

  // Check for every node if we have to start a new track. Condition is:
  // node is active
  // AND
  // node has no incoming arc OR exactly one incoming arc from cell division
  for(NodeIt n_it(graph); n_it != lemon::INVALID; ++n_it) {
    if(nodes_active_map[n_it][index_]) {
      std::vector<HypothesesGraph::Node> from_nodes;
      // check for active incoming arcs
      for(InArcIt ia_it(graph, n_it); ia_it != lemon::INVALID; ++ia_it) {
        if(arcs_active_map[ia_it][index_]){
          from_nodes.push_back(graph.source(ia_it));
        }
      }
      size_t num_in = from_nodes.size();
      bool parent_division_active = false;
      // set the parent node to lemon::INVALID. If there is an active parent
      // node this variable will be updated.
      HypothesesGraph::Node parent = lemon::INVALID;
      if (num_in==1) {
        // update the parent node
        parent = from_nodes.front();
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
            if(arcs_active_map[oa_it][index_]) active_out_arcs++;
          }
          if (active_out_arcs==2) parent_division_active=true;
        }
      } // end if (num_in==1)
      if ((num_in==0) or (parent_division_active)) {
        // initialize the track id and the variable in which the last node will
        // be stored
        HypothesesGraph::Node last_node;
        size_t track_id = tracks_.size();
        // create the track and get the last node of it and set the id
        tracks_.push_back(track_from_start_node(graph, n_it, index, last_node));
        tracks_.back().set_id(track_id);
        // update the maps "node_in_track" and "has_parent_node"
        node_in_track[graph.id(last_node)] = track_id;
        has_parent_node[track_id] = graph.id(parent);
      }
    }
  }
  // iterate over all tracks to set the connections
  for (
    Trackvector::iterator t_it = tracks_.begin();
    t_it != tracks_.end();
    t_it++
  ) {
    // get the parent node
    int parent_node_id = has_parent_node[t_it->get_id()];
    if (parent_node_id != -1) {
      // get the map iterator to the parent track and check if it exists
      std::map<int, size_t>::iterator parent_track_it;
      parent_track_it = node_in_track.find(parent_node_id);
      if (parent_track_it == node_in_track.end()) {
        LOG(logDEBUG) << "In constructor of tracking: inconsistent or unsupported graph";
        LOG(logDEBUG) << "Node " << parent_node_id << " doesn't terminate track (although it should)";
      } else {
        // set the references
        Track* parent_track_ref = &(tracks_[parent_track_it->second]);
        t_it->set_parent(parent_track_ref);
        parent_track_ref->set_child(&(*t_it));
      }
    }
  }
}

/*=============================================================================
  pure virtual classes
=============================================================================*/
////
//// class SubsetFeatureExtractor
////
const FeatureMatrix& SubsetFeatureExtractor::extract_matrix(
  const ConstTraxelRefVector& /*traxelrefs*/
) {
  throw std::runtime_error(
    "Feature extractor " + name() + " doesn't provide a matrix valued result"
  );
  return *(new FeatureMatrix);
}
const FeatureVector& SubsetFeatureExtractor::extract_vector(
  const ConstTraxelRefVector& /*traxelrefs*/
) {
  throw std::runtime_error(
    "Feature extractor " + name() + " doesn't provide a vector valued result"
  );
  return *(new FeatureVector);
}
const FeatureScalar& SubsetFeatureExtractor::extract_scalar(
  const ConstTraxelRefVector& /*traxelrefs*/
) {
  throw std::runtime_error(
    "Feature extractor " + name() + " doesn't provide a scalar valued result"
  );
  return *(new FeatureScalar);
}

/*=============================================================================
  specific functors
=============================================================================*/
////
//// class TrackValue
////
const feature_type& TrackValue::operator()(const Track& track) {
  return (*feature_aggregator_) ( (*track_feature_extractor_) (track) );
}

////
//// class TrackingValue
////
TrackingValue::TrackingValue(
  SubsetsOfInterest* subsets_of_interest,
  SubsetFeatureExtractor* subset_feature_extractor,
  SubsetFeatureAggregator* subset_feature_aggregator,
  FeatureAggregator* feature_aggregator
) :
  subsets_of_interest_(subsets_of_interest),
  subset_feature_extractor_(subset_feature_extractor),
  subset_feature_aggregator_(subset_feature_aggregator),
  feature_aggregator_(feature_aggregator) {
}

// TODO Rewrite TrackingValue
const feature_type& TrackingValue::operator()(const Tracking& /*tracking*/) {
//  const std::vector<Trackvector>& subsets = (*subsets_of_interest_)(tracking);
//  feature_arrays feature_matrix;
//  ret_ = 0.0;
//  for (
//    std::vector<Trackvector>::const_iterator subset_it = subsets.begin();
//    subset_it != subsets.end();
//    subset_it++
//  ) {
//    feature_matrix.push_back(
//      (*subset_feature_aggregator_)((*subset_feature_extractor_)(*subset_it))
//    );
//  }
//  ret_ = (*feature_aggregator_)(feature_matrix);
  return ret_;
}

////
//// class SubsetFeaturesIdentity
////
const std::string SubsetFeaturesIdentity::name_ = "SubsetFeaturesIdentity";

SubsetFeaturesIdentity::SubsetFeaturesIdentity(
  const std::vector<std::string>& feature_names
) : feature_names_(feature_names) {
}

SubsetFeaturesIdentity::SubsetFeaturesIdentity(
  const std::string& feature_name
) {
  feature_names_.resize(1);
  feature_names_[0] = feature_name;
}

const std::string& SubsetFeaturesIdentity::name() const {
  return name_;
}

const FeatureMatrix& SubsetFeaturesIdentity::extract_matrix(
  const ConstTraxelRefVector& traxelrefs
) {
  // get the size of the return matrix
  size_t x_size = traxelrefs.size();

  // get y-size
  size_t y_size = 0;

  // get the feature map of the first traxel and iterate over the feature names
  const FeatureMap& feature_map = traxelrefs.front()->features;
  for (
    std::vector<std::string>::const_iterator fname_it = feature_names_.begin();
    fname_it != feature_names_.end();
    ++fname_it
  ) {
    FeatureMap::const_iterator feature_map_it = feature_map.find(*fname_it);
    if (feature_map_it != feature_map.end()) {
      y_size += feature_map_it->second.size();
    } else {
      LOG(logDEBUG) << "In SubsetFeaturesIdentity: Feature \"" << *fname_it << "\" not found";
    }
  }
  
  // initialize the return matrix
  ret_matrix_.reshape(vigra::Shape2(x_size, y_size));

  // iterate over all traxel
  size_t column_index = 0;
  for(
    ConstTraxelRefVector::const_iterator tref_it = traxelrefs.begin();
    tref_it != traxelrefs.end();
    tref_it++, column_index++
  ) {
    // fetch the features which are stored in a map
    const FeatureMap feature_map = (*tref_it)->features;

    // get the current column as a view
    FeatureVectorView column = ret_matrix_.bind<0>(column_index);

    // reset the row_index
    size_t row_index = 0;

    // iterate over all feature names
    for(
      std::vector<std::string>::const_iterator fname_it = feature_names_.begin();
      fname_it != feature_names_.end();
      fname_it++
    ) {
      // check if the feature name exists
      FeatureMap::const_iterator fmap_it = feature_map.find(*fname_it);
      if ( fmap_it != feature_map.end()) {
        // iterate over the elements of the feature vector
        for (
          feature_array::const_iterator f_it = fmap_it->second.begin();
          f_it != fmap_it->second.end();
          f_it++, row_index++
        ) {
          if (row_index >= y_size) {
            throw std::runtime_error(
              "In SubsetFeaturesIdentity: Invalid row index"
            );
          }
          column(row_index) = *f_it;
        }
      } // if (fmap_it != feature_map.end())"
    } // for (fname_it = ... )
  } // for(traxels_it = .. )
  return ret_matrix_;
}


////
//// class TrackFeaturesDiff
////
const std::string TrackFeaturesDiff::name_ = "TrackFeaturesDiff";

TrackFeaturesDiff::TrackFeaturesDiff(
  const std::vector<std::string>& feature_names
) : feature_names_(feature_names) {
}

TrackFeaturesDiff::TrackFeaturesDiff(
  const std::string& feature_name
) {
  feature_names_.resize(1);
  feature_names_[0] = feature_name;
}

const std::string& TrackFeaturesDiff::name() const {
  return name_;
}

const feature_arrays& TrackFeaturesDiff::operator()(
  const Track& track
) {
  ret_.clear();
  // iterate over all traxel
  if (track.traxels_.size() >= 2) {
    Traxelvector::const_iterator t1_it, t2_it;
    t1_it = track.traxels_.begin()+1;
    t2_it = track.traxels_.begin();
    for(; t1_it != track.traxels_.end(); t1_it++, t2_it++) {
      // fetch the features which are stored in a map
      const FeatureMap& fmap1 = t1_it->features;
      const FeatureMap& fmap2 = t2_it->features;
      feature_array feature_vector;
      // iterate over all feature names
      for(
        std::vector<std::string>::const_iterator fname_it = feature_names_.begin();
        fname_it != feature_names_.end();
        fname_it++
      ) {
        // check if the feature name exists
        FeatureMap::const_iterator fmap1_it = fmap1.find(*fname_it);
        FeatureMap::const_iterator fmap2_it = fmap2.find(*fname_it);
        assert(fmap1_it != fmap1.end());
        assert(fmap2_it != fmap2.end());
        // check if the feature vectors have the same size
        assert((fmap1_it->second).size() == (fmap2_it->second).size());
        // append the differences to the feature vector
        feature_array::const_iterator farray1_it, farray2_it;
        farray1_it = (fmap1_it->second).begin();
        farray2_it = (fmap2_it->second).begin();
        for (
          ;
          farray1_it != (fmap1_it->second).end();
          farray1_it++, farray2_it++
        ) {
          feature_vector.push_back(*farray1_it - *farray2_it);
        }
      }
      // append feature_vector to the feature matrix
      ret_.push_back(feature_vector);
    }
  }
  return ret_;
}

////
//// class TrackFeaturesCurvature
////
const std::string TrackFeaturesCurvature::name_ = "TrackFeaturesCurvature";

TrackFeaturesCurvature::TrackFeaturesCurvature(
  const std::vector<std::string>& feature_names
) : feature_names_(feature_names) {
}

TrackFeaturesCurvature::TrackFeaturesCurvature(
  const std::string& feature_name
) {
  feature_names_.resize(1);
  feature_names_[0] = feature_name;
}

const std::string& TrackFeaturesCurvature::name() const {
  return name_;
}

const feature_arrays& TrackFeaturesCurvature::operator()(
  const Track& track
) {
  ret_.clear();
  // iterate over all traxel
  assert(track.traxels_.size() >= 3);
  Traxelvector::const_iterator t1_it, t2_it, t3_it;
  t1_it = track.traxels_.begin()+2;
  t2_it = track.traxels_.begin()+1;
  t3_it = track.traxels_.begin();
  for(; t1_it != track.traxels_.end(); t1_it++, t2_it++, t3_it++) {
    // fetch the features which are stored in a map
    const FeatureMap& fmap1 = t1_it->features;
    const FeatureMap& fmap2 = t2_it->features;
    const FeatureMap& fmap3 = t3_it->features;
    feature_array feature_vector;
    // iterate over all feature names
    for(
      std::vector<std::string>::const_iterator fname_it = feature_names_.begin();
      fname_it != feature_names_.end();
      fname_it++
    ) {
      // check if the feature name exists
      FeatureMap::const_iterator fmap1_it = fmap1.find(*fname_it);
      FeatureMap::const_iterator fmap2_it = fmap2.find(*fname_it);
      FeatureMap::const_iterator fmap3_it = fmap3.find(*fname_it);
      assert(fmap1_it != fmap1.end());
      assert(fmap2_it != fmap2.end());
      assert(fmap3_it != fmap3.end());
      // check if the feature vectors have the same size
      assert((fmap1_it->second).size() == (fmap2_it->second).size());
      assert((fmap1_it->second).size() == (fmap3_it->second).size());
      // append the curvature to the feature vector
      feature_array::const_iterator farray1_it, farray2_it, farray3_it;
      farray1_it = (fmap1_it->second).begin();
      farray2_it = (fmap2_it->second).begin();
      farray3_it = (fmap3_it->second).begin();
      for (
        ;
        farray1_it != (fmap1_it->second).end();
        farray1_it++, farray2_it++, farray3_it++
      ) {
        feature_vector.push_back(*farray1_it - 2 * *farray2_it + *farray3_it);
      }
    }
    // append feature_vector to the feature matrix
    ret_.push_back(feature_vector);
  }
  return ret_;
}

////
//// class SumAggregator
////
const std::string SumAggregator::name_ = "SumAggregator";

const std::string& SumAggregator::name() const {
  return name_;
}

const feature_type& SumAggregator::operator()(
  const feature_arrays& features
) {
  ret_ = 0;
  for (
    feature_arrays::const_iterator farray = features.begin();
    farray != features.end();
    farray++
  ) {
    for (
      feature_array::const_iterator f = farray->begin();
      f != farray->end();
      f++
    ) {
      ret_ += *f;
    }
  }
  return ret_;
}

////
//// class TrackSubsets
////
const std::string TrackSubsets::name_ = "TrackSubsets";

const std::string& TrackSubsets::name() const {
  return name_;
}

const std::vector<ConstTraxelRefVector>& TrackSubsets::operator()(
  const HypothesesGraph& graph
) {
  ret_.clear();

  // Check if the graph has the necessary attributes
  if (not graph.has_property(node_active())) {
    throw std::runtime_error(
      "Graph doesn't have a \"node_active\" property map"
    );
  }
  if (not graph.has_property(arc_active())) {
    throw std::runtime_error(
      "Graph doesn't have an \"arc_active\" property map"
    );
  }

  // check if we have a tracklet graph
  bool with_tracklets = graph.has_property(node_tracklet());

  // check if the graph is legal
  if (not (graph.has_property(node_traxel()) or with_tracklets)) {
    throw std::runtime_error(
      "HypothesesGraph has neither traxel nor tracklet property map"
    );
  }

  // Get the property maps
  node_active_map_type& node_active_map = graph.get(node_active());
  arc_active_map_type& arc_active_map = graph.get(arc_active());

  // Make maps from child to parent and parent to child
  typedef std::map<HypothesesGraph::Node, HypothesesGraph::Node> NodeNodeMap;
  NodeNodeMap parent;
  NodeNodeMap child;
  // Initialize
  for (NodeActiveIt n_it(node_active_map); n_it != lemon::INVALID; ++n_it) {
    parent[n_it] = n_it;
    child[n_it] = n_it;
  }

  // Set connections
  for (ArcActiveIt a_it(arc_active_map); a_it != lemon::INVALID; ++a_it) {
    // count the active arcs with the same source
    size_t out_arcs = 0;
    for (
      OutArcIt o_it(graph, graph.source(a_it));
      o_it != lemon::INVALID;
      ++o_it
    ) {
      out_arcs += (arc_active_map[o_it] ? 1 : 0);
    }
    // link those nodes if there are no other active arcs with the same source
    if (out_arcs == 1) {
      parent[graph.target(a_it)] = graph.source(a_it);
      child[graph.source(a_it)] = graph.target(a_it);
    }
  }
  
  // Compose return vector of traxel reference vectors
  for (
    NodeNodeMap::const_iterator nmap_it = parent.begin();
    nmap_it != parent.end();
    ++nmap_it
  ) {
    HypothesesGraph::Node current_node = nmap_it->first;
    LOG(logDEBUG4) << "Is parent node invalid?";
    LOG(logDEBUG4) << (parent[current_node] == lemon::INVALID);
    if (parent[current_node] == current_node) {
      // resize the return vector
      ret_.resize(ret_.size()+1);
      bool loop = true;
      // loop as long as the track isn't finished
      while (loop) {
        if (with_tracklets) {
          // get the traxel vector of this node
          const std::vector<Traxel>& t_vec = graph.get(node_tracklet())[current_node];
          for (
            std::vector<Traxel>::const_iterator t_it = t_vec.begin();
            t_it != t_vec.end();
            t_it++
          ) {
            ret_.back().push_back( &(*t_it) );
          }
        } else {
          ret_.back().push_back( &(graph.get(node_traxel())[current_node]) );
        }
        loop = current_node != child[current_node];
        current_node = child[current_node];
      }
    }
  }
  return ret_;
}


////
//// class DivisionSubsets
////
const std::string DivisionSubsets::name_ = "DivisionSubsets";

const std::string& DivisionSubsets::name() const {
  return name_;
}

const std::vector<ConstTraxelRefVector>& DivisionSubsets::operator()(
  const HypothesesGraph& graph
) {
  return operator()(graph, 1);
}
const std::vector<ConstTraxelRefVector>& DivisionSubsets::operator()(
  const HypothesesGraph& graph,
  size_t depth
) {
  ret_.clear();

  // Check if the graph has the necessary attributes
  if (not graph.has_property(node_active())) {
    throw std::runtime_error(
      "Graph doesn't have a \"node_active\" property map"
    );
  }
  if (not graph.has_property(arc_active())) {
    throw std::runtime_error(
      "Graph doesn't have an \"arc_active\" property map"
    );
  }

  // check if we have a tracklet graph
  bool with_tracklets = graph.has_property(node_tracklet());

  // check if the graph is legal
  if (not (graph.has_property(node_traxel()) or with_tracklets)) {
    throw std::runtime_error(
      "HypothesesGraph has neither traxel nor tracklet property map"
    );
  }

  // Get the property maps
  node_active_map_type& node_active_map = graph.get(node_active());
  arc_active_map_type& arc_active_map = graph.get(arc_active());

  // Find the divisions
  for (NodeActiveIt n_it(node_active_map); n_it != lemon::INVALID; ++n_it) {
    // Count the active outgoing arcs
    std::vector<HypothesesGraph::Arc> out_arcs;
    for (OutArcIt o_it(graph, n_it); o_it != lemon::INVALID; ++o_it) {
      if (arc_active_map[o_it]) {
        out_arcs.push_back(o_it);
      }
    }
    // Two outgoing arcs: division
    if (out_arcs.size() == 2) {
      if (with_tracklets) {
        // %TODO
      } else {
        // %TODO
      }
    }
  }
  return ret_;
}

////
//// class SubsetFeatureExtractorFromFE
////
const std::string SubsetFeatureExtractorFromFE::name_ = "SubsetFeatureExtractorFromFE";

const std::string& SubsetFeatureExtractorFromFE::name() const {
  return name_;
}

const feature_arrays& SubsetFeatureExtractorFromFE::operator()(
  const Trackvector& tracks
) {
  return (*track_feature_extractor_)(tracks.front());
}

////
//// class DivisionFeatureExtractor
////
// const std::string DivisionFeatureExtractor::name_ = "DivisionFeatureExtractor";
// 
// const std::string& DivisionFeatureExtractor::name() const {
//   return name_;
// }
// 
// DivisionFeatureExtractor::DivisionFeatureExtractor(
//   const std::string& feature_name,
//   size_t depth
// ) : features_identity_(feature_name), depth_(depth) {
// }
//  
// DivisionFeatureExtractor::DivisionFeatureExtractor(
//   const std::vector<std::string>& feature_names,
//   size_t depth
// ) : features_identity_(feature_names), depth_(depth) {
// }
// 
// const feature_arrays& DivisionFeatureExtractor::operator()(
//     const Trackvector& tracks
// ) {
//   Track tmp;
//   const Track& parent = tracks[0];
//   const Track& left_child = tracks[1];
//   const Track& right_child = tracks[2];
//   if(
//     (parent.get_length() >= depth_) and 
//     (left_child.get_length() >= depth_) and
//     (right_child.get_length() >= depth_)
//   ) {
//     Traxelvector::const_iterator p_it = parent.traxels_.end() - depth_;
//     Traxelvector::const_iterator l_it = left_child.traxels_.begin();
//     Traxelvector::const_iterator r_it = right_child.traxels_.begin();
// 
//     tmp.traxels_.resize(depth_ * 3);
//     for (size_t i = 0; i < depth_; i++, p_it++, l_it++, r_it++) {
//       tmp.traxels_[i] = *p_it;
//       tmp.traxels_[i+depth_] = *l_it;
//       tmp.traxels_[i+2*depth_] = *r_it;
//     }
//   }
//   return features_identity_(tmp);
// }

////
//// class SubsetAggregatorFromFA
////
const std::string SubsetAggregatorFromFA::name_ = "SubsetAggregatorFromFA";

const std::string& SubsetAggregatorFromFA::name() const {
  return name_;
}

const feature_array& SubsetAggregatorFromFA::operator()(
  const feature_arrays& features
) {
  ret_.clear();
  ret_.push_back((*feature_aggregator_)(features));
  return ret_;
}

////
//// class ChildRatioAggregator
////
const std::string ChildRatioAggregator::name_ = "ChildRatioAggregator";

const std::string& ChildRatioAggregator::name() const {
  return name_;
}

const feature_array& ChildRatioAggregator::operator()(
  const feature_arrays& features
) {
  ret_.resize(features.front().size());
  feature_array::iterator ret_it = ret_.begin();

  feature_arrays::const_iterator left_it = features.begin() + depth_;
  feature_arrays::const_iterator right_it = features.begin() + 2 * depth_;
  feature_array::const_iterator l_val_it, r_val_it;

  feature_array l_sum(features.front().size(), 0);
  feature_array r_sum(features.front().size(), 0);
  feature_array::iterator l_sum_it, r_sum_it;

  for(size_t i = 0; i < depth_; i++, left_it++, right_it++, ret_it++) {
    l_val_it = left_it->begin();
    r_val_it = right_it->begin();
    l_sum_it = l_sum.begin();
    r_sum_it = r_sum.begin();
    for(
      ;
      l_val_it != left_it->end();
      l_val_it++, r_val_it++, l_sum_it++, r_sum_it++
    ) {
      *l_sum_it += *l_val_it;
      *r_sum_it += *r_val_it;
    }
  }
  l_sum_it = l_sum.begin();
  r_sum_it = r_sum.begin();
  for(; l_sum_it != l_sum.end(); l_sum_it++, r_sum_it++, ret_it++) {
    if (*l_sum_it > *r_sum_it) {
      *ret_it = *r_sum_it / *l_sum_it;
    } else {
      *ret_it = *l_sum_it / *r_sum_it;
    }
  }
  return ret_;
}

////
//// class OutlierCountAggregator
////
const std::string OutlierCountAggregator::name_ = "OutlierCountAggregator";

const std::string& OutlierCountAggregator::name() const {
  return name_;
}

const feature_type& OutlierCountAggregator::operator()(
  const feature_arrays& features
) {
  std::vector<size_t> outlier_ids = outlier_calculator_->calculate(features);
  ret_ = static_cast<feature_type>(outlier_ids.size());
  // prevent divide by 0 error
  if (features.size()) {
    ret_ = ret_ / static_cast<feature_type>(features.size());
  }
  return ret_;
}

////
//// class OutlierBadnessAggregator
////
const std::string OutlierBadnessAggregator::name_ = "OutlierBadnessAggregator";

const std::string& OutlierBadnessAggregator::name() const {
  return name_;
}

const feature_type& OutlierBadnessAggregator::operator()(
  const feature_arrays& features
) {
  outlier_calculator_->calculate(features);
  feature_array outlier_badness = outlier_calculator_->get_measures();
  ret_ = 0.0;
  for (
    feature_array::const_iterator f_it = outlier_badness.begin();
    f_it != outlier_badness.end();
    f_it++
  ) {
    ret_ = *f_it > ret_ ? *f_it : ret_;
  }
  return ret_;
}


////
//// function to_arma_matrix
////
arma::Mat<feature_type> to_arma_matrix(const feature_arrays& features) {
  assert(features.empty() == false);
  size_t cols = features.size();
  size_t rows = features[0].size();
  arma::Mat<feature_type> ret(rows, cols);

  typename std::vector<feature_array>::const_iterator feature_array_it;
  feature_array_it = features.begin();
  for (
    size_t j=0;
    feature_array_it != features.end();
    feature_array_it++, j++
  ) {
    assert(feature_array_it->size() == rows);
    arma::Col<feature_type> column(*feature_array_it);
    ret.col(j) = column;
  }
  return ret;
}

////
//// class OutlierCalculator
////
const std::string OutlierCalculator::name_ = "OutlierCalculator";

const std::string& OutlierCalculator::name() const {
  return name_;
}

const feature_array& OutlierCalculator::get_measures() const {
  throw std::runtime_error(
    "OutlierCalculator \"" + name() + "\"doesn't provide a measure"
  );
  return measures_;
}

////
//// class MVNOutlierCalculator
////
const std::string MVNOutlierCalculator::name_ = "MVNOutlierCalculator";

MVNOutlierCalculator::MVNOutlierCalculator(const feature_type sigma_threshold) {
  sigma_threshold_ = sigma_threshold;
}

const std::string& MVNOutlierCalculator::name() const {
  return name_;
}

const feature_array& MVNOutlierCalculator::get_measures() const {
  return measures_;
}

const arma::Mat<feature_type>& MVNOutlierCalculator::get_covariance() const {
  return covariance_;
}

const arma::Mat<feature_type>& MVNOutlierCalculator::get_inverse_covariance() const {
  return inv_covariance_;
}

const arma::Col<feature_type>& MVNOutlierCalculator::get_mean() const {
  return mean_;
}

const std::vector<size_t>& MVNOutlierCalculator::calculate(
  const feature_arrays& features
) {
  measures_.clear();
  outlier_ids_.clear();
  mean_.clear();
  covariance_.clear();
  inv_covariance_.clear();

  bool good_data = true;
  if (features.size() == 0) {
    good_data = false;
  } else if (features.size() <= features[0].size()+1) {
    good_data = false;
  } 
  if (good_data) {
    // Get covariance and inverse covariance matrix
    arma::Mat<feature_type> features_mat(to_arma_matrix(features));
    arma::Mat<feature_type> features_mat_t(trans(features_mat));
    covariance_ = arma::cov(features_mat_t);
    bool invertible = arma::inv_sympd(inv_covariance_, covariance_);

    if (invertible) {
      // Get mean values
      mean_ = arma::mean(features_mat, 1);
  
      // Calculate the outliers
      outlier_ids_.clear();
      measures_.clear();
      feature_arrays::const_iterator features_it = features.begin();
      for(size_t id=0; features_it != features.end(); features_it++, id++) {
        arma::Col<feature_type> diff_vector(*features_it);
        diff_vector -= mean_;
        feature_type norm_residual = arma::dot(diff_vector, inv_covariance_*diff_vector);
        measures_.push_back(norm_residual);
        if (norm_residual >= sigma_threshold_) {
          outlier_ids_.push_back(id);
        }
      }
    }
  } // else
  return outlier_ids_;
}

} // namespace pgmlink
