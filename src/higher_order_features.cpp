#include "pgmlink/higher_order_features.h"

namespace pgmlink {

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
  const size_t index,
  HypothesesGraph::Node& last_node
) {
  bool tracklet_graph = graph.has_property(node_tracklet());
  node_timestep_type& timestep_map = graph.get(node_timestep());
  arc_active_map_type& arc_active_map = graph.get(arc_active_count());

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
      if(arc_active_map[oa_it][index]){
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
  std::map<HypothesesGraph::Node, Track*> node_in_track;
  // stores the parent node of a track with the track as the key
  std::map<Track*, HypothesesGraph::Node> has_parent_node;

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
            if(arc_active_map[oa_it][index_]) active_out_arcs++;
          }
          if (active_out_arcs==2) parent_division_active=true;
        }
      } // end if (num_in==1)
      if ((num_in==0) or (parent_division_active)) {
        // initialize the track id and the variable in which the last node will
        // be stored
        HypothesesGraph::Node last_node;
        size_t track_id = tracks_.size();
        // create the track and get the last node of it
        tracks_.push_back(track_from_start_node(graph, n_it, index, last_node));
        // get reference to track and set its id
        Track* track_ref = &(tracks_.back());
        track_ref->set_id(track_id);
        // update the maps "node_in_track" and "has_parent_node"
        node_in_track[last_node] = track_ref;
        has_parent_node[track_ref] = parent;
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
    HypothesesGraph::Node parent = has_parent_node[&(*t_it)];
    if (parent != lemon::INVALID) {
      // get the map iterator to the parent track and check if it exists
      std::map<HypothesesGraph::Node, Track*>::iterator parent_track_it;
      parent_track_it = node_in_track.find(parent);
      if (parent_track_it == node_in_track.end()) {
        LOG(logDEBUG) << "In constructor of tracking: inconsistent or unsupported graph";
      } else {
        // set the references
        Track* parent_track_ref = parent_track_it->second;
        t_it->set_parent(parent_track_ref);
        parent_track_ref->set_child(&(*t_it));
      }
    }
  }
}

/*=============================================================================
  specific functors
=============================================================================*/
////
//// class TrackValue
////
feature_type TrackValue::operator()(const Track& track) {
  return (*feature_aggregator_) ( (*track_feature_extractor_) (track) );
}

////
//// class TrackFeaturesIdentity
////
const std::string TrackFeaturesIdentity::name_ = "TrackFeaturesIdentity";

TrackFeaturesIdentity::TrackFeaturesIdentity(
  const std::vector<std::string>& feature_names
) : feature_names_(feature_names) {
}

TrackFeaturesIdentity::TrackFeaturesIdentity(
  const std::string& feature_name
) {
  feature_names_.resize(1);
  feature_names_[0] = feature_name;
}

const std::string& TrackFeaturesIdentity::name() const {
  return name_;
}

feature_arrays TrackFeaturesIdentity::operator()(
  const Track& track
) {
  feature_arrays feature_matrix;
  // iterate over all traxel
  for(
    Traxelvector::const_iterator traxel_it = track.traxels_.begin();
    traxel_it != track.traxels_.end();
    traxel_it++
  ) {
    // fetch the features which are stored in a map
    FeatureMap feature_map = traxel_it->features;
    feature_array feature_vector;
    // iterate over all feature names
    for(
      std::vector<std::string>::const_iterator fname_it = feature_names_.begin();
      fname_it != feature_names_.end();
      fname_it++
    ) {
      // assert if the feature name exists
      FeatureMap::const_iterator f = feature_map.find(*fname_it);
      assert(f != feature_map.end());
      // append the features to the feature vector
      feature_vector.insert(
        feature_vector.end(),
        (f->second).begin(),
        (f->second).end()
      );
    }
    // append feature_vector to the feature matrix
    feature_matrix.push_back(feature_vector);
  }
  return feature_matrix;
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

feature_arrays TrackFeaturesDiff::operator()(
  const Track& track
) {
  feature_arrays feature_matrix;
  // iterate over all traxel
  assert(track.traxels_.size() >= 2);
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
    feature_matrix.push_back(feature_vector);
  }
  return feature_matrix;
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

feature_arrays TrackFeaturesCurvature::operator()(
  const Track& track
) {
  feature_arrays feature_matrix;
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
    feature_matrix.push_back(feature_vector);
  }
  return feature_matrix;
}

////
//// class DivisionSubsets
////
const std::string DivisionSubsets::name_ = "DivisionSubsets";

const std::string& DivisionSubsets::name() const {
  return name_;
}

std::vector<Trackvector> DivisionSubsets::operator()(const Tracking& tracking) {
  const Trackvector& tracks = tracking.tracks_;
  std::vector<Trackvector> ret;
  for(
    Trackvector::const_iterator t_it = tracks.begin();
    t_it != tracks.end();
    t_it++
  ) {
    if ((t_it->children_).size() == 2) {
      Trackvector t;
      t.push_back(*t_it);
      t.push_back(*((t_it->children_)[0]));
      t.push_back(*((t_it->children_)[1]));
      ret.push_back(t);
    }
  }
  return ret;
}

////
//// class SubsetAggregatorFromFA
////
const std::string SubsetAggregatorFromFA::name_ = "SubsetAggregatorFromFA";

const std::string& SubsetAggregatorFromFA::name() const {
  return name_;
}

feature_array SubsetAggregatorFromFA::operator()(
  const feature_arrays& features
) {
  feature_array f(1, (*feature_aggregator_)(features));
  return f;
}

////
//// class OutlierCountAggregator
////
const std::string OutlierCountAggregator::name_ = "OutlierCountAggregator";

const std::string& OutlierCountAggregator::name() const {
  return name_;
}

feature_type OutlierCountAggregator::operator()(
  const feature_arrays& features
) {
  std::vector<size_t> outlier_ids = outlier_calculator_->calculate(features);
  return static_cast<feature_type>(outlier_ids.size());
}

////
//// class OutlierBadnessAggregator
////
const std::string OutlierBadnessAggregator::name_ = "OutlierBadnessAggregator";

const std::string& OutlierBadnessAggregator::name() const {
  return name_;
}

feature_type OutlierBadnessAggregator::operator()(
  const feature_arrays& features
) {
  outlier_calculator_->calculate(features);
  feature_array outlier_badness = outlier_calculator_->get_measures();
  feature_type ret = 0.0;
  for (
    feature_array::const_iterator f_it = outlier_badness.begin();
    f_it != outlier_badness.end();
    f_it++
  ) {
    ret = *f_it > ret ? *f_it : ret;
  }
  return ret;
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
  if (features.size() <= features[0].size()) {
    measures_.clear();
    outlier_ids_.clear();
    mean_.clear();
    covariance_.clear();
    inv_covariance_.clear();
  } else {
    // Get covariance and inverse covariance matrix
    arma::Mat<feature_type> features_mat(to_arma_matrix(features));
    arma::Mat<feature_type> features_mat_t(trans(features_mat));
    try {
      covariance_ = arma::cov(features_mat_t);
      inv_covariance_ = arma::inv_sympd(covariance_);

      LOG(logDEBUG4) << "In MVNOutlierCalculator: covariance matrix";
      LOG(logDEBUG4) << covariance_;
      LOG(logDEBUG4) << "In MVNOutlierCalculator: inverse covariance matrix";
      LOG(logDEBUG4) << inv_covariance_;

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
    } catch (std::exception& exception) {
      LOG(logDEBUG) << "In MVNOutlierCalculator: Too few data to calculate outliers" << std::endl;
    }
  } // else
  return outlier_ids_;
}

} // namespace pgmlink
