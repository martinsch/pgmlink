#include "pgmlink/tracks.h"

namespace pgmlink {

////
//// Class Track
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
//// Class Tracking
////
Tracking::Tracking() {
}

Tracking::Tracking(const Events& events, const TraxelStore& traxelstore) 
: events_(events), traxelstore_(traxelstore) {
  Tracking::build_tracks();
}

void Tracking::build_tracks() {
  Events::const_iterator events_it = events_.begin();
  for(size_t t = 0; events_it != events_.begin(); events_it++, t++) {
    Tracking::apply_events(*events_it, t);
  }
}

void Tracking::apply_events(
  const std::vector<Event>& events,
  const size_t time_step
) {
  for(
    typename std::vector<Event>::const_iterator event_it = events.begin();
    event_it != events.end();
    event_it++
  ) {
    Tracking::apply_event(*event_it, time_step);
  }
}

size_t Tracking::start_new_track(
  const size_t time_start,
  const size_t track_parent_id
) {
  size_t track_id = tracks_.size();
  Track new_track;
  new_track.set_id(track_id);
  new_track.set_time_start(time_start);
  new_track.set_parent_id(track_parent_id);

  tracks_.push_back(new_track);
  return track_id;
}

// TODO
void Tracking::apply_event(const Event& event, const size_t time_step) {
  switch (event.type) {
    case Event::Move: {
        assert(time_step > 0);
        size_t parent_traxel_id = event.traxel_ids[0];
        size_t child_traxel_id = event.traxel_ids[1];
        Traxel parent; // TODO get it out of the TraxelStore
        Traxel child; // TODO get it out of the TraxelStore
        size_t track_id;
        if (in_track_[time_step-1].count(parent_traxel_id)) {
          track_id = in_track_[time_step-1][parent_traxel_id];
        } else {
          track_id = start_new_track(time_step);
        }
        tracks_[track_id].traxels_.push_back(child);
        in_track_[time_step][child_traxel_id] = track_id;
      }
      break;
    case Event::Division: {
        assert(time_step > 0);
        size_t parent_traxel_id = event.traxel_ids[0];
        size_t lchild_traxel_id = event.traxel_ids[1];
        size_t rchild_traxel_id = event.traxel_ids[2];
  
        Traxel parent; // TODO get it out of the TraxelStore
        Traxel lchild; // TODO get it out of the TraxelStore
        Traxel rchild; // TODO get it out of the TraxelStore
        
        size_t ptrack_id=0;
        size_t ltrack_id, rtrack_id;
        if (in_track_[time_step-1].count(parent_traxel_id)) {
          ptrack_id = in_track_[time_step-1][parent_traxel_id];
          tracks_[ptrack_id].child_ids_[0] = tracks_.size();
          tracks_[ptrack_id].child_ids_[1] = tracks_.size() + 1;
        }
        ltrack_id = start_new_track(time_step, ptrack_id);
        rtrack_id = start_new_track(time_step, ptrack_id);
  
        tracks_[ltrack_id].traxels_.push_back(lchild);
        tracks_[rtrack_id].traxels_.push_back(rchild);

        in_track_[time_step][rchild_traxel_id] = rtrack_id;
        in_track_[time_step][lchild_traxel_id] = ltrack_id;
      }
      break;
    case Event::Appearance: {
        size_t traxel_id = event.traxel_ids[0];
        
        Traxel traxel; // TODO get it out of the TraxelStore

        size_t track_id = start_new_track(time_step);
        tracks_[track_id].traxels_.push_back(traxel);

        in_track_[time_step][traxel_id] = track_id;
      }
      break;
    case Event::Disappearance:
      break;
    case Event::Merger:
      break;
    case Event::ResolvedTo:
      break;
    case Event::MultiFrameMove:
      break;
  }
}

////
//// Class TrackFeatureExtractor
////
TrackFeatureExtractor::TrackFeatureExtractor(
  boost::shared_ptr<FeatureCalculator> calculator,
  const std::string& feature_name,
  const FeatureOrder order 
) :
  calculator_(calculator),
  extractor_(new FeatureExtractor(calculator, feature_name)),
  feature_name_(feature_name),
  order_(order) {
}

TrackFeatureExtractor::TrackFeatureExtractor(
  boost::shared_ptr<FeatureCalculator> calculator,
  boost::shared_ptr<FeatureAggregator> aggregator,
  const std::string& feature_name,
  const FeatureOrder order
) :
  calculator_(calculator),
  extractor_(new FeatureExtractor(calculator, feature_name)),
  aggregator_(aggregator),
  feature_name_(feature_name),
  order_(order) {
}

feature_arrays TrackFeatureExtractor::extract(const Track& track) const {
  assert(track.get_length - order_ > 0);

  feature_arrays ret;
  Traxelvector::const_iterator traxel_it;

  for (
    traxel_it = track.traxels_.begin();
    traxel_it+order_ != track.traxels_.end();
    traxel_it++
  ) {
    switch(order_) {
      case SINGLE:
        ret.push_back(extractor_->extract(*traxel_it));
        break;
      case PAIRWISE:
        ret.push_back(extractor_->extract(*traxel_it, *(traxel_it+1)));
        break;
      case TRIPLET:
        ret.push_back(extractor_->extract(*traxel_it, *(traxel_it+1), *(traxel_it+2)));
    }
  }
  return ret;
}

feature_array TrackFeatureExtractor::extract_vector(const Track& track) const {
  if (!aggregator_) {
    throw std::runtime_error (
      "No aggregator set in TrackFeatureExtractor"
    );
  }
  feature_arrays features = TrackFeatureExtractor::extract(track);
  return aggregator_->vector_valued(features);
}

feature_type TrackFeatureExtractor::extract_scalar(const Track& track) const {
  if (!aggregator_) {
    throw std::runtime_error (
      "No aggregator set in TrackFeatureExtractor"
    );
  }
  feature_arrays features = TrackFeatureExtractor::extract(track);
  return aggregator_->scalar_valued(features);
} 

FeatureOrder TrackFeatureExtractor::get_feature_order() const {
  return TrackFeatureExtractor::order_;
}

////
//// Class TrackingFeatureExtractor
////
TrackingFeatureExtractor::TrackingFeatureExtractor(
  boost::shared_ptr<TrackFeatureExtractor> track_feature_extractor,
  boost::shared_ptr<FeatureAggregator> feature_aggregator,
  TrackFeatureOrder track_feature_order
) : track_feature_extractor_(track_feature_extractor),
  feature_aggregator_(feature_aggregator),
  track_feature_order_(track_feature_order)
{

}

feature_type TrackingFeatureExtractor::extract(const Tracking& tracking) const {
  feature_type ret = 0.0;
  switch (track_feature_order_) {
    case VECTOR: {
      feature_arrays features(tracking.tracks_.size());
      feature_arrays::iterator f_it = features.begin();
      Trackvector::const_iterator tracks_it = tracking.tracks_.begin();
      for (; tracks_it != tracking.tracks_.end(); tracks_it++, f_it++) {
        *f_it = track_feature_extractor_->extract_vector(*tracks_it);
      }
      ret = feature_aggregator_->scalar_valued(features);
    }
    break;
    case SCALAR: {
      feature_array features(tracking.tracks_.size());
      feature_array::iterator f_it = features.begin();
      Trackvector::const_iterator tracks_it = tracking.tracks_.begin();
      for(; tracks_it != tracking.tracks_.end(); tracks_it++, f_it++) {
        *f_it = track_feature_extractor_->extract_scalar(*tracks_it);
      }
      ret = feature_aggregator_->scalar_valued(features);
    }
    break;
  }
  return ret;
}
  


} // namespace pgmlink
