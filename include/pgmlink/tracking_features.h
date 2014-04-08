/**
* @file
* @ingroup tracking
* @brief track datastructures
*
* This file provides the datastructure track and a feature calculator operating
* on tracks
*/
#ifndef PGMLINK_TRACKS_H
#define PGMLINK_TRACKS_H

// stl
#include <vector>

// pgmlink
#include "pgmlink/traxels.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/event.h"
#include "pgmlink/classifier_auxiliary.h"
#include "pgmlink/tracking_features_auxiliary.h"

// boost
#include <boost/shared_ptr.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

namespace pgmlink {

// forward declaration of class Track
class Track;

// typedefs
typedef std::vector<Traxel> Traxelvector;
typedef std::vector<Track> Trackvector;

enum FeatureOrder{
  SINGLE,
  PAIRWISE,
  TRIPLET
};

enum TrackFeatureOrder{
  VECTOR,
  SCALAR
};

////
//// class Track
////
class Track {
 public:
  Track();
  Track(
    const size_t id,
    const size_t time_start,
    const Traxelvector& traxels,
    const size_t parent_id = 0,
    const size_t left_child_id = 0,
    const size_t right_child_id = 0
  );

  void set_id(const size_t id);
  size_t get_id() const;

  void set_time_start(const size_t time_start);
  size_t get_time_start() const;

  void set_parent_id(const size_t parent_id);
  size_t get_parent_id() const;

  void set_child_ids(const size_t left, const size_t right);
  const std::vector<size_t>& get_child_ids() const;
  
  size_t get_length() const;

  size_t id_;
  size_t time_start_;
  Traxelvector traxels_;
  size_t parent_id_;
  std::vector<size_t> child_ids_;
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & archive, const unsigned int version) {
      (void) version;
    archive & id_;
    archive & time_start_;
    archive & traxels_;
    archive & parent_id_;
    archive & child_ids_;
  }
}; // class Track

////
//// Class Tracking
//// TODO function apply_events, error handling?
class Tracking {
 public:
  typedef std::vector<std::vector<Event> > Events;
  Tracking();
  Tracking(const HypothesesGraph& graph, const size_t index);
  Tracking(const Events& events, const TraxelStore& traxelstore);
  void build_tracks(const Events& events, const TraxelStore& traxelstore);

  Trackvector tracks_;
 private:
  void build_tracks();
  void apply_event(const Event& event, const size_t timestep);
  void apply_events(const std::vector<Event>& events, const size_t timestep);
  size_t start_new_track(const size_t time_start, const size_t parent_id=0);
  size_t index_;
  Events events_;
  TraxelStore traxelstore_;
  std::map<size_t, std::map<size_t, size_t> > in_track_;
}; // Class Tracking

////
//// class TrackFeatureExtractor
////
class TrackFeatureExtractor {
 public:
  TrackFeatureExtractor(
    boost::shared_ptr<FeatureCalculator> calculator,
    const std::string& feature_name,
    const FeatureOrder order
  );
  TrackFeatureExtractor(
    boost::shared_ptr<FeatureCalculator> calculator,
    boost::shared_ptr<FeatureAggregator> aggregator,
    const std::string& feature_name,
    const FeatureOrder order
  );
  virtual ~TrackFeatureExtractor() {};
  virtual feature_arrays extract(const Track& track) const;
  virtual feature_array extract_vector(const Track& track) const;
  virtual feature_type extract_scalar(const Track& track) const;
  FeatureOrder get_feature_order() const;

 protected:
  boost::shared_ptr<FeatureCalculator> calculator_;
  boost::shared_ptr<FeatureExtractor> extractor_;
  boost::shared_ptr<FeatureAggregator> aggregator_;
  std::string feature_name_;
  FeatureOrder order_;
}; // class TrackFeatureExtractor


////
//// Class TrackingFeatureExtractor
////
class TrackingFeatureExtractor {
 public:
  TrackingFeatureExtractor(
    boost::shared_ptr<TrackFeatureExtractor>,
    boost::shared_ptr<FeatureAggregator>,
    TrackFeatureOrder track_feature_order
  );
  virtual ~TrackingFeatureExtractor() {};
  virtual feature_type extract(const Tracking& tracking) const;
  TrackFeatureOrder get_track_feature_order() const;

  boost::shared_ptr<TrackFeatureExtractor> track_feature_extractor_;
  boost::shared_ptr<FeatureAggregator> feature_aggregator_;
  TrackFeatureOrder track_feature_order_;
};

////
//// Class OutlierCount
////
class OutlierCount {
 public:
  OutlierCount(const std::string& feature_name = "com");
  ~OutlierCount() {};
  size_t operator()(const Tracking& tracking) const;
  size_t operator()(const Track& track) const;
 protected:
  boost::shared_ptr<FeatureCalculator> feature_identity_;
  boost::shared_ptr<FeatureAggregator> outlier_count_aggregator_;
  boost::shared_ptr<FeatureAggregator> sum_aggregator_;
  boost::shared_ptr<TrackFeatureExtractor> track_outlier_count_;
  boost::shared_ptr<TrackingFeatureExtractor> tracking_outlier_count_;
}; // Class OutlierCount

////
//// Class DiffOutlierCount
////
class DiffOutlierCount {
 public:
  DiffOutlierCount(const std::string& feature_name = "com");
  ~DiffOutlierCount() {};
  size_t operator()(const Tracking& tracking) const;
  size_t operator()(const Track& track) const;
 protected:
  boost::shared_ptr<FeatureCalculator> feature_difference_;
  boost::shared_ptr<FeatureAggregator> outlier_count_aggregator_;
  boost::shared_ptr<FeatureAggregator> sum_aggregator_;
  boost::shared_ptr<TrackFeatureExtractor> track_outlier_count_;
  boost::shared_ptr<TrackingFeatureExtractor> tracking_outlier_count_;
}; // Class DiffOutlierCount

} // namespace pgmlink

#endif // PGMLINK_TRACKS_H
