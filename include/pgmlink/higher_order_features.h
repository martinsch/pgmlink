/**
* @file
* @ingroup tracking
* @brief calculation of higher order features
*
* This file provides an interface to calculate higher order features from a
* tracking.
*/

#ifndef PGMLINK_HIGHER_ORDER_FEATURES_H
#define PGMLINK_HIGHER_ORDER_FEATURES_H

// stl
#include <vector>

// pgmlink
#include "pgmlink/traxels.h" /* for traxels */
#include "pgmlink/hypotheses.h" /* for hypotheses graph */

// boost
#include <boost/serialization/serialization.hpp> /* for serialization */

namespace pgmlink {

// forward declaration of class Track
class Track;

/*=============================================================================
 type definitions
=============================================================================*/
typedef std::vector<Traxel> Traxelvector;
typedef std::vector<Track> Trackvector;

/*=============================================================================
 class definitions
=============================================================================*/

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
  void set_time_start(const size_t time_start);
  void set_parent_id(const size_t parent_id);
  void set_child_ids(const size_t left, const size_t right);
  size_t get_id() const;
  size_t get_time_start() const;
  size_t get_parent_id() const;
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
//// class Tracking
////
class Tracking{
 public:
  Tracking();
  Tracking(const HypothesesGraph& graph, const size_t index);

  Trackvector tracks_;
  size_t index_;
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & archive, const unsigned int version) {
    (void) version;
    archive & tracks_;
    archive & index_;
  }
}; // class Tracking

/*=============================================================================
  functor definitions
=============================================================================*/

////
//// class TrackValue
////
template<typename TrackFeatureExtractor_T, typename FeatureAggregator_T>
class TrackValue {
 public:
  TrackValue() : trackfeatureextractor_(), featureaggregator_() {};
  virtual feature_type operator()(const Track& track) const;

  TrackFeatureExtractor_T trackfeatureextractor_;
  FeatureAggregator_T featureaggregator_;
}; // class TrackValue

} // namespace pgmlink

#endif // PGMLINK_HIGHER_ORDER_FEATURES_H
