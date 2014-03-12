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
#include "pgmlink/classifier_auxiliary.h"
#include "pgmlink/track_features.h"

// boost
#include <boost/shared_ptr.hpp>

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

////
//// class Track
////
class Track {
 public:
  Track();
  Track(size_t id, size_t time_start, const Traxelvector& traxels);

  void set_id(const size_t id);
  size_t get_id() const;

  void set_time_start(const size_t time_start);
  size_t get_time_start() const;
  
  size_t get_length() const;

  size_t id_;
  size_t time_start_;
  Traxelvector traxels_;
}; // class Track

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

} // namespace pgmlink

#endif // PGMLINK_TRACKS_H
