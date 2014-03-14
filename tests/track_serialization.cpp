#define BOOST_TEST_MODULE outlier_detection_test

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "pgmlink/tracks.h"

using namespace pgmlink;
using namespace boost;

BOOST_AUTO_TEST_CASE( Track_serialization ) {
  // Set up the test data
  feature_type x_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {4., 5.}
  };
  Track track;
  for(size_t i = 0; i < 7; i++) {
    feature_array y(x_array[i], x_array[i]+2);
    Traxel traxel;
    traxel.features["feature_x"] = y;
    track.traxels_.push_back(traxel);
  }
  track.set_id(5);
  track.set_time_start(7);
  
  {
    // Create an output archive
    std::ofstream ofilestream("track_serialization_test");
  
    // Save the track
    archive::text_oarchive oarchive(ofilestream);
    oarchive << track;
  }
  Track loaded_track;
  {
    std::ifstream ifilestream("track_serialization_test");
    // Load the track
    
    archive::text_iarchive iarchive(ifilestream);
    iarchive >> loaded_track;
  }

  BOOST_CHECK_EQUAL(track.get_length(), loaded_track.get_length());
  BOOST_CHECK_EQUAL(track.get_id(), loaded_track.get_id());
  BOOST_CHECK_EQUAL(track.get_time_start(), loaded_track.get_time_start());

  BOOST_CHECK_EQUAL(track.get_id(), 5);
  BOOST_CHECK_EQUAL(track.get_time_start(), 7);
  
  for(size_t i=0; i < track.get_length(); i++) {
    bool traxel_equal = track.traxels_[i] == loaded_track.traxels_[i];
    BOOST_CHECK(traxel_equal);
  }

}

/*
BOOST_AUTO_TEST_CASE( load_Tracking ) {
  std::vector<std::vector<Event> > events;
  std::vector<Event> events_t0;
  Event e00; e00.type = Event::Appearance;
  e00.traxel_ids.push_back(3);
  events_t0.push_back(e00);
  events.push_back(events_t0);
  
  std::vector<Event> events_t1;
  Event e10; e10.type = Event::Move;
  e10.traxel_ids.push_back(3);
  e10.traxel_ids.push_back(2);
  events_t1.push_back(e10);
  events.push_back(events_t1);
  
  std::vector<Event> events_t2;
  Event e20; e20.type = Event::Division;
  e20.traxel_ids.push_back(2);
  e20.traxel_ids.push_back(2);
  e20.traxel_ids.push_back(5);
  events_t1.push_back(e10);
  events.push_back(events_t1);

  
  TraxelStore ts;
  Tracking tracking(events, ts);
  Trackvector::iterator t_it = tracking.tracks_.begin();
  for(; t_it != tracking.tracks_.end(); t_it++) {
    std::cout << t_it->get_id() << "\t";
    std::cout << t_it->get_length() << "\t";
    std::cout << t_it->get_time_start() << "\t";
    std::cout << t_it->parent_id_ << "\t";
    std::cout << (t_it->child_ids_)[0] << "\t";
    std::cout << (t_it->child_ids_)[1] << std::endl;
  }
  
}
*/
