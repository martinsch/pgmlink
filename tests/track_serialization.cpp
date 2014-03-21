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

BOOST_AUTO_TEST_CASE( load_Tracking_from_hypotheses_graph ) {
  HypothesesGraph graph;
  HypothesesGraph::Node n00 = graph.add_node(1);
  HypothesesGraph::Node n01 = graph.add_node(1);
  HypothesesGraph::Node n10 = graph.add_node(2);
  HypothesesGraph::Node n11 = graph.add_node(2);
  HypothesesGraph::Node n20 = graph.add_node(3);
  HypothesesGraph::Node n21 = graph.add_node(3);

  HypothesesGraph::Arc a13 = graph.addArc(n00, n10);
  HypothesesGraph::Arc a35 = graph.addArc(n10, n20);
  HypothesesGraph::Arc a36 = graph.addArc(n10, n21);
  HypothesesGraph::Arc a24 = graph.addArc(n01, n11);

  // corresponding traxels
  Traxel t00, t01, t10, t11, t20, t21;
  t00.Id = 1;
  t00.Timestep = 1;
  t00.features["com"].push_back(0.);
  t00.features["com"].push_back(0.);
  t01.Id = 2;
  t01.Timestep = 1;
  t01.features["com"].push_back(0.);
  t01.features["com"].push_back(1.);
  t10.Id = 3;
  t10.Timestep = 2;
  t10.features["com"].push_back(1.);
  t10.features["com"].push_back(0.);
  t11.Id = 4;
  t11.Timestep = 2;
  t11.features["com"].push_back(1.);
  t11.features["com"].push_back(1.);
  t20.Id = 5;
  t20.Timestep = 3;
  t20.features["com"].push_back(2.);
  t20.features["com"].push_back(0.);
  t21.Id = 6;
  t21.Timestep = 3;
  t21.features["com"].push_back(2.);
  t21.features["com"].push_back(1.);

  graph.add(node_traxel());
  graph.get(node_traxel()).set(n00, t00);
  graph.get(node_traxel()).set(n01, t01);
  graph.get(node_traxel()).set(n10, t10);
  graph.get(node_traxel()).set(n11, t11);
  graph.get(node_traxel()).set(n20, t20);
  graph.get(node_traxel()).set(n21, t21);

  graph.add(node_active_count());
  graph.get(node_active_count()).set(n00, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n01, std::vector<size_t>(1,0));
  graph.get(node_active_count()).set(n10, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n11, std::vector<size_t>(1,0));
  graph.get(node_active_count()).set(n20, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n21, std::vector<size_t>(1,1));

  graph.add(arc_active_count());
  graph.get(arc_active_count()).set(a13, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a35, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a36, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a24, std::vector<bool>(1,false));

  graph.add(division_active_count());
  graph.get(division_active_count()).set(n00, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n01, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n10, std::vector<bool>(1,true));
  graph.get(division_active_count()).set(n11, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n20, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n21, std::vector<bool>(1,false));
  
  Tracking tracking(graph, 0);
  Trackvector tracks = tracking.tracks_;
  BOOST_CHECK_EQUAL(tracks.size(), 3);
  size_t traxel_count = 
    tracks[0].get_length() + tracks[1].get_length() + tracks[2].get_length();
  BOOST_CHECK_EQUAL(traxel_count, 4);
 
}

BOOST_AUTO_TEST_CASE( load_Tracking_from_tracklet_hypotheses_graph ) {
  HypothesesGraph graph;
  HypothesesGraph::Node n00 = graph.add_node(1);
  HypothesesGraph::Node n01 = graph.add_node(1);
  HypothesesGraph::Node n10 = graph.add_node(3);
  HypothesesGraph::Node n11 = graph.add_node(3);

  HypothesesGraph::Arc a13 = graph.addArc(n00, n10);
  HypothesesGraph::Arc a14 = graph.addArc(n00, n11);

  // corresponding traxels
  Traxel t00, t01, t10, t11, t20, t21;
  t00.Id = 1;
  t00.Timestep = 1;
  t00.features["com"].push_back(0.);
  t00.features["com"].push_back(0.);
  t01.Id = 2;
  t01.Timestep = 1;
  t01.features["com"].push_back(0.);
  t01.features["com"].push_back(1.);
  t10.Id = 3;
  t10.Timestep = 2;
  t10.features["com"].push_back(1.);
  t10.features["com"].push_back(0.);
  t11.Id = 4;
  t11.Timestep = 2;
  t11.features["com"].push_back(1.);
  t11.features["com"].push_back(1.);
  t20.Id = 5;
  t20.Timestep = 3;
  t20.features["com"].push_back(2.);
  t20.features["com"].push_back(0.);
  t21.Id = 6;
  t21.Timestep = 3;
  t21.features["com"].push_back(2.);
  t21.features["com"].push_back(1.);

  graph.add(node_tracklet());
  Traxelvector tlet1;
  Traxelvector tlet2;
  Traxelvector tlet3;
  Traxelvector tlet4;
  tlet1.push_back(t00); tlet1.push_back(t10);
  tlet2.push_back(t01); tlet2.push_back(t11);
  tlet3.push_back(t20);
  tlet4.push_back(t21);
  graph.get(node_tracklet()).set(n00, tlet1);
  graph.get(node_tracklet()).set(n01, tlet2);
  graph.get(node_tracklet()).set(n10, tlet3);
  graph.get(node_tracklet()).set(n11, tlet4);

  graph.add(node_active_count());
  graph.get(node_active_count()).set(n00, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n01, std::vector<size_t>(1,0));
  graph.get(node_active_count()).set(n10, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n11, std::vector<size_t>(1,1));

  graph.add(arc_active_count());
  graph.get(arc_active_count()).set(a13, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a14, std::vector<bool>(1,true));

  graph.add(division_active_count());
  graph.get(division_active_count()).set(n00, std::vector<bool>(1,true));
  graph.get(division_active_count()).set(n01, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n10, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n11, std::vector<bool>(1,false));
  
  Tracking tracking(graph, 0);
  Trackvector tracks = tracking.tracks_;
  BOOST_CHECK_EQUAL(tracks.size(), 3);
  size_t traxel_count = 
    tracks[0].get_length() + tracks[1].get_length() + tracks[2].get_length();
  BOOST_CHECK_EQUAL(traxel_count, 4);
}

BOOST_AUTO_TEST_CASE( load_Tracking_from_hypotheses_graph_no_div_active ) {
  HypothesesGraph graph;
  HypothesesGraph::Node n00 = graph.add_node(1);
  HypothesesGraph::Node n01 = graph.add_node(1);
  HypothesesGraph::Node n10 = graph.add_node(2);
  HypothesesGraph::Node n11 = graph.add_node(2);
  HypothesesGraph::Node n20 = graph.add_node(3);
  HypothesesGraph::Node n21 = graph.add_node(3);

  HypothesesGraph::Arc a13 = graph.addArc(n00, n10);
  HypothesesGraph::Arc a35 = graph.addArc(n10, n20);
  HypothesesGraph::Arc a36 = graph.addArc(n10, n21);
  HypothesesGraph::Arc a24 = graph.addArc(n01, n11);

  // corresponding traxels
  Traxel t00, t01, t10, t11, t20, t21;
  t00.Id = 1;
  t00.Timestep = 1;
  t00.features["com"].push_back(0.);
  t00.features["com"].push_back(0.);
  t01.Id = 2;
  t01.Timestep = 1;
  t01.features["com"].push_back(0.);
  t01.features["com"].push_back(1.);
  t10.Id = 3;
  t10.Timestep = 2;
  t10.features["com"].push_back(1.);
  t10.features["com"].push_back(0.);
  t11.Id = 4;
  t11.Timestep = 2;
  t11.features["com"].push_back(1.);
  t11.features["com"].push_back(1.);
  t20.Id = 5;
  t20.Timestep = 3;
  t20.features["com"].push_back(2.);
  t20.features["com"].push_back(0.);
  t21.Id = 6;
  t21.Timestep = 3;
  t21.features["com"].push_back(2.);
  t21.features["com"].push_back(1.);

  graph.add(node_traxel());
  graph.get(node_traxel()).set(n00, t00);
  graph.get(node_traxel()).set(n01, t01);
  graph.get(node_traxel()).set(n10, t10);
  graph.get(node_traxel()).set(n11, t11);
  graph.get(node_traxel()).set(n20, t20);
  graph.get(node_traxel()).set(n21, t21);

  graph.add(node_active_count());
  graph.get(node_active_count()).set(n00, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n01, std::vector<size_t>(1,0));
  graph.get(node_active_count()).set(n10, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n11, std::vector<size_t>(1,0));
  graph.get(node_active_count()).set(n20, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n21, std::vector<size_t>(1,1));

  graph.add(arc_active_count());
  graph.get(arc_active_count()).set(a13, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a35, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a36, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a24, std::vector<bool>(1,false));
 
  Tracking tracking(graph, 0);
  Trackvector tracks = tracking.tracks_;
  BOOST_CHECK_EQUAL(tracks.size(), 3);
  size_t traxel_count = 
    tracks[0].get_length() + tracks[1].get_length() + tracks[2].get_length();
  BOOST_CHECK_EQUAL(traxel_count, 4);
 
}

BOOST_AUTO_TEST_CASE( load_Tracking_from_tracklet_hypotheses_graph_no_div_active ) {
  HypothesesGraph graph;
  HypothesesGraph::Node n00 = graph.add_node(1);
  HypothesesGraph::Node n01 = graph.add_node(1);
  HypothesesGraph::Node n10 = graph.add_node(3);
  HypothesesGraph::Node n11 = graph.add_node(3);

  HypothesesGraph::Arc a13 = graph.addArc(n00, n10);
  HypothesesGraph::Arc a14 = graph.addArc(n00, n11);

  // corresponding traxels
  Traxel t00, t01, t10, t11, t20, t21;
  t00.Id = 1;
  t00.Timestep = 1;
  t00.features["com"].push_back(0.);
  t00.features["com"].push_back(0.);
  t01.Id = 2;
  t01.Timestep = 1;
  t01.features["com"].push_back(0.);
  t01.features["com"].push_back(1.);
  t10.Id = 3;
  t10.Timestep = 2;
  t10.features["com"].push_back(1.);
  t10.features["com"].push_back(0.);
  t11.Id = 4;
  t11.Timestep = 2;
  t11.features["com"].push_back(1.);
  t11.features["com"].push_back(1.);
  t20.Id = 5;
  t20.Timestep = 3;
  t20.features["com"].push_back(2.);
  t20.features["com"].push_back(0.);
  t21.Id = 6;
  t21.Timestep = 3;
  t21.features["com"].push_back(2.);
  t21.features["com"].push_back(1.);

  graph.add(node_tracklet());
  Traxelvector tlet1;
  Traxelvector tlet2;
  Traxelvector tlet3;
  Traxelvector tlet4;
  tlet1.push_back(t00); tlet1.push_back(t10);
  tlet2.push_back(t01); tlet2.push_back(t11);
  tlet3.push_back(t20);
  tlet4.push_back(t21);
  graph.get(node_tracklet()).set(n00, tlet1);
  graph.get(node_tracklet()).set(n01, tlet2);
  graph.get(node_tracklet()).set(n10, tlet3);
  graph.get(node_tracklet()).set(n11, tlet4);

  graph.add(node_active_count());
  graph.get(node_active_count()).set(n00, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n01, std::vector<size_t>(1,0));
  graph.get(node_active_count()).set(n10, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n11, std::vector<size_t>(1,1));

  graph.add(arc_active_count());
  graph.get(arc_active_count()).set(a13, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a14, std::vector<bool>(1,true));
  
  Tracking tracking(graph, 0);
  Trackvector tracks = tracking.tracks_;
  BOOST_CHECK_EQUAL(tracks.size(), 3);
  size_t traxel_count = 
    tracks[0].get_length() + tracks[1].get_length() + tracks[2].get_length();
  BOOST_CHECK_EQUAL(traxel_count, 4);
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
