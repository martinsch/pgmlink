#define BOOST_TEST_MODULE reasoner_constracking_test

#include <vector>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/bind.hpp>

#include <lemon/color.h>
#include <lemon/graph_to_eps.h>

#include <lemon/core.h>
#include <lemon/concepts/digraph.h>
#include <lemon/list_graph.h>
#include <lemon/maps.h>

#include "pgmlink/graph.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/energy.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/traxels.h"
#include "pgmlink/track.h"

using namespace Tracking;
using namespace std;
using namespace boost;

BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_Merger ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	typedef HypothesesGraph::ArcIt ArcIt2;
	typedef HypothesesGraph::Arc Arc;
	typedef HypothesesGraph::NodeIt NodeIt;
	typedef HypothesesGraph::Node Node;
	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2      3       4
	//  o                       o
	//    |                    |
	//      ---- o ---- o ----
	//    |                    |
	//  o                       o
	TraxelStore ts;
	Traxel n11, n12, n21, n31, n41, n42;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	n11.Id = 1; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
	n11.features["com"] = com; n11.features["divProb"] = divProb;
	add(ts,n11);
	n12.Id = 3; n12.Timestep = 1; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1;
	n12.features["com"] = com; n12.features["divProb"] = divProb;
	add(ts,n12);
	n21.Id = 10; n21.Timestep = 2; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.1;
	n21.features["com"] = com; n21.features["divProb"] = divProb;
	add(ts,n21);
	n31.Id = 11; n31.Timestep = 3; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.1;
	n31.features["com"] = com; n31.features["divProb"] = divProb;
	add(ts,n31);
	n41.Id = 12; n41.Timestep = 4; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
	n41.features["com"] = com; n41.features["divProb"] = divProb;
	add(ts,n41);
	n42.Id = 13; n42.Timestep = 4; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
	n42.features["com"] = com; n42.features["divProb"] = divProb;
	add(ts,n42);


	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

	ConsTracking tracking = ConsTracking(
			  2, // max_number_objects
		      20, // max_neighbor_distance
			  0.3, // division_threshold
			  "none", // random_forest_filename
	  	      false, // cellness_by_random_forest
	  	      0, // forbidden_cost
	  	      true, // with_constraints
	  	      false, // fixed_detections
	  	      0.0 // ep_gap
	  	      );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts);

	size_t count_moves = 0;
	size_t t = 1;
	for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin(); it_t != events.end(); ++it_t) {
		// events:
		// t = 1: 2x move
		// t = 2: 1x move, 1x merger
		// t = 3: 2x move, 1x merger
		if (t==1 || t ==2) {
			BOOST_CHECK_EQUAL(it_t->size(),2);
		} else { // t == 3
			BOOST_CHECK_EQUAL(it_t->size(),3);
		}

		for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {
			Event e = *it;
			BOOST_CHECK_NE(e.type, Event::Disappearance);
			BOOST_CHECK_NE(e.type, Event::Appearance);
			BOOST_CHECK_NE(e.type, Event::Division);


			if (e.type == Event::Move) {
				++count_moves;
			}
		}
		++t;
	}
	BOOST_CHECK_EQUAL(count_moves, 5);
}


BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_Division ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	typedef HypothesesGraph::ArcIt ArcIt2;
	typedef HypothesesGraph::Arc Arc;
	typedef HypothesesGraph::NodeIt NodeIt;
	typedef HypothesesGraph::Node Node;
	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2      3
	//                  o
	//                |
	//  D ---- D ----
	//                |
	//                  o
	TraxelStore ts;
	Traxel n11, n21, n31, n32;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	n11.Id = 1; n11.Timestep = 1; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.8;
	n11.features["com"] = com; n11.features["divProb"] = divProb;
	add(ts,n11);
	n21.Id = 10; n21.Timestep = 2; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.7;
	n21.features["com"] = com; n21.features["divProb"] = divProb;
	add(ts,n21);
	n31.Id = 11; n31.Timestep = 3; com[0] = 1; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
	n31.features["com"] = com; n31.features["divProb"] = divProb;
	add(ts,n31);
	n32.Id = 13; n32.Timestep = 3; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1;
	n32.features["com"] = com; n32.features["divProb"] = divProb;
	add(ts,n32);


	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

	ConsTracking tracking = ConsTracking(
			  2, // max_number_objects
		      20, // max_neighbor_distance
			  0.3, // division_threshold
			  "none", // random_forest_filename
	  	      false, // cellness_by_random_forest
	  	      0, // forbidden_cost
	  	      true, // with_constraints
	  	      false, // fixed_detections
	  	      0.01 // ep_gap
	  	      );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts);

	size_t count_moves = 0;
	size_t count_divisions = 0;
	size_t t = 1;
	for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin(); it_t != events.end(); ++it_t) {
		// events:
		// t = 1: 1x move
		// t = 2: 1x division
		BOOST_CHECK_EQUAL(it_t->size(),1);

		for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {
			Event e = *it;
			BOOST_CHECK_NE(e.type, Event::Disappearance);
			BOOST_CHECK_NE(e.type, Event::Appearance);

			if (e.type == Event::Move) {
				++count_moves;
			} else if (e.type == Event::Division) {
				++count_divisions;
			}
		}
		++t;
	}
	BOOST_CHECK_EQUAL(count_moves, 1);
	BOOST_CHECK_EQUAL(count_divisions, 1);
}


BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_SimpleMove ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	typedef HypothesesGraph::ArcIt ArcIt2;
	typedef HypothesesGraph::Arc Arc;
	typedef HypothesesGraph::NodeIt NodeIt;
	typedef HypothesesGraph::Node Node;
	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2
	//  o ------ o
	TraxelStore ts;
	Traxel n11, n21;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	n11.Id = 1; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
	n11.features["com"] = com; n11.features["divProb"] = divProb;
	add(ts,n11);
	n21.Id = 10; n21.Timestep = 2; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
	n21.features["com"] = com; n21.features["divProb"] = divProb;
	add(ts,n21);

	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

	ConsTracking tracking = ConsTracking(
			  2, // max_number_objects
		      20, // max_neighbor_distance
			  0.3, // division_threshold
			  "none", // random_forest_filename
	  	      false, // cellness_by_random_forest
	  	      0, // forbidden_cost
	  	      true, // with_constraints
	  	      false, // fixed_detections
	  	      0.0 // ep_gap
	  	      );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts);

	BOOST_CHECK_EQUAL(events.size(),1);
	size_t count_moves = 0;
	for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin(); it_t != events.end(); ++it_t) {
		for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {
			Event e = *it;
			if (e.type == Event::Move) {
				++count_moves;
			}
		}
	}
	BOOST_CHECK_EQUAL(count_moves, 1);
}


BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_Merger_Volume ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	typedef HypothesesGraph::ArcIt ArcIt2;
	typedef HypothesesGraph::Arc Arc;
	typedef HypothesesGraph::NodeIt NodeIt;
	typedef HypothesesGraph::Node Node;
	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2      3       4
	//  o                       o
	//    |                    |
	//      ---- o ---- o ----
	//    |                    |
	//  o                       o
	TraxelStore ts;
	Traxel n11, n12, n21, n31, n41, n42;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	feature_array count(feature_array::difference_type(1));
	n11.Id = 1; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 2;
	n11.features["com"] = com; n11.features["divProb"] = divProb; n11.features["count"] = count;
	add(ts,n11);
	n12.Id = 3; n12.Timestep = 1; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1; count[0] = 100;
	n12.features["com"] = com; n12.features["divProb"] = divProb; n12.features["count"] = count;
	add(ts,n12);
	n21.Id = 10; n21.Timestep = 2; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.1; count[0] = 90;
	n21.features["com"] = com; n21.features["divProb"] = divProb; n21.features["count"] = count;
	add(ts,n21);
	n31.Id = 11; n31.Timestep = 3; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.1; count[0] = 90;
	n31.features["com"] = com; n31.features["divProb"] = divProb; n31.features["count"] = count;
	add(ts,n31);
	n41.Id = 12; n41.Timestep = 4; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 110;
	n41.features["com"] = com; n41.features["divProb"] = divProb; n41.features["count"] = count;
	add(ts,n41);
	n42.Id = 13; n42.Timestep = 4; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 60;
	n42.features["com"] = com; n42.features["divProb"] = divProb; n42.features["count"] = count;
	add(ts,n42);


	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

	ConsTracking tracking = ConsTracking(
			  2, // max_number_objects
		      20, // max_neighbor_distance
			  0.3, // division_threshold
			  "none", // random_forest_filename
	  	      true, // detectino_by_volume
	  	      0, // forbidden_cost
	  	      true, // with_constraints
	  	      false, // fixed_detections
	  	      0.0, // ep_gap
	  	      double(99) // avg_obj_size
	  	      );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts);

	size_t count_moves = 0;
	size_t t = 1;
	for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin(); it_t != events.end(); ++it_t) {
		// events:
		// t = 1: 1x move
		// t = 2: 1x move
		// t = 3: 1x move
		BOOST_CHECK_EQUAL(it_t->size(),1);

		for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {
			Event e = *it;
			BOOST_CHECK_NE(e.type, Event::Disappearance);
			BOOST_CHECK_NE(e.type, Event::Appearance);
			BOOST_CHECK_NE(e.type, Event::Division);
			BOOST_CHECK_NE(e.type, Event::Merger);


			if (e.type == Event::Move) {
				++count_moves;
			}
		}
		++t;
	}
	BOOST_CHECK_EQUAL(count_moves, 3);
}


BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_Disappearance ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	typedef HypothesesGraph::ArcIt ArcIt2;
	typedef HypothesesGraph::Arc Arc;
	typedef HypothesesGraph::NodeIt NodeIt;
	typedef HypothesesGraph::Node Node;
	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2
	//  1
	//    |
	//    |
	//    |
	//      ---- 1
	//    |
	//  1
	TraxelStore ts;
	Traxel n11, n12, n21;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	feature_array count(feature_array::difference_type(1));
	n11.Id = 1; n11.Timestep = 1; com[0] = 18; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n11.features["com"] = com; n11.features["divProb"] = divProb; n11.features["count"] = count;
	add(ts,n11);
	n12.Id = 3; n12.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n12.features["com"] = com; n12.features["divProb"] = divProb; n12.features["count"] = count;
	add(ts,n12);
	n21.Id = 10; n21.Timestep = 2; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n21.features["com"] = com; n21.features["divProb"] = divProb; n21.features["count"] = count;
	add(ts,n21);

	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

	ConsTracking tracking = ConsTracking(
			  2, // max_number_objects
		      20, // max_neighbor_distance
			  0.3, // division_threshold
			  "none", // random_forest_filename
	  	      true, // detection_by_volume
	  	      0, // forbidden_cost
	  	      true, // with_constraints
	  	      false, // fixed_detections
	  	      0.0, // ep_gap
	  	      double(1), // avg_obj_size
	  	      false, // with_appearance
	  	      true   // with_disappearance
	  	      );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts);

	size_t num_events = 0;
	size_t t = 1;
	BOOST_CHECK_EQUAL(events.size(),1); // 1 timesteps
	for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin(); it_t != events.end(); ++it_t) {
		for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {
			Event e = *it;
			++num_events;
			if (e.type == Event::Move) {
				BOOST_CHECK_EQUAL(e.traxel_ids[0],3);
				BOOST_CHECK_EQUAL(e.traxel_ids[1],10);
			} else if (e.type == Event::Disappearance) {
				BOOST_CHECK_EQUAL(e.traxel_ids[0],1);
			} else {
				BOOST_CHECK(false);
			}
		}
		++t;
	}
	BOOST_CHECK_EQUAL(num_events,2); // 2 events
}


BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_AppearanceAndDisappearance ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	typedef HypothesesGraph::ArcIt ArcIt2;
	typedef HypothesesGraph::Arc Arc;
	typedef HypothesesGraph::NodeIt NodeIt;
	typedef HypothesesGraph::Node Node;
	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2        3
	//
	//  1                 1
	//    |             |
	//    |             |
	//    |             |
	//      ---- 1 ----
	//    |             |
	//  1                 1
	TraxelStore ts;
	Traxel n11, n12, n21, n31, n32;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	feature_array count(feature_array::difference_type(1));
	n11.Id = 1; n11.Timestep = 1; com[0] = 18; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n11.features["com"] = com; n11.features["divProb"] = divProb; n11.features["count"] = count;
	add(ts,n11);
	n12.Id = 3; n12.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n12.features["com"] = com; n12.features["divProb"] = divProb; n12.features["count"] = count;
	add(ts,n12);
	n21.Id = 10; n21.Timestep = 2; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n21.features["com"] = com; n21.features["divProb"] = divProb; n21.features["count"] = count;
	add(ts,n21);
	n31.Id = 11; n31.Timestep = 3; com[0] = 18; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n31.features["com"] = com; n31.features["divProb"] = divProb; n31.features["count"] = count;
	add(ts,n31);
	n32.Id = 15; n32.Timestep = 3; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n32.features["com"] = com; n32.features["divProb"] = divProb; n32.features["count"] = count;
	add(ts,n32);

	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

	ConsTracking tracking = ConsTracking(
			  2, // max_number_objects
		      20, // max_neighbor_distance
			  0.3, // division_threshold
			  "none", // random_forest_filename
	  	      true, // detection_by_volume
	  	      0, // forbidden_cost
	  	      true, // with_constraints
	  	      false, // fixed_detections
	  	      0.0, // ep_gap
	  	      double(1), // avg_obj_size
	  	      true, // with_appearance
	  	      true   // with_disappearance
	  	      );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts);

	size_t t = 1;
	size_t num_events = 0;
	for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin(); it_t != events.end(); ++it_t) {
		for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {
			Event e = *it;
			++num_events;
			if (e.type == Event::Move) {
				if (t==1) {
					BOOST_CHECK_EQUAL(e.traxel_ids[0],3);
					BOOST_CHECK_EQUAL(e.traxel_ids[1],10);
				} else if (t==2) {
					BOOST_CHECK_EQUAL(e.traxel_ids[0],10);
					BOOST_CHECK_EQUAL(e.traxel_ids[1],15);
				} else {
					BOOST_CHECK(false);
				}
			} else if (e.type == Event::Disappearance) {
				BOOST_CHECK_EQUAL(e.traxel_ids[0],1);
			} else if (e.type == Event::Appearance) {
				BOOST_CHECK_EQUAL(e.traxel_ids[0],11);
			} else {
				BOOST_CHECK(false);
			}
		}
		++t;
	}
	BOOST_CHECK_EQUAL(num_events, 4);
}



BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_FixDetections ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	typedef HypothesesGraph::ArcIt ArcIt2;
	typedef HypothesesGraph::Arc Arc;
	typedef HypothesesGraph::NodeIt NodeIt;
	typedef HypothesesGraph::Node Node;
	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2
	//  1
	//    |
	//    |
	//      ---- 1
	//    |
	//  2
	TraxelStore ts;
	Traxel n11, n12, n21;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	feature_array count(feature_array::difference_type(1));
	n11.Id = 1; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n11.features["com"] = com; n11.features["divProb"] = divProb; n11.features["count"] = count;
	add(ts,n11);
	n12.Id = 3; n12.Timestep = 1; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1; count[0] = 10;
	n12.features["com"] = com; n12.features["divProb"] = divProb; n12.features["count"] = count;
	add(ts,n12);
	n21.Id = 10; n21.Timestep = 2; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1; count[0] = 1;
	n21.features["com"] = com; n21.features["divProb"] = divProb; n21.features["count"] = count;
	add(ts,n21);

	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

	ConsTracking tracking = ConsTracking(
			  2, // max_number_objects
		      20, // max_neighbor_distance
			  0.3, // division_threshold
			  "none", // random_forest_filename
	  	      true, // detection_by_volume
	  	      0, // forbidden_cost
	  	      true, // with_constraints
	  	      true, // fixed_detections
	  	      0.0, // ep_gap
	  	      double(1), // avg_obj_size
			  false, // with_appearance
			  true   // with_disappearance
	  	      );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts);

	size_t t = 1;
	for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin(); it_t != events.end(); ++it_t) {
		for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {
			Event e = *it;
			if (e.type == Event::Move) {
				BOOST_CHECK_EQUAL(e.traxel_ids[0],1);
				BOOST_CHECK_EQUAL(e.traxel_ids[1],10);
			}
			if (e.type == Event::Disappearance) {
				BOOST_CHECK_EQUAL(e.traxel_ids[0],3);
			}
			if (e.type == Event::Merger) {
				BOOST_CHECK_EQUAL(e.traxel_ids[0],3);
				BOOST_CHECK_EQUAL(e.traxel_ids[1],2);
			}
		}
		++t;
	}
}
