#define BOOST_TEST_MODULE reasoner_constracking_test

#include <vector>
#include <iostream>
#include <set>

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
#include "pgmlink/feature.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/traxels.h"
#include "pgmlink/tracking.h"
#include "pgmlink/field_of_view.h"

using namespace pgmlink;
using namespace std;
using namespace boost;

BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_Merger ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

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

	FieldOfView fov(0, 0, 0, 0, 4, 5, 5, 5); // tlow, xlow, ylow, zlow, tup, xup, yup, zup
	ConsTracking tracking = ConsTracking(2, // max_number_objects
					     false, // detection_by_volume
					     double(1.1), // avg_obj_size
					     20, // max_neighbor_distance
					     true, //with_divisions
					     0.3, // division_threshold
					     "none", // random_forest_filename
					     fov
		  	      );


	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.0, // ep_gap
							    false, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    1500., // disappearance_cost,
							    1500., // appearance_cost
							    false, //with_merger_resolution
							    3, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );

	size_t count_moves = 0;
    size_t t = 1;
    for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin()+1; it_t != events.end(); ++it_t) {        
		// events:
		// t = 1: 2x move
		// t = 2: 1x move, 1x merger
		// t = 3: 2x move, 1x merger
		if (t==2 || t ==3) {
			BOOST_CHECK_EQUAL(it_t->size(),2);
		} else { // t == 1
			BOOST_CHECK_EQUAL(it_t->size(),3);
		}

        for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {            
			Event e = *it;
            cout << t << ": " << e << endl;
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

	FieldOfView fov(0, 0, 0, 0, 4, 5, 5, 5); // tlow, xlow, ylow, zlow, tup, xup, yup, zup
	ConsTracking tracking = ConsTracking(
					     2, // max_number_objects
					     false, // detection_by_volume
					     double(1.1), // avg_obj_size
					     20, // max_neighbor_distance
					     true, //with_divisions
					     0.3, // division_threshold
					     "none", // random_forest_filename
					     fov
					     );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,

							    0, // forbidden_cost
							    0.0, // ep_gap
							    false, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    1500., // disappearance_cost,
							    1500., // appearance_cost
							    false, //with_merger_resolution
							    3, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );

	size_t count_moves = 0;
	size_t count_divisions = 0;
	size_t t = 1;
    for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin()+1; it_t != events.end(); ++it_t) {
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

	FieldOfView fov(0, 0, 0, 0, 4, 5, 5, 5); // tlow, xlow, ylow, zlow, tup, xup, yup, zup
	ConsTracking tracking = ConsTracking(
					     2, // max_number_objects
					     false, // detection_by_volume
					     double(1.1), // avg_obj_size
					     20, // max_neighbor_distance
					     true, //with_divisions
					     0.3, // division_threshold
					     "none", // random_forest_filename
					     fov
				  );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.0, // ep_gap
							    false, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    1500., // disappearance_cost,
							    1500., // appearance_cost
							    false, //with_merger_resolution
							    3, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );

    BOOST_CHECK_EQUAL(events.size(),2);
	size_t count_moves = 0;
    for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin()+1; it_t != events.end(); ++it_t) {
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

	FieldOfView fov(0, 0, 0, 0, 4, 5, 5, 5); // tlow, xlow, ylow, zlow, tup, xup, yup, zup
	ConsTracking tracking = ConsTracking(
					     2, // max_number_objects
					     true, // detection_by_volume
					     double(99), // avg_obj_size
					     20, // max_neighbor_distance
					     true, //with_divisions
					     0.3, // division_threshold
					     "none", // random_forest_filename
					     fov
					     );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.0, // ep_gap
							    false, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    1500., // disappearance_cost,
							    1500., // appearance_cost
							    false, //with_merger_resolution
							    3, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
);

	size_t count_moves = 0;
	size_t t = 1;
    for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin()+1; it_t != events.end(); ++it_t) {
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

	FieldOfView fov(0, 0, 0, 0, 4, 5, 5, 5); // tlow, xlow, ylow, zlow, tup, xup, yup, zup
	ConsTracking tracking = ConsTracking(
					     2, // max_number_objects
					     true, // detection_by_volume
					     double(1), // avg_obj_size
					     20, // max_neighbor_distance
					     true, //with_divisions
					     0.3, // division_threshold
					     "none", // random_forest_filename
					     fov
					     );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.0, // ep_gap
							    false, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    1500., // disappearance_cost,
							    1500., // appearance_cost
							    false, //with_merger_resolution
							    3, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );

	size_t num_events = 0;
	size_t t = 1;
    BOOST_CHECK_EQUAL(events.size(),2); // 2 timesteps
    for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin()+1; it_t != events.end(); ++it_t) {
		for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {
			Event e = *it;
			++num_events;
			if (e.type == Event::Move) {
				BOOST_CHECK_EQUAL(e.traxel_ids[0],3);
				BOOST_CHECK_EQUAL(e.traxel_ids[1],10);
//			} else if (e.type == Event::Disappearance) {
//				BOOST_CHECK_EQUAL(e.traxel_ids[0],1);
			} else {
				BOOST_CHECK(false);
			}
		}
		++t;
	}
	BOOST_CHECK_EQUAL(num_events,1); // 1 events
}


BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_AppearanceAndDisappearance ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

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

	FieldOfView fov(0, 0, 0, 0, 4, 5, 5, 5); // tlow, xlow, ylow, zlow, tup, xup, yup, zup
	ConsTracking tracking = ConsTracking(
					     2, // max_number_object
					     true, // detection_by_volume
					     double(1.), // avg_obj_size
					     20, // max_neighbor_distance
					     true, //with_divisions
					     0.3, // division_threshold
					     "none", // random_forest_filename
					     fov
					     );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.0, // ep_gap
							    false, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    1500., // disappearance_cost,
							    1500., // appearance_cost
							    false, //with_merger_resolution
							    3, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );

	size_t t = 1;
	size_t num_events = 0;
    for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin()+1; it_t != events.end(); ++it_t) {
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
//			} else if (e.type == Event::Disappearance) {
//				BOOST_CHECK_EQUAL(e.traxel_ids[0],1);
//			} else if (e.type == Event::Appearance) {
//				BOOST_CHECK_EQUAL(e.traxel_ids[0],11);
			} else {
				BOOST_CHECK(false);
			}
		}
		++t;
	}
	BOOST_CHECK_EQUAL(num_events, 2);
}



BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_Appearance ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2        3
	//
	//  1 ------ 1 ------ 1
	//  1 ------ 1 ------ 1
	//  2 ------ 2 ------ 2
	//           1 ------ 1
	//                    1
	TraxelStore ts;
	Traxel n11, n12, n13, n21, n22, n23, n24, n31, n32, n33, n34, n35;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	feature_array count(feature_array::difference_type(1));
	n11.Id = 1; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n11.features["com"] = com; n11.features["divProb"] = divProb; n11.features["count"] = count;
	add(ts,n11);
	n21.Id = 2; n21.Timestep = 2; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n21.features["com"] = com; n21.features["divProb"] = divProb; n21.features["count"] = count;
	add(ts,n21);
	n31.Id = 3; n31.Timestep = 3; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n31.features["com"] = com; n31.features["divProb"] = divProb; n31.features["count"] = count;
	add(ts,n31);

	n12.Id = 4; n12.Timestep = 1; com[0] = 5; com[1] = 5; com[2] = 5; divProb[0] = 0.1; count[0] = 1;
	n12.features["com"] = com; n12.features["divProb"] = divProb; n12.features["count"] = count;
	add(ts,n12);
	n22.Id = 5; n22.Timestep = 2; com[0] = 5; com[1] = 5; com[2] = 5; divProb[0] = 0.1; count[0] = 1;
	n22.features["com"] = com; n22.features["divProb"] = divProb; n22.features["count"] = count;
	add(ts,n22);
	n32.Id = 6; n32.Timestep = 3; com[0] = 5; com[1] = 5; com[2] = 5; divProb[0] = 0.1; count[0] = 1;
	n32.features["com"] = com; n32.features["divProb"] = divProb; n32.features["count"] = count;
	add(ts,n32);

	n13.Id = 7; n13.Timestep = 1; com[0] = 10; com[1] = 10; com[2] = 10; divProb[0] = 0.1; count[0] = 10;
	n13.features["com"] = com; n13.features["divProb"] = divProb; n13.features["count"] = count;
	add(ts,n13);
	n23.Id = 8; n23.Timestep = 2; com[0] = 10; com[1] = 10; com[2] = 10; divProb[0] = 0.1; count[0] = 10;
	n23.features["com"] = com; n23.features["divProb"] = divProb; n23.features["count"] = count;
	add(ts,n23);
	n33.Id = 9; n33.Timestep = 3; com[0] = 10; com[1] = 10; com[2] = 10; divProb[0] = 0.1; count[0] = 10;
	n33.features["com"] = com; n33.features["divProb"] = divProb; n33.features["count"] = count;
	add(ts,n33);

	n24.Id = 10; n24.Timestep = 2; com[0] = 15; com[1] = 15; com[2] = 15; divProb[0] = 0.1; count[0] = 1;
	n24.features["com"] = com; n24.features["divProb"] = divProb; n24.features["count"] = count;
	add(ts,n24);
	n34.Id = 11; n34.Timestep = 3; com[0] = 15; com[1] = 15; com[2] = 15; divProb[0] = 0.1; count[0] = 1;
	n34.features["com"] = com; n34.features["divProb"] = divProb; n34.features["count"] = count;
	add(ts,n34);

	n35.Id = 12; n35.Timestep = 3; com[0] = 20; com[1] = 20; com[2] = 20; divProb[0] = 0.1; count[0] = 1;
	n35.features["com"] = com; n35.features["divProb"] = divProb; n35.features["count"] = count;
	add(ts,n35);


	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

	FieldOfView fov(0, 0, 0, 0, 4, 5, 5, 5); // tlow, xlow, ylow, zlow, tup, xup, yup, zup
	ConsTracking tracking = ConsTracking(
					     2, // max_number_objects
					     true, // detection_by_volume
					     double(1), // avg_obj_size
					     2, // max_neighbor_distance
					     true, //with_divisions
					     0.3, // division_threshold
					     "none", // random_forest_filename
					     fov
					     );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.0, // ep_gap
							    false, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    10., // disappearance_cost,
							    10., // appearance_cost
							    false, //with_merger_resolution
							    3, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );

    size_t t = 0;
	size_t num_events = 0;
	size_t moves = 0;
	size_t apps = 0;
	size_t mergers = 0;
    for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin(); it_t != events.end(); ++it_t) {
		for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {
			Event e = *it;
            cout << e << endl;
			++num_events;
			BOOST_CHECK_NE(e.type, Event::Disappearance);
			BOOST_CHECK_NE(e.type, Event::Division);
			if (e.type == Event::Move) {
                ++moves;
			} else if (e.type == Event::Appearance) {
				++apps;
				BOOST_CHECK(e.traxel_ids[0] == 10 || e.traxel_ids[0] == 12);
			} else if (e.type == Event::Merger) {
				++mergers;
                if (t == 0) {
                    BOOST_CHECK(e.traxel_ids[0] == 7);
                } else if (t==1){
                    BOOST_CHECK(e.traxel_ids[0] == 8);
                } else if (t==2) {
                    BOOST_CHECK(e.traxel_ids[0] == 9);
                }
			}
		}
        ++t;
	}
    BOOST_CHECK_EQUAL(num_events, 12);
	BOOST_CHECK_EQUAL(moves,7);
	BOOST_CHECK_EQUAL(apps,2);
    BOOST_CHECK_EQUAL(mergers,3);
}

BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_AppearanceSimple ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2        3
	//
	//  1 ------ 1 ------ 1
	//  1 ------ 1-\----- 1
	TraxelStore ts;
	Traxel n11, n21, n31;
	Traxel n12, n22, n32;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	feature_array count(feature_array::difference_type(1));
	n11.Id = 1; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n11.features["com"] = com; n11.features["divProb"] = divProb; n11.features["count"] = count;
	add(ts,n11);
	n21.Id = 2; n21.Timestep = 2; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n21.features["com"] = com; n21.features["divProb"] = divProb; n21.features["count"] = count;
	add(ts,n21);
	n31.Id = 3; n31.Timestep = 3; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n31.features["com"] = com; n31.features["divProb"] = divProb; n31.features["count"] = count;
	add(ts,n31);

	n12.Id = 4; n12.Timestep = 1; com[0] = 5; com[1] = 5; com[2] = 5; divProb[0] = 0.1; count[0] = 1;
	n12.features["com"] = com; n12.features["divProb"] = divProb; n12.features["count"] = count;
	add(ts,n12);
	n22.Id = 5; n22.Timestep = 2; com[0] = 5; com[1] = 5; com[2] = 5; divProb[0] = 0.1; count[0] = 1;
	n22.features["com"] = com; n22.features["divProb"] = divProb; n22.features["count"] = count;
	add(ts,n22);
	n32.Id = 6; n32.Timestep = 3; com[0] = 2; com[1] = 1; com[2] = 1; divProb[0] = 0.1; count[0] = 1;
	n32.features["com"] = com; n32.features["divProb"] = divProb; n32.features["count"] = count;
	add(ts,n32);

	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

	FieldOfView fov(0, 0, 0, 0, 4, 5, 5, 5); // tlow, xlow, ylow, zlow, tup, xup, yup, zup
	ConsTracking tracking = ConsTracking(
					     1, // max_number_objects
					     true, // detection_by_volume
					     double(1), // avg_obj_size
					     100, // max_neighbor_distance
					     true, //with_divisions
					     0.3, // division_threshold
					     "none", // random_forest_filename
					     fov
					     );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.0, // ep_gap
							    false, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    1500., // disappearance_cost,
							    1500., // appearance_cost
							    false, //with_merger_resolution
							    3, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );

	size_t t = 1;
	size_t num_events = 0;
    for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin()+1; it_t != events.end(); ++it_t) {
		for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {
			Event e = *it;
			++num_events;
			BOOST_CHECK_NE(e.type, Event::Disappearance);
			BOOST_CHECK_NE(e.type, Event::Division);
			BOOST_CHECK_NE(e.type, Event::Appearance);
			BOOST_CHECK_NE(e.type, Event::Merger);
			if (e.type == Event::Move) {
				if (t == 1) {
					BOOST_CHECK((e.traxel_ids[0] == 1 && e.traxel_ids[1] == 2) ||
							(e.traxel_ids[0] == 4 && e.traxel_ids[1] == 5));
				} else if (t== 2) {
					BOOST_CHECK((e.traxel_ids[0] == 2 && e.traxel_ids[1] == 3) ||
							(e.traxel_ids[0] == 5 && e.traxel_ids[1] == 6));
				}
			}
		}
		++t;
	}
	BOOST_CHECK_EQUAL(num_events, 4);
}



BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_Tracklets ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2       3
	//  1
	//    |
	//    |
	//      ---- 1 ----  1
	//    |
	//  2
	//  1 ----   2 ----- 1
	//  0 ------ 1 ----- 0
	TraxelStore ts;
	Traxel n11, n12, n21, n31, n13, n22, n32;
	Traxel n14, n24, n34;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	feature_array count(feature_array::difference_type(1));
	n11.Id = 1; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n11.features["com"] = com; n11.features["divProb"] = divProb; n11.features["count"] = count;
	add(ts,n11);
	n12.Id = 3; n12.Timestep = 1; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1; count[0] = 3;
	n12.features["com"] = com; n12.features["divProb"] = divProb; n12.features["count"] = count;
	add(ts,n12);
	n21.Id = 10; n21.Timestep = 2; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1; count[0] = 1;
	n21.features["com"] = com; n21.features["divProb"] = divProb; n21.features["count"] = count;
	add(ts,n21);
	n31.Id = 11; n31.Timestep = 3; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1; count[0] = 1;
	n31.features["com"] = com; n31.features["divProb"] = divProb; n31.features["count"] = count;
	add(ts,n31);

	n13.Id = 12; n13.Timestep = 1; com[0] = 100; com[1] = 100; com[2] = 100; divProb[0] = 0.1; count[0] = 1;
	n13.features["com"] = com; n13.features["divProb"] = divProb; n13.features["count"] = count;
	add(ts,n13);
	n22.Id = 13; n22.Timestep = 2; com[0] = 100; com[1] = 100; com[2] = 100; divProb[0] = 0.1; count[0] = 3;
	n22.features["com"] = com; n22.features["divProb"] = divProb; n22.features["count"] = count;
	add(ts,n22);
	n32.Id = 14; n32.Timestep = 3; com[0] = 100; com[1] = 100; com[2] = 100; divProb[0] = 0.1; count[0] = 1;
	n32.features["com"] = com; n32.features["divProb"] = divProb; n32.features["count"] = count;
	add(ts,n32);

	n14.Id = 15; n14.Timestep = 1; com[0] = 200; com[1] = 100; com[2] = 100; divProb[0] = 0.1; count[0] = 0.1;
	n14.features["com"] = com; n14.features["divProb"] = divProb; n14.features["count"] = count;
	add(ts,n14);
	n24.Id = 16; n24.Timestep = 2; com[0] = 200; com[1] = 100; com[2] = 100; divProb[0] = 0.1; count[0] = 1;
	n24.features["com"] = com; n24.features["divProb"] = divProb; n24.features["count"] = count;
	add(ts,n24);
	n34.Id = 17; n34.Timestep = 3; com[0] = 200; com[1] = 100; com[2] = 100; divProb[0] = 0.1; count[0] = 0.1;
	n34.features["com"] = com; n34.features["divProb"] = divProb; n34.features["count"] = count;
	add(ts,n34);

	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

	FieldOfView fov(0, 0, 0, 0, 4, 5, 5, 5); // tlow, xlow, ylow, zlow, tup, xup, yup, zup
	ConsTracking tracking = ConsTracking(
					     2, // max_number_objects
					     true, // detection_by_volume
					     double(1.1), // avg_obj_size
					     20, // max_neighbor_distance
					     true, //with_divisions
					     0.3, // division_threshold
					     "none", // random_forest_filename
					     fov
					     );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.0, // ep_gap
							    true, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    1500., // disappearance_cost,
							    1500., // appearance_cost
							    false, //with_merger_resolution
							    3, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );

	size_t moves = 0;
	size_t t = 1;
    for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin()+1; it_t != events.end(); ++it_t) {
		for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {
			Event e = *it;
			if (e.type == Event::Move) {
				++moves;
				BOOST_CHECK( (e.traxel_ids[0] == 1 && e.traxel_ids[1] == 10) ||
						 (e.traxel_ids[0] == 10 && e.traxel_ids[1] == 11) ||
						 (e.traxel_ids[0] == 12 && e.traxel_ids[1] == 13) ||
						 (e.traxel_ids[0] == 13 && e.traxel_ids[1] == 14) ||
						 (e.traxel_ids[0] == 3 && e.traxel_ids[1] == 10) ||
						 (e.traxel_ids[0] == 15 && e.traxel_ids[1] == 16) ||
						 (e.traxel_ids[0] == 16 && e.traxel_ids[1] == 17));
			}
			if (e.type == Event::Disappearance) {
				BOOST_CHECK(e.traxel_ids[0] == 3 ||
						e.traxel_ids[0] == 1);
			}
//			if (e.type == Event::Merger) {
//				BOOST_CHECK(e.traxel_ids[0] == 3 ||
//							e.traxel_ids[0] == 12 ||
//							e.traxel_ids[0] == 14);
//			}
		}
		++t;
	}
	BOOST_CHECK_EQUAL(moves, 6);
}


BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_Merger3 ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2       3
	//  1				1
	//    |		       |
	//      ---- 2 ----  		1
	//    |			   |	  |
	//  1				1 ----
	//  					  |
	//							0
	TraxelStore ts;
	Traxel n11, n12, n21, n31, n32, n41, n42, n43;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	feature_array detProb(feature_array::difference_type(3));
	feature_array count(feature_array::difference_type(1));
	feature_array coordinates(feature_array::difference_type(2*3));

	n11.Id = 11; n11.Timestep = 1; com[0] = 854; com[1] = 718; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.01; detProb[1] = 0.98; detProb[2] = 0.01;
	n11.features["com"] = com; n11.features["divProb"] = divProb; n11.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n11.features["coordinates"] = coordinates;
	add(ts,n11);

	n12.Id = 12; n12.Timestep = 1; com[0] = 846; com[1] = 751; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.01; detProb[1] = 0.98; detProb[2] = 0.01;
	n12.features["com"] = com; n12.features["divProb"] = divProb; n12.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n12.features["coordinates"] = coordinates;
	add(ts,n12);

	n21.Id = 21; n21.Timestep = 2; com[0] = 845; com[1] = 732; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.01; detProb[1] = 0.01; detProb[2] = 0.98;
	n21.features["com"] = com; n21.features["divProb"] = divProb; n21.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n21.features["coordinates"] = coordinates;
	add(ts,n21);

	n31.Id = 31; n31.Timestep = 3; com[0] = 853; com[1] = 724; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.01; detProb[1] = 0.98; detProb[2] = 0.01;
	n31.features["com"] = com; n31.features["divProb"] = divProb; n31.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n31.features["coordinates"] = coordinates;
	add(ts,n31);

	n32.Id = 32; n32.Timestep = 3; com[0] = 835; com[1] = 747; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.01; detProb[1] = 0.98; detProb[2] = 0.01;
	n32.features["com"] = com; n32.features["divProb"] = divProb; n32.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n32.features["coordinates"] = coordinates;
	add(ts,n32);

	n41.Id = 41; n41.Timestep = 4; com[0] = 825; com[1] = 753; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.01; detProb[1] = 0.98; detProb[2] = 0.01;
	n41.features["com"] = com; n41.features["divProb"] = divProb; n41.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n41.features["coordinates"] = coordinates;
	add(ts,n41);

	n42.Id = 42; n42.Timestep = 4; com[0] = 842; com[1] = 761; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.98; detProb[1] = 0.01; detProb[2] = 0.01;
	n42.features["com"] = com; n42.features["divProb"] = divProb; n42.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n42.features["coordinates"] = coordinates;
	add(ts,n42);

	n43.Id = 43; n43.Timestep = 4; com[0] = 853; com[1] = 724; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.01; detProb[1] = 0.98; detProb[2] = 0.01;
	n43.features["com"] = com; n43.features["divProb"] = divProb; n43.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n43.features["coordinates"] = coordinates;
	add(ts,n43);


	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

    FieldOfView fov(0, 0, 0, 0, 4, 1000, 1000, 1); // tlow, xlow, ylow, zlow, tup, xup, yup, zup

    ConsTracking tracking = ConsTracking(
					 2, // max_number_objects
					 false, // detection_by_volume
					 double(1.1), // avg_obj_size
					 99999, // max_neighbor_distance
					 true, //with_divisions
					 0.1, // division_threshold
					 "none", // random_forest_filename
					 fov
					 );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.05, // ep_gap
							    true, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    500., // disappearance_cost,
							    500., // appearance_cost
							    true, //with_merger_resolution
							    2, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );

	size_t t = 1;
	size_t apps = 0;
	size_t disapps = 0;
	size_t num_events = 0;
    for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin()+1; it_t != events.end(); ++it_t) {
		for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {
			Event e = *it;
			if (e.type == Event::Appearance) {
				++apps;
			}
			if (e.type == Event::Disappearance) {
				++disapps;
			}
			cout << e << endl;
			++num_events;
		}
		++t;
	}

	BOOST_CHECK_EQUAL(apps, 0);
	BOOST_CHECK_EQUAL(disapps, 0);
	BOOST_CHECK_GE(num_events, 3);
}

namespace{
void constructTraxel(Traxel& n, size_t id, size_t timestep, double pDiv, std::vector<int> center,
		std::vector<int> shift, std::vector<double> pDet) {
	feature_array com(feature_array::difference_type(3));
	feature_array com_corrected(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	feature_array detProb(feature_array::difference_type(pDet.size()));
	feature_array count(feature_array::difference_type(1));
	feature_array coordinates(feature_array::difference_type(2*3));
	n.Id = id; n.Timestep = timestep;
	divProb[0] = pDiv;
	for(size_t i = 0; i<com.size(); ++i) {
		com[i] = center[i];
		com_corrected[i] = center[i] + shift[i];
	}
	for(size_t i = 0; i<pDet.size(); ++i) {
		detProb[i] = pDet[i];
	}
	n.features["com"] = com;
	n.features["divProb"] = divProb;
	n.features["detProb"] = detProb;
	n.features["com_corrected"] = com_corrected;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n.features["coordinates"] = coordinates;
}

template <class T>
void pushBackMultiple(std::vector<T>& vec, const T v1, const T v2) {
	vec.push_back(v1);
	vec.push_back(v2);
}

template <class T>
void pushBackMultiple(std::vector<T>& vec, const T v1, const T v2, const T v3) {
	vec.push_back(v1);
	vec.push_back(v2);
	vec.push_back(v3);
}

}

BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_TranslationVector_Traxels ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2       3
	//  11 ------21 -----31
	//    |
	//     ----- 22 -----32
	//
	//  12 ------23 -----33
	TraxelStore ts;
	Traxel n11, n12, n21, n22, n23, n31, n32, n33;
	std::vector<int> center;
	std::vector<int> shift;
	std::vector<double> pDet;

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 100, 100, 0);
	pushBackMultiple(shift, 100,100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n11,11,1,0.00,center,shift,pDet);
	add(ts,n11);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 200, 200, 0);
	pushBackMultiple(shift, 100,100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n12,12,1,0.00,center,shift,pDet);
	add(ts,n12);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 100, 100, 0);
	pushBackMultiple(shift, -100,-100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n21,21,2,0.00,center,shift,pDet);
	add(ts,n21);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 200, 200, 0);
	pushBackMultiple(shift, -100,-100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n22,22,2,0.00,center,shift,pDet);
	add(ts,n22);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 300, 300, 0);
	pushBackMultiple(shift, -100,-100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n23,23,2,0.00,center,shift,pDet);
	add(ts,n23);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 0, 0, 0);
	pushBackMultiple(shift, 100,100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n31,31,3,0.00,center,shift,pDet);
	add(ts,n31);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 100, 100, 0);
	pushBackMultiple(shift, -100,-100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n32,32,3,0.00,center,shift,pDet);
	add(ts,n32);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 200, 200, 0);
	pushBackMultiple(shift, 100,100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n33,33,3,0.00,center,shift,pDet);
	add(ts,n33);


	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

    FieldOfView fov(0, 0, 0, 0, 3, 1000, 1000, 1); // tlow, xlow, ylow, zlow, tup, xup, yup, zup

    ConsTracking tracking = ConsTracking(
					 1, // max_number_objects
					 false, // detection_by_volume
					 double(1.1), // avg_obj_size
					 301, // max_neighbor_distance
					 false, //with_divisions
					 0.5, // division_threshold
					 "none", // random_forest_filename
					 fov
					 );

    std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.05, // ep_gap
							    false, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    10., // disappearance_cost,
							    10., // appearance_cost
							    false, //with_merger_resolution
							    2, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );


    size_t t = 1;
	BOOST_CHECK_EQUAL(events[t].size(), 3);
	for (std::vector<Event>::const_iterator it = events[t].begin(); it!=events[t].end(); ++it) {
			Event e = *it;
			if (e.type == Event::Move && e.traxel_ids[0] == 11) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 22);
			} else if (e.type == Event::Move && e.traxel_ids[0] == 12) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 23);
			} else if (e.type == Event::Appearance) {
				BOOST_CHECK_EQUAL(e.traxel_ids[0], 21);
			} else {
				cout << "unexpected event: " << e;
				BOOST_CHECK(false);
			}
	}
    t=2;
	BOOST_CHECK_EQUAL(events[t].size(), 3);
	for (std::vector<Event>::const_iterator it = events[t].begin(); it!=events[t].end(); ++it) {
			Event e = *it;
			if (e.type == Event::Move && e.traxel_ids[0] == 22) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 32);
			} else if (e.type == Event::Move && e.traxel_ids[0] == 23) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 33);
			} else if (e.type == Event::Move && e.traxel_ids[0] == 21) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 31);
			} else {
				cout << "unexpected event: " << e;
				BOOST_CHECK(false);
			}
	}

}


BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_TranslationVector_Tracklets ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2       3
	//  11 ------21 -----31
	//    |
	//     ----- 22 -----32
	//
	//  12 ------23 -----33
	TraxelStore ts;
	Traxel n11, n12, n21, n22, n23, n31, n32, n33;
	std::vector<int> center;
	std::vector<int> shift;
	std::vector<double> pDet;

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 100, 100, 0);
	pushBackMultiple(shift, 100,100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n11,11,1,0.00,center,shift,pDet);
	add(ts,n11);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 200, 200, 0);
	pushBackMultiple(shift, 100,100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n12,12,1,0.00,center,shift,pDet);
	add(ts,n12);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 100, 100, 0);
	pushBackMultiple(shift, -100,-100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n21,21,2,0.00,center,shift,pDet);
	add(ts,n21);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 200, 200, 0);
	pushBackMultiple(shift, -100,-100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n22,22,2,0.00,center,shift,pDet);
	add(ts,n22);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 300, 300, 0);
	pushBackMultiple(shift, -100,-100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n23,23,2,0.00,center,shift,pDet);
	add(ts,n23);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 0, 0, 0);
	pushBackMultiple(shift, 200,200,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n31,31,3,0.00,center,shift,pDet);
	add(ts,n31);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 100, 100, 0);
	pushBackMultiple(shift, -100,-100,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n32,32,3,0.00,center,shift,pDet);
	add(ts,n32);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 300, 300, 0);
	pushBackMultiple(shift, -200,-200,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n33,33,3,0.00,center,shift,pDet);
	add(ts,n33);


	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

    FieldOfView fov(0, 0, 0, 0, 3, 1000, 1000, 1); // tlow, xlow, ylow, zlow, tup, xup, yup, zup

    ConsTracking tracking = ConsTracking(
					 1, // max_number_objects
					 false, // detection_by_volume
					 double(1.1), // avg_obj_size
					 301, // max_neighbor_distance
					 false, //with_divisions
					 0.5, // division_threshold
					 "none", // random_forest_filename
					 fov
					 );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.05, // ep_gap
							    true, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    10., // disappearance_cost,
							    10., // appearance_cost
							    false, //with_merger_resolution
							    2, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );


    size_t t = 1;
	BOOST_CHECK_EQUAL(events[t].size(), 3);
	for (std::vector<Event>::const_iterator it = events[t].begin(); it!=events[t].end(); ++it) {
			Event e = *it;
			if (e.type == Event::Move && e.traxel_ids[0] == 11) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 22);
			} else if (e.type == Event::Move && e.traxel_ids[0] == 12) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 23);
			} else if (e.type == Event::Appearance) {
				BOOST_CHECK_EQUAL(e.traxel_ids[0], 21);
			} else {
				cout << "unexpected event: " << e;
				BOOST_CHECK(false);
			}
	}
    t=2;
	BOOST_CHECK_EQUAL(events[t].size(), 3);
	for (std::vector<Event>::const_iterator it = events[t].begin(); it!=events[t].end(); ++it) {
			Event e = *it;
			if (e.type == Event::Move && e.traxel_ids[0] == 22) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 32);
			} else if (e.type == Event::Move && e.traxel_ids[0] == 23) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 33);
			} else if (e.type == Event::Move && e.traxel_ids[0] == 21) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 31);
			} else {
				cout << "unexpected event: " << e;
				BOOST_CHECK(false);
			}
	}

}


BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_Merger4 ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2       3
	//  1
	//    |
	//      ---- 2 ---- 2
	//    |
	//  1
	TraxelStore ts;
	Traxel n11, n12, n21, n31;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	feature_array detProb(feature_array::difference_type(3));
	feature_array count(feature_array::difference_type(1));
	feature_array coordinates(feature_array::difference_type(2*3));

	n11.Id = 11; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.01; detProb[1] = 0.98; detProb[2] = 0.01;
	n11.features["com"] = com; n11.features["divProb"] = divProb; n11.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n11.features["coordinates"] = coordinates;
	add(ts,n11);

	n12.Id = 12; n12.Timestep = 1; com[0] = 2; com[1] = 2; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.01; detProb[1] = 0.98; detProb[2] = 0.01;
	n12.features["com"] = com; n12.features["divProb"] = divProb; n12.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n12.features["coordinates"] = coordinates;
	add(ts,n12);

	n21.Id = 21; n21.Timestep = 2; com[0] = 1; com[1] = 1; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.01; detProb[1] = 0.01; detProb[2] = 0.98;
	n21.features["com"] = com; n21.features["divProb"] = divProb; n21.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n21.features["coordinates"] = coordinates;
	add(ts,n21);

	n31.Id = 31; n31.Timestep = 3; com[0] = 1; com[1] = 1; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.01; detProb[1] = 0.01; detProb[2] = 0.98;
	n31.features["com"] = com; n31.features["divProb"] = divProb; n31.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n31.features["coordinates"] = coordinates;
	add(ts,n31);


	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

    FieldOfView fov(0, 0, 0, 0, 4, 1000, 1000, 1); // tlow, xlow, ylow, zlow, tup, xup, yup, zup

    ConsTracking tracking = ConsTracking(
					 2, // max_number_objects
					 false, // detection_by_volume
					 double(1.1), // avg_obj_size
					 99999, // max_neighbor_distance
					 true, //with_divisions
					 0.1, // division_threshold
					 "none", // random_forest_filename
					 fov
					 );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.0, // ep_gap
							    true, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    1500., // disappearance_cost,
							    1500., // appearance_cost
							    true, //with_merger_resolution
							    2, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );

	size_t t = 1;
	size_t apps = 0;
	size_t disapps = 0;
	size_t moves = 0;
	size_t mergers = 0;
	size_t multis = 0;
	size_t num_events = 0;
    for (std::vector< std::vector<Event> >::const_iterator it_t = events.begin()+1; it_t != events.end(); ++it_t) {
		for (std::vector<Event>::const_iterator it = (*it_t).begin(); it!=(*it_t).end(); ++it) {
			Event e = *it;
			if (e.type == Event::Appearance) {
				++apps;
			}
			if (e.type == Event::Disappearance) {
				++disapps;
			}
			if (e.type == Event::MultiFrameMove) {
				++multis;
			}
			if (e.type == Event::Merger) {
				++mergers;
			}
			if (e.type == Event::Move) {
				++moves;
			}
			cout << e << endl;
			++num_events;
		}
		++t;
	}

	BOOST_CHECK_EQUAL(apps, 0);
	BOOST_CHECK_EQUAL(disapps, 0);
	BOOST_CHECK_EQUAL(moves, 3);
    BOOST_CHECK_EQUAL(mergers, 2);
	BOOST_CHECK_EQUAL(multis, 0);
}




BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_MergerResolvingDivision ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	//  t=1      2
	//  1 ------ 2
	//      |
	//  1D<
	//      |
	//  1 ------ 2
	TraxelStore ts;
	Traxel n11, n12, n13, n21, n22;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	feature_array detProb(feature_array::difference_type(3));
	feature_array count(feature_array::difference_type(1));
	feature_array coordinates(feature_array::difference_type(2*3));

	n11.Id = 11; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.0; detProb[1] = 1.; detProb[2] = 0.0;
	n11.features["com"] = com; n11.features["divProb"] = divProb; n11.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n11.features["coordinates"] = coordinates;
	add(ts,n11);

	n12.Id = 12; n12.Timestep = 1; com[0] = 2; com[1] = 2; com[2] = 0; divProb[0] = 1.;
	detProb[0] = 0.0; detProb[1] = 1.; detProb[2] = 0.0;
	n12.features["com"] = com; n12.features["divProb"] = divProb; n12.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n12.features["coordinates"] = coordinates;
	add(ts,n12);

	n13.Id = 13; n13.Timestep = 1; com[0] = 4; com[1] = 4; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.00; detProb[1] = 1.; detProb[2] = 0.0;
	n13.features["com"] = com; n13.features["divProb"] = divProb; n13.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 1;
	n13.features["coordinates"] = coordinates;
	add(ts,n13);

	n21.Id = 21; n21.Timestep = 2; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.0; detProb[1] = 0.0; detProb[2] = 1.;
	n21.features["com"] = com; n21.features["divProb"] = divProb; n21.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 0;
	n21.features["coordinates"] = coordinates;
	add(ts,n21);

	n22.Id = 22; n22.Timestep = 2; com[0] = 4; com[1] = 4; com[2] = 0; divProb[0] = 0.01;
	detProb[0] = 0.0; detProb[1] = 0.0; detProb[2] = 1.;
	n22.features["com"] = com; n22.features["divProb"] = divProb; n22.features["detProb"] = detProb;
	coordinates[0] = com[0]; coordinates[1] = com[1]; coordinates[2] = com[2];
	coordinates[3] = com[0] + 1; coordinates[4] = com[1] + 1; coordinates[5] = com[2] + 0;
	n22.features["coordinates"] = coordinates;
	add(ts,n22);

	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

    FieldOfView fov(0, 0, 0, 0, 4, 1000, 1000, 1); // tlow, xlow, ylow, zlow, tup, xup, yup, zup

    ConsTracking tracking = ConsTracking(
					 2, // max_number_objects
					 false, // detection_by_volume
					 double(1.1), // avg_obj_size
					 99999, // max_neighbor_distance
					 true, //with_divisions
					 0.1, // division_threshold
					 "none", // random_forest_filename
					 fov,
					 "none"
					 );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.0, // ep_gap
							    true, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    1500., // disappearance_cost,
							    1500., // appearance_cost
							    true, //with_merger_resolution
							    2, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );

    size_t t = 1;
    BOOST_CHECK_EQUAL(events[t].size(), 5);
    size_t num_mergers = 0;
	for (std::vector<Event>::const_iterator it = events[t].begin(); it!=events[t].end(); ++it) {
        Event e = *it;
        if (e.type == Event::Move && e.traxel_ids[0] == 11) {
            BOOST_CHECK_EQUAL(e.traxel_ids[1], 21);
        } else if (e.type == Event::Move && e.traxel_ids[0] == 13) {
            BOOST_CHECK_EQUAL(e.traxel_ids[1], 22);
        } else if (e.type == Event::Division && e.traxel_ids[0] == 12) {
            set<unsigned> division_set(e.traxel_ids.begin()+1, e.traxel_ids.end());
            set<unsigned> comparison_set;
            comparison_set.insert(21);
            comparison_set.insert(22);
            BOOST_CHECK_EQUAL_COLLECTIONS(division_set.begin(),
                                         division_set.end(),
                                         comparison_set.begin(),
                                         comparison_set.end());
        } else if (e.type == Event::Merger) {
            ++num_mergers;
            BOOST_CHECK(e.traxel_ids[0] == 21 || e.traxel_ids[0] == 22);
        } else  {
            cout << "unexpected event: " << e;
            BOOST_CHECK(false);
        }
    }
    BOOST_CHECK_EQUAL(num_mergers, 2);

}




BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_TranslationVector2 ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	TraxelStore ts;
	Traxel n11, n12;
	Traxel n21, n22;
	Traxel n31, n32;
	Traxel n41, n42;
	std::vector<int> center;
	std::vector<int> shift;
	std::vector<double> pDet;

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 90, 267, 0);
	pushBackMultiple(shift, -7,3,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n11,11,1,0.00,center,shift,pDet);
	add(ts,n11);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 113, 312, 0);
	pushBackMultiple(shift, -7,3,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n12,12,1,0.00,center,shift,pDet);
	add(ts,n12);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 83, 270, 0);
	pushBackMultiple(shift, 7,-3,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n21,21,2,0.00,center,shift,pDet);
	add(ts,n21);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 106, 315, 0);
	pushBackMultiple(shift, 7,-3,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n22,22,2,0.00,center,shift,pDet);
	add(ts,n22);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 90, 267, 0);
	pushBackMultiple(shift, 27,46,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n31,31,3,0.00,center,shift,pDet);
	add(ts,n31);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 113, 312, 0);
	pushBackMultiple(shift, 27,46,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n32,32,3,0.00,center,shift,pDet);
	add(ts,n32);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 117, 313, 0);
	pushBackMultiple(shift, 0,0,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n41,41,4,0.00,center,shift,pDet);
	add(ts,n41);

	center.clear(); shift.clear(); pDet.clear();
	pushBackMultiple(center, 140, 358, 0);
	pushBackMultiple(shift, 0,0,0);
	pushBackMultiple(pDet, 0.,1.);
	constructTraxel(n42,42,4,0.00,center,shift,pDet);
	add(ts,n42);

	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

    FieldOfView fov(0, 0, 0, 0, 3, 1000, 1000, 1); // tlow, xlow, ylow, zlow, tup, xup, yup, zup

    ConsTracking tracking = ConsTracking(
					 1, // max_number_objects
					 false, // detection_by_volume
					 double(1.1), // avg_obj_size
					 301, // max_neighbor_distance
					 false, //with_divisions
					 0.5, // division_threshold
					 "none", // random_forest_filename
					 fov
					 );

	std::cout << "Run Conservation tracking" << std::endl;
	std::cout << std::endl;
	std::vector< std::vector<Event> > events = tracking(ts,
							    0, // forbidden_cost
							    0.05, // ep_gap
							    true, // with_tracklets
							    10.0, //division_weight
							    10.0, //transition_weight
							    10., // disappearance_cost,
							    10., // appearance_cost
							    false, //with_merger_resolution
							    2, //n_dim
							    5, //transition_parameter
							    0 //border_width for app/disapp costs
							    );


    size_t t = 1;
	BOOST_CHECK_EQUAL(events[t].size(), 2);
	for (std::vector<Event>::const_iterator it = events[t].begin(); it!=events[t].end(); ++it) {
			Event e = *it;
			if (e.type == Event::Move && e.traxel_ids[0] == 11) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 21);
			} else if (e.type == Event::Move && e.traxel_ids[0] == 12) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 22);
			} else {
				cout << "unexpected event: " << e;
				BOOST_CHECK(false);
			}
	}
    t=2;
	BOOST_CHECK_EQUAL(events[t].size(), 2);
	for (std::vector<Event>::const_iterator it = events[t].begin(); it!=events[t].end(); ++it) {
			Event e = *it;
			if (e.type == Event::Move && e.traxel_ids[0] == 21) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 31);
			} else if (e.type == Event::Move && e.traxel_ids[0] == 22) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 32);
			} else {
				cout << "unexpected event: " << e;
				BOOST_CHECK(false);
			}
	}
    t=3;
	BOOST_CHECK_EQUAL(events[t].size(), 2);
	for (std::vector<Event>::const_iterator it = events[t].begin(); it!=events[t].end(); ++it) {
			Event e = *it;
			if (e.type == Event::Move && e.traxel_ids[0] == 31) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 41);
			} else if (e.type == Event::Move && e.traxel_ids[0] == 32) {
				BOOST_CHECK_EQUAL(e.traxel_ids[1], 42);
			} else {
				cout << "unexpected event: " << e;
				BOOST_CHECK(false);
			}
	}

}
