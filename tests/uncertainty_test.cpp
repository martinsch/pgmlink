#define BOOST_TEST_MODULE hypotheses_test

#include <vector>
#include <string>
#include <iostream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <lemon/core.h>
#include <lemon/concepts/digraph.h>
#include <lemon/list_graph.h>
#include <lemon/maps.h>

#include "pgmlink/hypotheses.h"
#include "pgmlink/traxels.h"
#include "pgmlink/tracking.h"
#include "pgmlink/reasoner_constracking.h"

using namespace pgmlink;
using namespace std;
using namespace boost;

BOOST_AUTO_TEST_CASE( uncertainty ) {
  
	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	typedef HypothesesGraph::ArcIt ArcIt2;
	typedef HypothesesGraph::Arc Arc;
	typedef HypothesesGraph::NodeIt NodeIt;
	typedef HypothesesGraph::Node Node;
	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;
	
	//  t=1       2       3
	//   o                 o
	//    |               |
	//     ------ o ------
	//    |               |
	//   o                 o
	TraxelStore ts;
	Traxel n11, n12, n21, n31, n32;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	feature_array detProb(feature_array::difference_type(2));
	//detProb[2]=0;
	n11.Id = 11; n11.Timestep = 1; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0; detProb[0] = 0.4;detProb[1]=0.6;
	n11.features["com"] = com; n11.features["divProb"] = divProb; n11.features["detProb"] = detProb;
	add(ts,n11);
	n12.Id = 12; n12.Timestep = 1; com[0] = 3; com[1] = 2; com[2] = 3; divProb[0] = 0; detProb[0] = 0.6;detProb[1]=0.4;
	n12.features["com"] = com; n12.features["divProb"] = divProb; n12.features["detProb"] = detProb;
	add(ts,n12);
	
	n21.Id = 21; n21.Timestep = 2; com[0] = 2; com[1] = 2; com[2] = 3; divProb[0] = 0.39; detProb[0] = 0.1;detProb[1]=0.9;
	n21.features["com"] = com; n21.features["divProb"] = divProb; n21.features["detProb"] = detProb;
	add(ts,n21);
	
	n31.Id = 31; n31.Timestep = 3; com[0] = 2; com[1] = 1; com[2] = 1; divProb[0] = 0; detProb[0] = 0.6;detProb[1]=0.4;
	n31.features["com"] = com; n31.features["divProb"] = divProb; n31.features["detProb"] = detProb;
	add(ts,n31);
	n32.Id = 32; n32.Timestep = 3; com[0] = 3; com[1] = 1; com[2] = 1; divProb[0] = 0; detProb[0] = 0.2;detProb[1]=0.8;
	n32.features["com"] = com; n32.features["divProb"] = divProb; n32.features["detProb"] = detProb;
	add(ts,n32);

	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

	FieldOfView fov(0, 0, 0, 0, 3, 5, 5, 5); // tlow, xlow, ylow, zlow, tup, xup, yup, zup
	ConsTracking tracking = ConsTracking(
				  1, // max_number_objects
	              20, // max_neighbor_distance
				  0.3, // division_threshold
				  "none", // random_forest_filename
	              false, // detection_by_volume
	              0, // forbidden_cost
	              0.0, // ep_gap
	              double(1.1), // avg_obj_size
				  false, // with_tracklets
				  10.0, //division_weight
				  10.0, //transition_weight
				  true, //with_divisions
				  10., // disappearance_cost,
				  10., // appearance_cost
				  false, //with_merger_resolution
				  3, //n_dim
				  5, //transition_parameter
				  0, //border_width for app/disapp costs
	              fov
		  	      );
	tracking(ts, TimestepIdCoordinateMapPtr(), 100,0,1);
}

// EOF

