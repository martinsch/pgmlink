#define BOOST_TEST_MODULE field_of_view_test

#include <iostream>
#include <vector>
#include <utility>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "field_of_view.h"

using namespace Tracking;
using namespace std;
using namespace boost;

BOOST_AUTO_TEST_CASE( FieldOfView_test ) {
  FieldOfView fov;

  //
  // bounding box
  //
  BOOST_CHECK(fov.contains( 0,0,0,0 ));
  BOOST_CHECK(!fov.contains( 3,2,2,7 ));
  BOOST_CHECK(!fov.contains( 3,2,39,7));

  vector<double> should_lb;
  should_lb.push_back(0);
  should_lb.push_back(0);
  should_lb.push_back(0);  
  should_lb.push_back(0);

  vector<double> should_ub;
  should_ub.push_back(0);
  should_ub.push_back(0);
  should_ub.push_back(0);  
  should_ub.push_back(0);

  BOOST_CHECK_EQUAL_COLLECTIONS(fov.lower_bound().begin(), fov.lower_bound().end(), should_lb.begin(), should_lb.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(fov.upper_bound().begin(), fov.upper_bound().end(), should_ub.begin(), should_ub.end());

  // field of view: (t 3, x 1, y 1, z 1) -- (t 5, x 10, y 12, z 6)
  fov.set_boundingbox( 3,1,1,1, 5,10,12,6 );

  should_lb.clear();
  should_lb.push_back(3);
  should_lb.push_back(1);
  should_lb.push_back(1);  
  should_lb.push_back(1);

  should_ub.clear();
  should_ub.push_back(5);
  should_ub.push_back(10);
  should_ub.push_back(12);  
  should_ub.push_back(6);

  BOOST_CHECK_EQUAL_COLLECTIONS(fov.lower_bound().begin(), fov.lower_bound().end(), should_lb.begin(), should_lb.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(fov.upper_bound().begin(), fov.upper_bound().end(), should_ub.begin(), should_ub.end());

  BOOST_CHECK(!fov.contains( 0,0,0,0 ));
  BOOST_CHECK(fov.contains( 3,2,2,6 ));
  BOOST_CHECK(!fov.contains( 3,2,39,7));


  //
  // margin
  //
  BOOST_CHECK_EQUAL(fov.spatial_margin( 0,0,0,0 ), 1);
  BOOST_CHECK_EQUAL(fov.spatial_margin( 3,2,1,1 ), 0);
  BOOST_CHECK_EQUAL(fov.spatial_margin( 3,9,11,5 ), 1);
  BOOST_CHECK_EQUAL(fov.temporal_margin( 0,0,0,0 ), 3);
  BOOST_CHECK_EQUAL(fov.temporal_margin( 3,2,2,7 ), 0);
} 
