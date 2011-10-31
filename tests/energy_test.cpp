#define BOOST_TEST_MODULE energy_test

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "energy.h"

using namespace Tracking;
using namespace std;
using namespace boost;

BOOST_AUTO_TEST_CASE( GeometryDivision_operator )
{
    // perpare mock objects
    Traxel ancestor, descendant1, descendant2;
    feature_array com1(feature_array::difference_type(3));
    feature_array com2(feature_array::difference_type(3));
    feature_array com3(feature_array::difference_type(3));

    com1[0] = 1;
    com1[1] = 2;
    com1[2] = 0;
    ancestor.features["com"] = com1;

    com2[0] = 2;
    com2[1] = 3;
    com2[2] = 7;
    descendant1.features["com"] = com2;

    com3[0] = 1;
    com3[1] = 0.5;
    com3[2] = 1.5;
    descendant2.features["com"] = com3;

    Traxels prev, curr;
    prev[0] = ancestor;
    curr[1] = descendant1;
    curr[2] = descendant2;

    GeometryDivision energy;
    energy(ancestor, descendant1, descendant2, prev, curr);
    //BOOST_CHECK_EQUAL(e, 0);
}

// EOF
