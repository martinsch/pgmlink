#define BOOST_TEST_MODULE energy_test

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/function.hpp>

#include "pgmlink/feature.h"
#include "pgmlink/field_of_view.h"

using namespace pgmlink;
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
    energy(ancestor, descendant1, descendant2);
    //BOOST_CHECK_EQUAL(e, 0);
}


BOOST_AUTO_TEST_CASE( SpatialDistanceToBorder )
{
    Traxel t1, t2, t3, t4, t5, t6;
    feature_array com(feature_array::difference_type(3));

    com[0] = 1;
    com[1] = 1;
    com[2] = 0;
    t1.features["com"] = com;

	com[0] = 1;
	com[1] = 10;
	com[2] = 0;
	t2.features["com"] = com;

	com[0] = 5;
	com[1] = 5;
	com[2] = 0;
	t3.features["com"] = com;

	com[0] = 8;
	com[1] = 8;
	com[2] = 0;
	t4.features["com"] = com;

	com[0] = 8;
	com[1] = 5;
	com[2] = 0;
	t5.features["com"] = com;

	com[0] = 0.5;
	com[1] = 7;
	com[2] = 0;
	t6.features["com"] = com;

	FieldOfView fov(0, 0, 0, 0, 1, 10, 10, 0); // tlow, xlow, ylow, zlow, tup, xup, yup, zup

	// absolute margin
	double cost = 100;
	double border_width = 2;
	boost::function<double(const Traxel&)> cost_fn = SpatialBorderAwareWeight(cost,
																			border_width,
																			false, //relative margin to border
																			fov);


	BOOST_CHECK_EQUAL(cost_fn(t1), 50.);
	BOOST_CHECK_EQUAL(cost_fn(t2), 0.);
	BOOST_CHECK_EQUAL(cost_fn(t3), 100.);
	BOOST_CHECK_EQUAL(cost_fn(t4), 100.);
	BOOST_CHECK_EQUAL(cost_fn(t5), 100.);
	BOOST_CHECK_EQUAL(cost_fn(t6), 25.);


	// relative margin
	cost = 100;
	border_width = 0.2;
	cost_fn = SpatialBorderAwareWeight(cost,
									border_width,
									true, //relative margin to border
									fov);


	BOOST_CHECK_EQUAL(cost_fn(t1), 50.);
	BOOST_CHECK_EQUAL(cost_fn(t2), 0.);
	BOOST_CHECK_EQUAL(cost_fn(t3), 100.);
	BOOST_CHECK_EQUAL(cost_fn(t4), 100.);
	BOOST_CHECK_EQUAL(cost_fn(t5), 100.);
	BOOST_CHECK_EQUAL(cost_fn(t6), 25.);
}

// EOF
