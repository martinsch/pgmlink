#define BOOST_TEST_MODULE traxels_test

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/tuple/tuple.hpp>

#include "traxels.h"
#include "field_of_view.h"

using namespace Tracking;
using namespace std;
using namespace boost;

BOOST_AUTO_TEST_CASE( Traxel_assignment_op )
{
  Traxel t1,t2;
  t1.Id = 1;
  t2.Id = 2;
  t2 = t1;
  BOOST_CHECK_EQUAL(t1.Id, t2.Id);
}

BOOST_AUTO_TEST_CASE( Traxel_XYZ )
{
  Traxel t;
  feature_array com(3);
  com[0] = 34; com[1] = 45; com[2] = 12.3;
  t.features["com"] = com;
  BOOST_CHECK_EQUAL(t.X(), com[0]);
  BOOST_CHECK_EQUAL(t.Y(), com[1]);
  BOOST_CHECK_EQUAL(t.Z(), com[2]);
}

BOOST_AUTO_TEST_CASE( Traxel_resolution_scaling )
{
  ComLocator* locator = new ComLocator();
  locator->x_scale = 0.5;
  locator->z_scale = 12.3;
  Traxel t;
  t.set_locator(locator);
  feature_array com(3);
  com[0] = 34; com[1] = 45; com[2] = 12.3;
  t.features["com"] = com;
  BOOST_CHECK_EQUAL(t.X(), 0.5 * com[0]);
  BOOST_CHECK_EQUAL(t.Y(), com[1]);
  BOOST_CHECK_EQUAL(t.Z(), 12.3 * com[2]);
}

BOOST_AUTO_TEST_CASE( Traxel_distance_to )
{
    // prepare mock objects
    Traxel from, to;
    feature_array com1(3);
    feature_array com2(3);

    com1[0] = 1;
    com1[1] = 2;
    com1[2] = 0;
    from.features["com"] = com1;

    com2[0] = 2;
    com2[1] = 3;
    com2[2] = 7;
    to.features["com"] = com2;

    double distance = from.distance_to(to);
    BOOST_CHECK_CLOSE( distance, 7.1414284285, 0.01);
 }

BOOST_AUTO_TEST_CASE( Traxel_angle )
{
    // prepare mock objects
    Traxel vertex, leg1, leg2;
    feature_array com1(3);
    feature_array com2(3);
    feature_array com3(3);

    com1[0] = 1;
    com1[1] = 2;
    com1[2] = 0;
    vertex.features["com"] = com1;

    com2[0] = 2;
    com2[1] = 3;
    com2[2] = 7;
    leg1.features["com"] = com2;

    com3[0] = 1;
    com3[1] = 0.5;
    com3[2] = 1.5;
    leg2.features["com"] = com3;

    double angle = vertex.angle(leg1, leg2);
    BOOST_CHECK_CLOSE( angle, 0.93466, 0.01);
    angle = vertex.angle(leg2, leg1);
    BOOST_CHECK_CLOSE( angle, 0.93466, 0.01);

    // and another one
    com1[0] = 1;
    com1[1] = 2;
    com1[2] = 0;
    vertex.features["com"] = com1;

    com2[0] = 2;
    com2[1] = 3;
    com2[2] = 0;
    leg1.features["com"] = com2;

    com3[0] = 3;
    com3[1] = 1;
    com3[2] = 0;
    leg2.features["com"] = com3;

    angle = vertex.angle(leg1, leg2);
    BOOST_CHECK_CLOSE( angle, 1.2490, 0.01);
    angle = vertex.angle(leg2, leg1);
    BOOST_CHECK_CLOSE( angle, 1.2490, 0.01);

    // check angle > 90°
    com1[0] = 1;
    com1[1] = 0;
    com1[2] = 0;
    vertex.features["com"] = com1;

    com2[0] = 2;
    com2[1] = 0;
    com2[2] = 0;
    leg1.features["com"] = com2;

    com3[0] = 0;
    com3[1] = 0;
    com3[2] = 0;
    leg2.features["com"] = com3;

    angle = vertex.angle(leg1, leg2);
    BOOST_CHECK_CLOSE( angle, 3.1415, 0.01);
    angle = vertex.angle(leg2, leg1);
    BOOST_CHECK_CLOSE( angle, 3.1415, 0.01);

}

BOOST_AUTO_TEST_CASE( Traxel_strict_weak_ordering )
{
    Traxel t1, t2;
    t1.Timestep = 1;
    t2.Timestep = 2;
    t1.Id = 1;
    t2.Id = 2;
    BOOST_CHECK(t1 < t2);
    BOOST_CHECK(!(t2 < t1));

    t1.Timestep = 2;
    t2.Timestep = 1;
    t1.Id = 1;
    t2.Id = 2;
    BOOST_CHECK(t2 < t1);
    BOOST_CHECK(!(t1 < t2));

    t1.Timestep = 1;
    t2.Timestep = 1;
    t1.Id = 1;
    t2.Id = 1;
    BOOST_CHECK(!(t1 < t2 || t2 < t1));
}

BOOST_AUTO_TEST_CASE( global_fun_nested_vec_from )
{
  Traxel t1,t21,t22,t4;
  t1.Timestep = 1;
  t1.Id = 0;
  t21.Timestep = 2;
  t21.Id = 0;
  t22.Timestep = 2;
  t22.Id = 1;
  t4.Timestep = 4;
  TraxelStore t;
  t.get<by_timestep>().insert(t1);
  t.get<by_timestep>().insert(t21);
  t.get<by_timestep>().insert(t22);
  t.get<by_timestep>().insert(t4);
  
  vector<vector<Traxel> > nested_vec;
  nested_vec = nested_vec_from(t);
  BOOST_REQUIRE_EQUAL(nested_vec.size(), 4);
  BOOST_CHECK_EQUAL(nested_vec[0].size(), 1);
  BOOST_CHECK_EQUAL(nested_vec[1].size(), 2);
  BOOST_CHECK_EQUAL(nested_vec[2].size(), 0);
  BOOST_CHECK_EQUAL(nested_vec[3].size(), 1);
}

BOOST_AUTO_TEST_CASE( global_fun_filter_by_fov )
{
  FieldOfView fov( 1,1,1,1, 2,2,2,2 );
  Traxel outside, inside, on_the_border, in_space_out_time, out_space_in_time;
  feature_array com1(3), com2(3), com3(3);
  com1[0] = 0.5; com1[1] = 2.5; com1[2] = 1;
  com2[0] = 1.5; com2[1] = 1.5; com2[2] = 1.5;
  com3[0] = 2;   com3[1] = 2;   com3[2] = 2;

  outside.features["com"] = com1;
  outside.Timestep = 3;
  outside.Id = 0;

  inside.features["com"] = com2;
  inside.Timestep = 2;
  inside.Id = 1;

  on_the_border.features["com"] = com3;
  on_the_border.Timestep = 1;
  on_the_border.Id = 2;

  in_space_out_time.features["com"] = com2;
  in_space_out_time.Timestep = 0;
  in_space_out_time.Id = 3;

  out_space_in_time.features["com"] = com1;
  out_space_in_time.Timestep = 2;
  out_space_in_time.Id = 4;

  Traxel traxels[] = {outside, inside, on_the_border, in_space_out_time, out_space_in_time};
  TraxelStore ts, ts_out;
  add(ts, traxels, traxels + (sizeof(traxels) / sizeof(Traxel)));
  BOOST_CHECK_EQUAL(ts.size(), 5);
  BOOST_CHECK_EQUAL(ts_out.size(), 0);  

  size_t n = filter_by_fov(ts, ts_out, fov);
  BOOST_CHECK_EQUAL(ts.size(), 5);
  BOOST_CHECK_EQUAL(n, 2);
  BOOST_CHECK_EQUAL(ts_out.size(), 2);
  BOOST_CHECK_EQUAL(ts_out.get<by_timeid>().count(tuple<int, unsigned int>(2,1)), 1);
  BOOST_CHECK_EQUAL(ts_out.get<by_timeid>().count(tuple<int, unsigned int>(1,2)), 1);
}
// EOF
