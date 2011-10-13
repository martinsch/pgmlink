#define BOOST_TEST_MODULE graphical_model_test

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <vigra/multi_array.hxx>

#include "energy.h"
#include "graphical_model.h"
#include "traxels.h"
#include "event.h"

using namespace Tracking;
using namespace std;
using namespace boost;
using namespace vigra;

BOOST_AUTO_TEST_CASE( opengm_factor_table_ordering )
{
  OpengmMrf mrf;
  mrf.Space()->addDimension(2);
  mrf.Space()->addDimension(2);
  mrf.Space()->addDimension(2);
  size_t vi[] = {0,1,2};
  size_t vv[] = {1,2,3,4,5,6,7,8};
  OpengmMrf::ogmFactor f(*(mrf.Space()), vi, vi+3, vv, vv+8);
  BOOST_CHECK_EQUAL(f(0,0,0), 1);
  BOOST_CHECK_EQUAL(f(1,0,0), 2);
  BOOST_CHECK_EQUAL(f(0,1,0), 3);
  BOOST_CHECK_EQUAL(f(1,1,0), 4);
  BOOST_CHECK_EQUAL(f(0,0,1), 5);
  BOOST_CHECK_EQUAL(f(1,0,1), 6);
  BOOST_CHECK_EQUAL(f(0,1,1), 7);
  BOOST_CHECK_EQUAL(f(1,1,1), 8);
}
// EOF

