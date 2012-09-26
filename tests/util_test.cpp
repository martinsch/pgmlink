#define BOOST_TEST_MODULE util_test

#include <iostream>
#include <vector>
#include <utility>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "pgmlink/util.h"

using namespace Tracking;
using namespace std;
using namespace boost;

BOOST_AUTO_TEST_CASE( sort_indices_test ) {
  vector<int> v;
  v.push_back( 9 );
  v.push_back( 3 );
  v.push_back( 1 );
  v.push_back( 1 );
  v.push_back( 7 );
  vector<size_t> inds;
  
  indexsorter::sort_indices(v.begin(), v.end(), inds);

  BOOST_CHECK_EQUAL(v.size(), 5);  
  BOOST_CHECK_EQUAL(v[0], 9);
  BOOST_CHECK_EQUAL(v[1], 3);
  BOOST_CHECK_EQUAL(v[2], 1);
  BOOST_CHECK_EQUAL(v[3], 1);
  BOOST_CHECK_EQUAL(v[4], 7);

  BOOST_CHECK_EQUAL(inds.size(), 5);  
  BOOST_CHECK_EQUAL(inds[0], 2);
  BOOST_CHECK_EQUAL(inds[1], 3);
  BOOST_CHECK_EQUAL(inds[2], 1);
  BOOST_CHECK_EQUAL(inds[3], 4);
  BOOST_CHECK_EQUAL(inds[4], 0);
} 

BOOST_AUTO_TEST_CASE( reorder_test ) {
  vector<int> v;
  v.push_back( 9 );
  v.push_back( 3 );
  v.push_back( 1 );
  v.push_back( 1 );
  v.push_back( 7 );

  vector<size_t> inds;
  inds.push_back(2);
  inds.push_back(3);
  inds.push_back(1);
  inds.push_back(4);
  inds.push_back(0);

  indexsorter::reorder( v, inds );
  BOOST_CHECK_EQUAL( v.size(), 5 );
  BOOST_CHECK_EQUAL(v[0], 1);
  BOOST_CHECK_EQUAL(v[1], 1);
  BOOST_CHECK_EQUAL(v[2], 3);
  BOOST_CHECK_EQUAL(v[3], 7);
  BOOST_CHECK_EQUAL(v[4], 9);
}
