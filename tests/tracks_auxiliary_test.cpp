#define BOOST_TEST_MODULE outlier_detection_test

#include <cmath>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/shared_ptr.hpp>

#include "pgmlink/track_features_auxiliary.h"

using namespace pgmlink;
using namespace boost;

BOOST_AUTO_TEST_CASE( TotalDiffAggregator_test ) {
  // Set up some test data
  feature_type setI_array[3][2] = {
    {0., 1.},
    {4., 8.},
    {2., 3.}
  };
  feature_arrays setI;
  for(size_t i = 0; i < 7; i++) {
    feature_array y(setI_array[i], setI_array[i]+2);
    setI.push_back(y);
  }
  feature_type setII_array[2][1] = {{0}, {2}};
  feature_arrays setII;
  for(size_t i = 0; i < 7; i++) {
    feature_array y(setII_array[i], setII_array[i]+1);
    setII.push_back(y);
  }

  TotalDiffAggregator totaldiff;
  feature_array vectorI = totaldiff.vector_valued(setI);
  feature_type scalarI = totaldiff.scalar_valued(setII);

  feature_array vectorII = totaldiff.vector_valued(setII);
  feature_type scalarII = totaldiff.scalar_valued(setII);

  feature_array vectorIref(2); vectorIref[0] = 2; vectorIref[1] = 2;
  feature_type scalarIref = 2 * sqrt(2);

  feature_array vectorIIref(1); vectorIIref[0] = 2;
  feature_type scalarIIref = 2;

  BOOST_CHECK_EQUAL(vectorIref[0], vectorI[0]);
  BOOST_CHECK_EQUAL(vectorIref[1], vectorI[1]);
  BOOST_CHECK_EQUAL(scalarIref, scalarI);
  
  BOOST_CHECK_EQUAL(vectorIIref[0], vectorII[0]);
  BOOST_CHECK_EQUAL(scalarIIref, scalarII);

}
