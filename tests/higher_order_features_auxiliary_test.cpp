#define BOOST_TEST_MODULE higher_order_features_auxiliary_test

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "pgmlink/log.h"
#include "pgmlink/higher_order_features_auxiliary.h"

using namespace pgmlink;
// using namespace boost;

BOOST_AUTO_TEST_CASE( MVNOutlierCalculator_calculate ) {
  // Set up the test data with outlier x6
  feature_type x_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {9., 8.}
  };
  feature_arrays x;
  for(size_t i = 0; i < 7; i++) {
    feature_array y(x_array[i], x_array[i]+2);
    x.push_back(y);
  }

  // Create outlier detection
  LOG(logINFO) << "Set up the outlier calculator";
  MVNOutlierCalculator mvnoutlier;
  LOG(logINFO) << "Calculate the outliers";
  std::vector<size_t> outlier(mvnoutlier.calculate(x));

  BOOST_CHECK_EQUAL(outlier.size(), 1);
  BOOST_CHECK_EQUAL(outlier[0], 6);
}

BOOST_AUTO_TEST_CASE( OutlierCountAggregator_calculate ) {
  // Set up the test data with outlier x6
  feature_type x_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {9., 8.}
  };
  feature_arrays x;
  for(size_t i = 0; i < 7; i++) {
    feature_array y(x_array[i], x_array[i]+2);
    x.push_back(y);
  }

  // Create outlier detection
  LOG(logINFO) << "Set up the outlier count aggregator";
  OutlierCountAggregator outliercount;
  LOG(logINFO) << "Calculate the outlier count";
  size_t outlier = outliercount(x);

  BOOST_CHECK_EQUAL(outlier, 1);
}
