#define BOOST_TEST_MODULE higher_order_features_auxiliary_test

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "pgmlink/log.h"
#include "pgmlink/higher_order_features.h"

using namespace pgmlink;
// using namespace boost;

BOOST_AUTO_TEST_CASE( TrackFeaturesIdentity_operator ) {
  LOG(logINFO) << "test case: TrackFeaturesIdentity_operator";
  feature_type fx_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {9., 8.}
  };
  feature_type fy_array[7][1] = {
    {13.},
    {14.},
    {13.},
    {14.},
    {15.},
    {15.},
    {15.}
  };
  Track track;
  for(size_t i = 0; i < 7; i++) {
    feature_array fx(fx_array[i], fx_array[i]+2);
    feature_array fy(fy_array[i], fy_array[i]+1);
    Traxel traxel;
    traxel.features["feature_x"] = fx;
    traxel.features["feature_y"] = fy;
    track.traxels_.push_back(traxel);
  }

  LOG(logINFO) << "  test \"TrackFeaturesIdentity\" with string as constructor argument";
  TrackFeaturesIdentity fx_identity("feature_x");
  feature_arrays feature_x = fx_identity(track);
  BOOST_CHECK_EQUAL(feature_x.size(), 7);
  BOOST_CHECK_EQUAL(feature_x[0].size(), 2);
  BOOST_CHECK_EQUAL(feature_x[0][0], 3.);
  BOOST_CHECK_EQUAL(feature_x[6][1], 8.);

  LOG(logINFO) << "  test \"TrackFeaturesIdentity\" with vector as constructor argument";
  std::vector<std::string> feature_names;
  feature_names.push_back("feature_x");
  feature_names.push_back("feature_y");
  TrackFeaturesIdentity fv_identity(feature_names);
  feature_arrays feature_v = fv_identity(track);
  BOOST_CHECK_EQUAL(feature_v.size(), 7);
  BOOST_CHECK_EQUAL(feature_v[0].size(), 3);
  BOOST_CHECK_EQUAL(feature_v[0][0],  3.);
  BOOST_CHECK_EQUAL(feature_v[0][2], 13.);
  BOOST_CHECK_EQUAL(feature_v[6][0],  9.);
  BOOST_CHECK_EQUAL(feature_v[6][2], 15.);
}

BOOST_AUTO_TEST_CASE( TrackFeaturesDiff_operator ) {
  LOG(logINFO) << "test case: TrackFeaturesDiff_operator";
  feature_type fx_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {9., 8.}
  };
  feature_type fy_array[7][1] = {
    {13.},
    {14.},
    {13.},
    {14.},
    {15.},
    {15.},
    {15.}
  };
  Track track;
  for(size_t i = 0; i < 7; i++) {
    feature_array fx(fx_array[i], fx_array[i]+2);
    feature_array fy(fy_array[i], fy_array[i]+1);
    Traxel traxel;
    traxel.features["feature_x"] = fx;
    traxel.features["feature_y"] = fy;
    track.traxels_.push_back(traxel);
  }

  LOG(logINFO) << "  test \"TrackFeaturesDiff\" with string as constructor argument";
  TrackFeaturesDiff fx_diff("feature_x");
  feature_arrays feature_x = fx_diff(track);
  BOOST_CHECK_EQUAL(feature_x.size(), 6);
  BOOST_CHECK_EQUAL(feature_x[0].size(), 2);
  BOOST_CHECK_EQUAL(feature_x[0][0], 1.);
  BOOST_CHECK_EQUAL(feature_x[5][1], 3.);

  LOG(logINFO) << "  test \"TrackFeaturesDiff\" with vector as constructor argument";
  std::vector<std::string> feature_names;
  feature_names.push_back("feature_x");
  feature_names.push_back("feature_y");
  TrackFeaturesDiff fv_diff(feature_names);
  feature_arrays feature_v = fv_diff(track);
  BOOST_CHECK_EQUAL(feature_v.size(), 6);
  BOOST_CHECK_EQUAL(feature_v[0].size(), 3);
  BOOST_CHECK_EQUAL(feature_v[0][0], 1.);
  BOOST_CHECK_EQUAL(feature_v[0][2], 1.);
  BOOST_CHECK_EQUAL(feature_v[5][0], 4.);
  BOOST_CHECK_EQUAL(feature_v[5][2], 0.);
}

BOOST_AUTO_TEST_CASE( TrackFeaturesCurvature_operator ) {
  LOG(logINFO) << "test case: TrackFeaturesCurvature_operator";
  feature_type fx_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {9., 8.}
  };
  feature_type fy_array[7][1] = {
    {13.},
    {14.},
    {13.},
    {14.},
    {15.},
    {15.},
    {15.}
  };
  Track track;
  for(size_t i = 0; i < 7; i++) {
    feature_array fx(fx_array[i], fx_array[i]+2);
    feature_array fy(fy_array[i], fy_array[i]+1);
    Traxel traxel;
    traxel.features["feature_x"] = fx;
    traxel.features["feature_y"] = fy;
    track.traxels_.push_back(traxel);
  }

  LOG(logINFO) << "  test \"TrackFeaturesCurvature\" with string as constructor argument";
  TrackFeaturesCurvature fx_curvature("feature_x");
  feature_arrays feature_x = fx_curvature(track);
  BOOST_CHECK_EQUAL(feature_x.size(), 5);
  BOOST_CHECK_EQUAL(feature_x[0].size(), 2);
  BOOST_CHECK_EQUAL(feature_x[0][0], -2.);
  BOOST_CHECK_EQUAL(feature_x[4][1],  2.);

  LOG(logINFO) << "  test \"TrackFeaturesCurvature\" with vector as constructor argument";
  std::vector<std::string> feature_names;
  feature_names.push_back("feature_x");
  feature_names.push_back("feature_y");
  TrackFeaturesCurvature fv_curvature(feature_names);
  feature_arrays feature_v = fv_curvature(track);
  BOOST_CHECK_EQUAL(feature_v.size(), 5);
  BOOST_CHECK_EQUAL(feature_v[0].size(), 3);
  BOOST_CHECK_EQUAL(feature_v[0][0], -2.);
  BOOST_CHECK_EQUAL(feature_v[0][2], -2.);
  BOOST_CHECK_EQUAL(feature_v[4][0],  4.);
  BOOST_CHECK_EQUAL(feature_v[4][2],  0.);
}

BOOST_AUTO_TEST_CASE( MVNOutlierCalculator_calculate ) {
  LOG(logINFO) << "test case: MVNOutlierCalculator_calculate";
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
  LOG(logINFO) << "  set up the outlier calculator";
  MVNOutlierCalculator mvnoutlier;
  LOG(logINFO) << "  test outliers";
  std::vector<size_t> outlier(mvnoutlier.calculate(x));
  BOOST_CHECK_EQUAL(outlier.size(), 1);
  BOOST_CHECK_EQUAL(outlier[0], 6);
}

BOOST_AUTO_TEST_CASE( OutlierCountAggregator_calculate ) {
  LOG(logINFO) << "test case: OutlierCountAggregator_calculate";
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
  LOG(logINFO) << "  Set up the outlier count aggregator";
  OutlierCountAggregator outliercount;
  LOG(logINFO) << "  Calculate the outlier count";
  size_t outlier = outliercount(x);

  BOOST_CHECK_EQUAL(outlier, 1);
}

BOOST_AUTO_TEST_CASE( OutlierBadnessAggregator_calculate ) {
  LOG(logINFO) << "test case: OutlierBadnessAggregator_calculate";
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
  LOG(logINFO) << "  Set up the outlier badness aggregator";
  OutlierBadnessAggregator outlierbadness;
  LOG(logINFO) << "  Calculate the outlier badness";
  feature_type value = outlierbadness(x);
  BOOST_CHECK( value > 3 );
}
