#define BOOST_TEST_MODULE outlier_detection_test

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/shared_ptr.hpp>

#include "pgmlink/tracking_features.h"
#include "pgmlink/traxels.h"

using namespace pgmlink;
using namespace boost;

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
  MVNOutlierCalculator mvnoutlier;
  std::vector<size_t> outlier(mvnoutlier.calculate(x));

  BOOST_CHECK_EQUAL(outlier.size(), 1);
  BOOST_CHECK_EQUAL(outlier[0], 6);
}

/*
BOOST_AUTO_TEST_CASE( Outlier_calculate ) {
  // Set up a track with outlier in feature "feature_x"
  feature_type x_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {9., 8.}
  };
  Track track;
  for(size_t i = 0; i < 7; i++) {
    feature_array y(x_array[i], x_array[i]+2);
    Traxel traxel;
    traxel.features["feature_x"] = y;
    track.traxels_.push_back(traxel);
  }
  
  // Set up feature calculator of first and second order
  boost::shared_ptr<FeatureCalculator> Identity(new IdentityCalculator());
  boost::shared_ptr<FeatureCalculator> Difference(new VectorDifferenceCalculator());
  // Set up the track feature extractors for the feature "feature_x"
  boost::shared_ptr<TrackFeatureExtractor> FeatureX(new TrackFeatureExtractor(Identity, "feature_x", SINGLE));
  boost::shared_ptr<TrackFeatureExtractor> FeatureXDiff(new TrackFeatureExtractor(Difference, "feature_x", PAIRWISE));
  
  // Set up the outlier detection
  Outlier OutlierX(FeatureX);
  Outlier OutlierXDiff(FeatureXDiff);

  // Calculate the outliers
  std::vector<size_t> outlierx = OutlierX.calculate(track);
  std::vector<size_t> outlierxdiff = OutlierXDiff.calculate(track);

  BOOST_CHECK_EQUAL(outlierx.size(), 1);
  BOOST_CHECK_EQUAL(outlierx[0], 6);
  BOOST_CHECK_EQUAL(outlierxdiff.size(), 1);
  BOOST_CHECK_EQUAL(outlierxdiff[0], 5);
}
*/

BOOST_AUTO_TEST_CASE( OutlierCount_test ) {
  // Set up a track with outlier in feature "feature_x"
  feature_type x_array1[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {9., 8.}
  };
  Track track1;
  for(size_t i = 0; i < 7; i++) {
    feature_array y(x_array1[i], x_array1[i]+2);
    Traxel traxel;
    traxel.features["feature_x"] = y;
    track1.traxels_.push_back(traxel);
  }

  feature_type x_array2[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {4., 5.}
  };
  Track track2;
  for(size_t i = 0; i < 7; i++) {
    feature_array y(x_array2[i], x_array2[i]+2);
    Traxel traxel;
    traxel.features["feature_x"] = y;
    track2.traxels_.push_back(traxel);
  }
  Tracking tracking;
  tracking.tracks_.push_back(track1);
  tracking.tracks_.push_back(track2);

  OutlierCount outliercount("feature_x");
  BOOST_CHECK_EQUAL(outliercount(track1), 1);
  BOOST_CHECK_EQUAL(outliercount(track2), 0);
  BOOST_CHECK_EQUAL(outliercount(tracking), 1);
}
