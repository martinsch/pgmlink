#define BOOST_TEST_MODULE outlier_detection_test

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/shared_ptr.hpp>

#include "pgmlink/outlier_detection.h"
#include "pgmlink/tracks.h"
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

BOOST_AUTO_TEST_CASE( Outlier_badness ) {
  // Set up a track with outlier in feature "feature_x"
  feature_type f1_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {9., 8.}
  };
  feature_type f2_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {5., 5.}
  };
  Track track1;
  Track track2;
  for(size_t i = 0; i < 7; i++) {
    feature_array x1(f1_array[i], f1_array[i]+2);
    Traxel traxel1;
    traxel1.features["feature_x"] = x1;
    track1.traxels_.push_back(traxel1);
    
    feature_array x2(f2_array[i], f2_array[i]+2);
    Traxel traxel2;
    traxel2.features["feature_x"] = x2;
    track2.traxels_.push_back(traxel2);
  }
  
  // Set up feature calculator of first and second order
  boost::shared_ptr<FeatureCalculator> Identity(new IdentityCalculator());
  boost::shared_ptr<FeatureCalculator> Difference(new VectorDifferenceCalculator());
  boost::shared_ptr<FeatureAggregator> OutlierBadness(new OutlierBadnessAggregator());
  // Set up the track feature extractor for the feature "feature_x"
  TrackFeatureExtractor FeatureXOutlierBadness(Identity, OutlierBadness, "feature_x", SINGLE);
  TrackFeatureExtractor FeatureXDiffOutlierBadness(Difference, OutlierBadness, "feature_x", PAIRWISE);

  // Calculate position and velocity outlier badness for both tracks
  feature_type posbadness1 = FeatureXOutlierBadness.extract_scalar(track1);
  feature_type posbadness2 = FeatureXOutlierBadness.extract_scalar(track2);
  feature_type vbadness1 =  FeatureXDiffOutlierBadness.extract_scalar(track1);
  feature_type vbadness2 =  FeatureXDiffOutlierBadness.extract_scalar(track2);

  BOOST_CHECK(posbadness1 > posbadness2);
  BOOST_CHECK(vbadness1 > vbadness2);
}





