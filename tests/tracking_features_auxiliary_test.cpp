#define BOOST_TEST_MODULE outlier_detection_test

#include <cmath>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/shared_ptr.hpp>

#include "pgmlink/tracking_features.h"
#include "pgmlink/classifier_auxiliary.h"

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
  for(size_t i = 0; i < 3; i++) {
    feature_array y(setI_array[i], setI_array[i]+2);
    setI.push_back(y);
  }
  feature_type setII_array[2][1] = {{0}, {2}};
  feature_arrays setII;
  for(size_t i = 0; i < 2; i++) {
    feature_array y(setII_array[i], setII_array[i]+1);
    setII.push_back(y);
  }

  TotalDiffAggregator totaldiff;
  feature_array vectorI = totaldiff.vector_valued(setI);
  feature_type scalarI = totaldiff.scalar_valued(setI);

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

BOOST_AUTO_TEST_CASE( MinAggregator_test ) {
  // Set up some test data
  feature_type setI_array[3][2] = {
    {0., 1.},
    {4., 8.},
    {2., 3.}
  };
  feature_arrays setI;
  for(size_t i = 0; i < 3; i++) {
    feature_array y(setI_array[i], setI_array[i]+2);
    setI.push_back(y);
  }
  feature_type setII_array[2][1] = {{0}, {2}};
  feature_arrays setII;
  for(size_t i = 0; i < 2; i++) {
    feature_array y(setII_array[i], setII_array[i]+1);
    setII.push_back(y);
  }

  MinAggregator min;
  feature_array vectorI = min.vector_valued(setI);
  feature_type scalarI = min.scalar_valued(setI);

  feature_array vectorII = min.vector_valued(setII);
  feature_type scalarII = min.scalar_valued(setII);

  feature_array vectorIref(2); vectorIref[0] = 0.; vectorIref[1] = 1.;
  feature_type scalarIref = 0.;

  feature_array vectorIIref(1); vectorIIref[0] = 0.;
  feature_type scalarIIref = 0.;

  BOOST_CHECK_EQUAL(vectorIref[0], vectorI[0]);
  BOOST_CHECK_EQUAL(vectorIref[1], vectorI[1]);
  BOOST_CHECK_EQUAL(scalarIref, scalarI);
  
  BOOST_CHECK_EQUAL(vectorIIref[0], vectorII[0]);
  BOOST_CHECK_EQUAL(scalarIIref, scalarII);

}

BOOST_AUTO_TEST_CASE( MaxAggregator_test ) {
  // Set up some test data
  feature_type setI_array[3][2] = {
    {0., 1.},
    {4., 8.},
    {2., 3.}
  };
  feature_arrays setI;
  for(size_t i = 0; i < 3; i++) {
    feature_array y(setI_array[i], setI_array[i]+2);
    setI.push_back(y);
  }
  feature_type setII_array[2][1] = {{0}, {2}};
  feature_arrays setII;
  for(size_t i = 0; i < 2; i++) {
    feature_array y(setII_array[i], setII_array[i]+1);
    setII.push_back(y);
  }

  MaxAggregator max;
  feature_array vectorI = max.vector_valued(setI);
  feature_type scalarI = max.scalar_valued(setI);

  feature_array vectorII = max.vector_valued(setII);
  feature_type scalarII = max.scalar_valued(setII);

  feature_array vectorIref(2); vectorIref[0] = 4.; vectorIref[1] = 8.;
  feature_type scalarIref = 8.;

  feature_array vectorIIref(1); vectorIIref[0] = 2.;
  feature_type scalarIIref = 2.;

  BOOST_CHECK_EQUAL(vectorIref[0], vectorI[0]);
  BOOST_CHECK_EQUAL(vectorIref[1], vectorI[1]);
  BOOST_CHECK_EQUAL(scalarIref, scalarI);
  
  BOOST_CHECK_EQUAL(vectorIIref[0], vectorII[0]);
  BOOST_CHECK_EQUAL(scalarIIref, scalarII);

}

BOOST_AUTO_TEST_CASE( MeanAggregator_test ) {
  // Set up some test data
  feature_type setI_array[3][2] = {
    {0., 1.},
    {4., 8.},
    {2., 3.}
  };
  feature_arrays setI;
  for(size_t i = 0; i < 3; i++) {
    feature_array y(setI_array[i], setI_array[i]+2);
    setI.push_back(y);
  }
  feature_type setII_array[2][1] = {{0}, {2}};
  feature_arrays setII;
  for(size_t i = 0; i < 2; i++) {
    feature_array y(setII_array[i], setII_array[i]+1);
    setII.push_back(y);
  }

  MeanAggregator mean;
  feature_array vectorI = mean.vector_valued(setI);
  feature_type scalarI = mean.scalar_valued(setI);

  feature_array vectorII = mean.vector_valued(setII);
  feature_type scalarII = mean.scalar_valued(setII);

  feature_array vectorIref(2); vectorIref[0] = 2.; vectorIref[1] = 4.;
  feature_type scalarIref = 3.;

  feature_array vectorIIref(1); vectorIIref[0] = 1.;
  feature_type scalarIIref = 1.;

  BOOST_CHECK_EQUAL(vectorIref[0], vectorI[0]);
  BOOST_CHECK_EQUAL(vectorIref[1], vectorI[1]);
  BOOST_CHECK_EQUAL(scalarIref, scalarI);
  
  BOOST_CHECK_EQUAL(vectorIIref[0], vectorII[0]);
  BOOST_CHECK_EQUAL(scalarIIref, scalarII);

}

BOOST_AUTO_TEST_CASE( OutlierCountAggregator_test ) {
  // Set up some test data
  feature_type set11_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {9., 8.}
  };
  feature_type set12_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {4., 5.}
  };
  feature_type set21_array[7][1] = {
    {3.},
    {4.},
    {3.},
    {4.},
    {5.},
    {5.},
    {9.}
  };
  feature_type set22_array[7][1] = {
    {3.},
    {4.},
    {3.},
    {4.},
    {5.},
    {5.},
    {5.}
  };
  feature_arrays set11;
  feature_arrays set12;
  for(size_t i = 0; i < 7; i++) {
    feature_array y1(set11_array[i], set11_array[i]+2);
    set11.push_back(y1);
    feature_array y2(set12_array[i], set12_array[i]+2);
    set12.push_back(y2);
  }
  feature_arrays set21;
  feature_arrays set22;
  for(size_t i = 0; i < 7; i++) {
    feature_array y1(set21_array[i], set21_array[i]+1);
    set21.push_back(y1);
    feature_array y2(set22_array[i], set22_array[i]+1);
    set22.push_back(y2);
  }

  OutlierCountAggregator outliercount;

  feature_type scalar11 = outliercount.scalar_valued(set11);
  feature_type scalar12 = outliercount.scalar_valued(set12);
  feature_type scalar21 = outliercount.scalar_valued(set21);
  feature_type scalar22 = outliercount.scalar_valued(set22);

  feature_type scalar11ref = 1.;
  feature_type scalar12ref = 0.;
  feature_type scalar21ref = 1.;
  feature_type scalar22ref = 0.;

  BOOST_CHECK_EQUAL(scalar11ref, scalar11);
  BOOST_CHECK_EQUAL(scalar12ref, scalar12);
  BOOST_CHECK_EQUAL(scalar21ref, scalar21);
  BOOST_CHECK_EQUAL(scalar22ref, scalar22);

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

BOOST_AUTO_TEST_CASE( TrackingFeatureExtractor_outliercount ) {
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
    {4., 4.}
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
  Trackvector tracks;
  tracks.push_back(track1); tracks.push_back(track2);
  Tracking tracking;
  tracking.tracks_ = tracks;
  
  // Set up feature calculator and the feature aggregator
  boost::shared_ptr<FeatureCalculator> Identity(new IdentityCalculator());
  boost::shared_ptr<FeatureAggregator> OutlierCount(
    new OutlierCountAggregator()
  );
  boost::shared_ptr<FeatureAggregator> Sum(
    new SumAggregator()
  );

  // Set up the track feature extractor for the feature "feature_x"
  boost::shared_ptr<TrackFeatureExtractor> TrackOutlier(
    new TrackFeatureExtractor(Identity, OutlierCount, "feature_x", SINGLE)
  );

  // Set up the tracking feature extractor
  TrackingFeatureExtractor OutlierCountSum(TrackOutlier, Sum, SCALAR);

  std::cout << "Outlier count in track 1: ";
  std::cout << TrackOutlier->extract_scalar(track1) << std::endl;
  std::cout << "Outlier count in track 2: ";
  std::cout << TrackOutlier->extract_scalar(track2) << std::endl;

  BOOST_CHECK_EQUAL(OutlierCountSum.extract(tracking), 1);

}
