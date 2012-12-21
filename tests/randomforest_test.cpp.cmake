#define BOOST_TEST_MODULE randomforest_test

#include <vector>
#include <iostream>

#include <boost/test/unit_test.hpp>
//#include <boost/test/floating_point_comparison.hpp>
//#include <boost/shared_ptr.hpp>

#include "pgmlink/randomforest.h"


BOOST_AUTO_TEST_CASE( checkCreateFeatureVector )
{
   // Set up one Traxel with some data
    pgmlink::Traxel tr;
    tr.Id = 1;
    pgmlink::feature_array feat1;
    pgmlink::feature_array feat2;
    pgmlink::feature_array feat3;
    pgmlink::feature_array feat4;

    feat1.push_back(1.5);
    for(int i = 0; i < 5; i++)
      feat2.push_back(float(i)*5.);
    for(int i = 0; i < 10; i++)
      feat3.push_back(float(i)*5.);
    for(int i = 0; i < 1000; i++)
      feat4.push_back(float(i)*5.);

    tr.features["feat1"] = feat1;
    tr.features["feat2"] = feat2;
    tr.features["feat3"] = feat2;
    tr.features["feat4"] = feat4;

    // 4 features were added
    BOOST_CHECK_EQUAL( tr.features.size(), 4 );

    {
        std::vector<std::string> sel;
        sel.push_back("feat1");
        sel.push_back("feat2");
        sel.push_back("feat4");
        vigra::MultiArray<2,float> fta = pgmlink::RF::createFeatureVector(tr,sel);

        BOOST_CHECK_EQUAL( fta.shape(0), 1 );
        BOOST_CHECK_EQUAL( fta.shape(1), 1006 );

        BOOST_CHECK_EQUAL( fta(0,0), (float)1.5 );

        BOOST_CHECK_EQUAL( fta(0,1), (float)0. );
        BOOST_CHECK_EQUAL( fta(0,2), (float)5. );
        BOOST_CHECK_EQUAL( fta(0,5), (float)20. );

        BOOST_CHECK_EQUAL( fta(0,1005), (float)4995. );
    }
    {
        std::vector<std::string> sel;
        sel.push_back("feat2");
        sel.push_back("feat1");
        vigra::MultiArray<2,float> fta = pgmlink::RF::createFeatureVector(tr,sel);

        BOOST_CHECK_EQUAL( fta.shape(0), 1 );
        BOOST_CHECK_EQUAL( fta.shape(1), 6 );

        BOOST_CHECK_EQUAL( fta(0,0), (float)0. );
        BOOST_CHECK_EQUAL( fta(0,1), (float)5. );
        BOOST_CHECK_EQUAL( fta(0,4), (float)20. );

        BOOST_CHECK_EQUAL( fta(0,5), (float)1.5 );
    }
    {

        std::vector<std::string> sel;
        sel.push_back("feat_invalid");
        vigra::MultiArray<2,float> fta = pgmlink::RF::createFeatureVector(tr,sel);

        BOOST_CHECK_EQUAL( fta.shape(0), 1 );
        BOOST_CHECK_EQUAL( fta.shape(1), 0 );
    }
    {
        std::vector<std::string> sel;
        sel.push_back("feat2");
        sel.push_back("feat_invalid");
        sel.push_back("feat1");
        vigra::MultiArray<2,float> fta = pgmlink::RF::createFeatureVector(tr,sel);

        BOOST_CHECK_EQUAL( fta.shape(0), 1 );
        BOOST_CHECK_EQUAL( fta.shape(1), 6 );

        BOOST_CHECK_EQUAL( fta(0,0), (float)0. );
        BOOST_CHECK_EQUAL( fta(0,1), (float)5. );
        BOOST_CHECK_EQUAL( fta(0,4), (float)20. );

        BOOST_CHECK_EQUAL( fta(0,5), (float)1.5 );
    }

}


BOOST_AUTO_TEST_CASE( checkGetRandomForest )
{
    // Read a Random Forest from the file "xorforest.h5"
    vigra::RandomForest<pgmlink::RF::RF_LABEL_TYPE> rf ( pgmlink::RF::getRandomForest("@PROJECT_SOURCE_DIR@/tests/xorforest.h5"));

    // Only two classes
    BOOST_CHECK_EQUAL( rf.class_count(), 2 );

    // Two features
    BOOST_CHECK_EQUAL( rf.feature_count(), 2 );

    // XOR test data
    vigra::MultiArray<2,float> test (pgmlink::RF::matrix_shape(4,2));
    test(0,0) = 0; test(0,1) = 0;
    test(1,0) = 1; test(1,1) = 0;
    test(2,0) = 0; test(2,1) = 1;
    test(3,0) = 1; test(3,1) = 1;

    vigra::MultiArray<2,pgmlink::RF::RF_LABEL_TYPE> pred (pgmlink::RF::matrix_shape(4,1));

    rf.predictLabels(test,pred);

    // test XOR behaviour
    BOOST_CHECK_EQUAL( pred(0,0), 0 );
    BOOST_CHECK_EQUAL( pred(1,0), 1 );
    BOOST_CHECK_EQUAL( pred(2,0), 1 );
    BOOST_CHECK_EQUAL( pred(3,0), 0 );
}


BOOST_AUTO_TEST_CASE( checkGetProbabilities )
{
    // Read a Random Forest from the file "xorforest.h5"
    vigra::RandomForest<pgmlink::RF::RF_LABEL_TYPE> rf ( pgmlink::RF::getRandomForest("@PROJECT_SOURCE_DIR@/tests/xorforest.h5"));

    // Only two classes
    BOOST_CHECK_EQUAL( rf.class_count(), 2 );

    // Two features
    BOOST_CHECK_EQUAL( rf.feature_count(), 2 );

    // XOR test data
    vigra::MultiArray<2,float> test (pgmlink::RF::matrix_shape(4,2));
    test(0,0) = 0; test(0,1) = 0;
    test(1,0) = 1; test(1,1) = 0;
    test(2,0) = 0; test(2,1) = 1;
    test(3,0) = 1; test(3,1) = 1;

    vigra::MultiArray<2,double> pred = pgmlink::RF::getProbabilities(test,rf);

    // test xor behaviour
    BOOST_CHECK( pred(0,0) > 0.5 );
    BOOST_CHECK( pred(1,0) < 0.5 );
    BOOST_CHECK( pred(2,0) < 0.5 );
    BOOST_CHECK( pred(3,0) > 0.5 );
}


BOOST_AUTO_TEST_CASE( checkPredictTracklets)
{
    pgmlink::feature_array f;
    pgmlink::Traxels ts;

    // Create 4 tracklets with 2 features each
    pgmlink::Traxel t1;
    t1.Id = 1;
    f.push_back(0.);
    t1.features["first"] = f;
    t1.features["second"] = f;
    ts[1] = t1;

    pgmlink::Traxel t2;
    t2.Id = 2;
    f[0] = 0;
    t2.features["first"] = f;
    f[0] = 1;
    t2.features["second"] = f;
    ts[2] = t2;

    pgmlink::Traxel t3;
    t3.Id = 3;
    f[0] = 1;
    t3.features["first"] = f;
    f[0] = 0;
    t3.features["second"] = f;
    ts[3] = t3;

    pgmlink::Traxel t4;
    t4.Id = 4;
    f[0] = 1;
    t4.features["first"] = f;
    t4.features["second"] = f;
    ts[4] = t4;


    // Read a Random Forest from the file "xorforest.h5"
    vigra::RandomForest<pgmlink::RF::RF_LABEL_TYPE> rf ( pgmlink::RF::getRandomForest("@PROJECT_SOURCE_DIR@/tests/xorforest.h5"));

    // Only two classes
    BOOST_CHECK_EQUAL( rf.class_count(), 2 );

    // Two features
    BOOST_CHECK_EQUAL( rf.feature_count(), 2 );


    std::vector<std::string> sel;
    sel.push_back("first"); sel.push_back("second");

    pgmlink::RF::predictTracklets(ts,rf,sel,1,"prediction");

    // check if tracklets have the desired feature
    BOOST_CHECK( ts[1].features.find("prediction") != ts[1].features.end());
    BOOST_CHECK( ts[2].features.find("prediction") != ts[2].features.end());
    BOOST_CHECK( ts[3].features.find("prediction") != ts[3].features.end());
    BOOST_CHECK( ts[4].features.find("prediction") != ts[4].features.end());

    // test XOR behaviour features["prediction"][0] contains probability
    // that features "first" and "second" are different.
    BOOST_CHECK( ts[1].features["prediction"][0] < 0.5 );
    BOOST_CHECK( ts[2].features["prediction"][0] > 0.5 );
    BOOST_CHECK( ts[3].features["prediction"][0] > 0.5 );
    BOOST_CHECK( ts[4].features["prediction"][0] < 0.5 );
}


BOOST_AUTO_TEST_CASE( checkLoadTracklets )
{
    // load the sample hdf5 file
    std::string filename = "@PROJECT_SOURCE_DIR@/tests/loadtest.h5";

    pgmlink::Traxels ts = pgmlink::RF::loadTracklets(filename);

    BOOST_CHECK( ts.size() == 5 );

    // check if every valid tracklet was found
    BOOST_CHECK( ts.find(1) != ts.end() );
    BOOST_CHECK( ts.find(2) != ts.end() );
    BOOST_CHECK( ts.find(3) != ts.end() );
    BOOST_CHECK( ts.find(4) != ts.end() );
    BOOST_CHECK( ts.find(5) == ts.end() );
    BOOST_CHECK( ts.find(6) != ts.end() );

    // check if ID's are correct
    BOOST_CHECK( ts[1].Id == 1 );
    BOOST_CHECK( ts[2].Id == 2 );
    BOOST_CHECK( ts[3].Id == 3 );
    BOOST_CHECK( ts[4].Id == 4 );
    BOOST_CHECK( ts[6].Id == 6 );


    // check some feature values
    BOOST_CHECK( ts[1].features["volume"][0] == 0 );
    BOOST_CHECK( ts[2].features["position"][1] == 1 );
    BOOST_CHECK( ts[3].features["com"][2] == 2 );
    BOOST_CHECK( ts[4].features["bbox"][3] == 3 );
    BOOST_CHECK( ts[6].features["intensity"][0] == 0 );
    BOOST_CHECK( ts[1].features["intminmax"][7] == 7 );
    BOOST_CHECK( ts[2].features["pair"][2] == 2 );
    BOOST_CHECK( ts[3].features["sgf"][4] == 4 );
    BOOST_CHECK( ts[4].features["pc"][1] == 1 );
}
