#define BOOST_TEST_MODULE feature_extraction_test


// stl
#include <string>
#include <stdexcept>
#include <iostream>
#include <utility>

// boost
#include <boost/shared_ptr.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator.h"
#include "pgmlink/feature_extraction.h"
#include "pgmlink/traxels.h"

namespace fe = pgmlink::feature_extraction;

BOOST_AUTO_TEST_CASE( FeatureCalculator_Test )
{
  fe::FeatureCalculator* calc = new fe::FeatureCalculator;
  fe::FeatureCalculator* squd = new fe::ElementWiseSquaredDistanceCalculator;
  pgmlink::feature_array f1( 2, 0. );
  pgmlink::feature_array f2( 2, 1. );
  pgmlink::feature_array f3 (2, 1. );

  BOOST_CHECK_THROW( calc->calculate( f1 ), std::runtime_error );
  BOOST_CHECK_THROW( calc->calculate( f1, f1 ), std::runtime_error );
  BOOST_CHECK_THROW( calc->calculate( f1, f1, f1 ), std::runtime_error );

  BOOST_CHECK_THROW( squd->calculate( f1 ), std::runtime_error );
  BOOST_CHECK_THROW( squd->calculate( f1, f1, f1 ), std::runtime_error );


  BOOST_CHECK_EQUAL( calc->name(), "" );
  BOOST_CHECK_EQUAL( squd->name(), "squared distance" );

  BOOST_CHECK( *calc == *calc );
  BOOST_CHECK( *calc != *squd );
  
  BOOST_CHECK( ! ( *calc == *squd ) );
  BOOST_CHECK( ! ( *calc != *calc ) );

  pgmlink::feature_array r1 = squd->calculate( f1, f2 );
  pgmlink::feature_array r2 = squd->calculate( f2, f1 );
  BOOST_CHECK_EQUAL_COLLECTIONS( r1.begin(), r1.end(), f3.begin(), f3.end() );
  BOOST_CHECK_EQUAL_COLLECTIONS( r2.begin(), r2.end(), f3.begin(), f3.end() );

  delete calc;
  delete squd;
}


BOOST_AUTO_TEST_CASE( FeatureExtractor_Test_2_Traxels )
{
  pgmlink::feature_array f1( 2, 2. );
  pgmlink::feature_array f2( 2, 1. );
  pgmlink::feature_array f3( 2, 1. );

  pgmlink::Traxel t1;
  pgmlink::Traxel t2;

  t1.features["some_feature"] = f1;
  t2.features["some_feature"] = f2;

  boost::shared_ptr<fe::FeatureCalculator> calc( new fe::ElementWiseSquaredDistanceCalculator );

  fe::FeatureExtractor e1( calc, "some_feature" );
  fe::FeatureExtractor e2( calc, "else_feature" );

  pgmlink::feature_array res = e1.extract( t1, t2 );
  BOOST_CHECK_EQUAL_COLLECTIONS( res.begin(), res.end(), f3.begin(), f3.end());

  BOOST_CHECK_THROW(e2.extract(t1, t2), std::runtime_error);
}


BOOST_AUTO_TEST_CASE( CalculatorLookup ) {
  boost::shared_ptr<fe::FeatureCalculator> calc = fe::helpers::CalculatorLookup::extract_calculator( "ElementWiseSquaredDistance" );
  BOOST_CHECK_EQUAL( calc->name(), "squared distance" );
  BOOST_CHECK_THROW( fe::helpers::CalculatorLookup::extract_calculator("SomeOtherCalculator"), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( MultipleFeatureExtraction ) {
  
  fe::MultipleFeatureExtraction ex;
  fe::MultipleFeatureExtraction::FeatureList unary_flist, pairwise_flist, ternary_flist;
  unary_flist["ElementWiseSquaredDistance"].push_back( "feat" );

  fe::MultipleFeatureExtraction::CombinedFeatureMap unary_result, pairwise_result, ternary_result;
  pairwise_result[std::make_pair("ElementWiseSquaredDistance", "feat")] = pgmlink::feature_array(2, 2.);
  
  pgmlink::Traxel t1;
  pgmlink::Traxel t2;
  pgmlink::Traxel t3;

  pgmlink::feature_array f1( 2, 1. );
  pgmlink::feature_array f2( 2, 3. );
  pgmlink::feature_array f3( 2, 4. );

  t1.features["feat"] = f1;
  t2.features["feat"] = f2;
  t3.features["feat"] = f3;

  // unary and ternary to be done!!!
  // fe::MultipleFeatureExtraction::CombinedFeatureMap unary_cmap = ex( unary_flist, t1 );
  fe::MultipleFeatureExtraction::CombinedFeatureMap pairwise_cmap = ex( pairwise_flist, t1, t2 );
  // fe::MultipleFeatureExtraction::CombinedFeatureMap ternary_cmap = ex( ternary_flist, t1, t2, t3 );

  
  // binary
  BOOST_REQUIRE( pairwise_cmap.size() == pairwise_flist.size() );
  for ( fe::MultipleFeatureExtraction::CombinedFeatureMap::const_iterator res = pairwise_cmap.begin(), comp = pairwise_result.begin();
        res != pairwise_cmap.end();
        ++res ) {
    BOOST_REQUIRE( res->first == comp->first );
    BOOST_CHECK_EQUAL_COLLECTIONS( res->second.begin(), res->second.end(), comp->second.begin(), comp->second.end() );
  }

  
  // convenience checks
  // binary
  fe::MultipleFeatureExtraction::CombinedFeatureMap pairwise_convenient = fe::helpers::convenience_feature_extraction( pairwise_flist, t1, t2 );
  BOOST_REQUIRE( pairwise_convenient.size() == pairwise_flist.size() );
  for ( fe::MultipleFeatureExtraction::CombinedFeatureMap::const_iterator res = pairwise_convenient.begin(), comp = pairwise_result.begin();
        res != pairwise_convenient.end();
        ++res ) {
    BOOST_REQUIRE( res->first == comp->first );
    BOOST_CHECK_EQUAL_COLLECTIONS( res->second.begin(), res->second.end(), comp->second.begin(), comp->second.end() );
  }
}

