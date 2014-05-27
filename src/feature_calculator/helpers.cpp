// stl
#include <map>
#include <string>

// pgmlink
#include "pgmlink/feature_calculator.h"
#include "pgmlink/feature_calculator/helpers.h"


namespace pgmlink {

namespace feature_extraction {



namespace helpers {


namespace {
std::map<std::string, boost::shared_ptr<FeatureCalculator> > define_features() {
  // put here all the available features:
  std::map<std::string, boost::shared_ptr<FeatureCalculator> > feature_map;
  
  boost::shared_ptr<FeatureCalculator> calc( new ElementWiseSquaredDistanceCalculator );
  feature_map.insert(std::make_pair( "ElementWiseSquaredDistance", calc) );

  calc = boost::shared_ptr<FeatureCalculator>( new IdentityCalculator );
  feature_map.insert( std::make_pair( "Identity", calc ) );

  calc = boost::shared_ptr<FeatureCalculator>( new AbsoluteDifferenceCalculator );
  feature_map.insert( std::make_pair( "AbsoluteDifference", calc ) );

  calc = boost::shared_ptr<FeatureCalculator>( new SquareRootSquaredDifferenceCalculator );
  feature_map.insert( std::make_pair( "SquareRooteSquaredDifference", calc ) );

  calc = boost::shared_ptr<FeatureCalculator>( new RatioCalculator );
  feature_map.insert( std::make_pair( "Ratio", calc ) );

  calc = boost::shared_ptr<FeatureCalculator>( new AsymmetricRatioCalculator );
  feature_map.insert( std::make_pair( "AsymmetricRatio", calc ) );

  calc = boost::shared_ptr<FeatureCalculator>( new SquaredDifferenceCalculator );
  feature_map.insert( std::make_pair( "SquaredDifference", calc ) );

  return feature_map;
}


} /* namespace */

////
//// class CalculatorLookup
////
boost::shared_ptr<FeatureCalculator> CalculatorLookup::extract_calculator(const std::string& name) {
  std::map<std::string, boost::shared_ptr<FeatureCalculator> >::const_iterator res = calculator_map_.find(name);
  if (res == calculator_map_.end()) {
    throw std::runtime_error("Calculator " + name + " not available!");
  }
  return res->second;
}


const std::map<std::string, boost::shared_ptr<FeatureCalculator> > CalculatorLookup::calculator_map_ = define_features();


} // namespace helpers

} // namespace feature_extraction

} // namespace pgmlink
