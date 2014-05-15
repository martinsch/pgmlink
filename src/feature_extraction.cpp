// stl
#include <string>
#include <stdexcept>
#include <map>

// boost
#include <boost/shared_ptr.hpp>

// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/feature_extraction.h"
#include "pgmlink/feature_calculator.h"
#include "pgmlink/feature_calculator/base.h"

namespace pgmlink {

namespace feature_extraction {


// class FeatureCalculator;


////
//// class FeatureExtractor
////
FeatureExtractor::FeatureExtractor(boost::shared_ptr<FeatureCalculator> calculator, const std::string& feature_name)
    : calculator_(calculator), feature_name_(feature_name) {
  
}


FeatureExtractor::~FeatureExtractor() {

}


feature_array FeatureExtractor::extract(const Traxel& t1) const {
  LOG(logDEBUG4) << "FeatureExtractor::extract: feature " << feature_name_;
  FeatureMap::const_iterator f1 = t1.features.find(feature_name_);
  if ( f1 == t1.features.end() ) {
    throw std::runtime_error("Feature " + feature_name_ + " not present in traxel.");
  }
  feature_array features1 = f1->second;
  LOG(logDEBUG4) << "FeatureExtractor::extract: feature[0] =  " << features1[0];
  return calculator_->calculate(features1);
}


feature_array FeatureExtractor::extract(const Traxel& t1, const Traxel& t2) const {
  LOG(logDEBUG4) << "FeatureExtractor::extract: feature " << feature_name_;
  FeatureMap::const_iterator f1 = t1.features.find(feature_name_);
  FeatureMap::const_iterator f2 = t2.features.find(feature_name_);
  if ( f1 == t1.features.end() || f2 == t2.features.end() ) {
    throw std::runtime_error("Feature " + feature_name_ + " not present in traxel.");
  }
  feature_array features1 = f1->second;
  feature_array features2 = f2->second;
  return calculator_->calculate(features1, features2);
}


feature_array FeatureExtractor::extract(const Traxel& t1, const Traxel& t2, const Traxel& t3) const {
  FeatureMap::const_iterator f1 = t1.features.find(feature_name_);
  FeatureMap::const_iterator f2 = t2.features.find(feature_name_);
  FeatureMap::const_iterator f3 = t3.features.find(feature_name_);
  if ( f1 == t1.features.end() || f2 == t2.features.end() || f3 == t3.features.end() ) {
    throw std::runtime_error("Feature " + feature_name_ + " not present in traxel.");
  }
  feature_array features1 = f1->second;
  feature_array features2 = f2->second;
  feature_array features3 = f3->second;
  return calculator_->calculate(features1, features2, features3);
}


boost::shared_ptr<FeatureCalculator> FeatureExtractor::calculator() const {
  return calculator_;
}

std::string FeatureExtractor::name() const {
  return calculator_->name() + "<" + feature_name_ + ">";
}


////
//// class MultipleFeaturesExtraction
////
MultipleFeatureExtraction::CombinedFeatureMap MultipleFeatureExtraction::operator() ( const FeatureList& features,
                                                                                      const Traxel& trax ) const {
  CombinedFeatureMap res;
  for ( FeatureList::const_iterator outer_features = features.begin(); outer_features != features.end(); ++outer_features) {
    boost::shared_ptr<FeatureCalculator> calc = helpers::CalculatorLookup::extract_calculator(outer_features->first );
    for ( std::vector<std::string>::const_iterator inner_features = outer_features->second.begin();
          inner_features != outer_features->second.end();
          ++inner_features ) {
      FeatureExtractor extractor( calc, *inner_features );
      res[std::make_pair(outer_features->first, *inner_features)] = extractor.extract( trax );
    }
  }
  return res;
}


MultipleFeatureExtraction::CombinedFeatureMap MultipleFeatureExtraction::operator() ( const FeatureList& features,
                                                                                      const Traxel& trax1,
                                                                                      const Traxel& trax2 ) const {
  CombinedFeatureMap res;
  for ( FeatureList::const_iterator outer_features = features.begin(); outer_features != features.end(); ++outer_features) {
    boost::shared_ptr<FeatureCalculator> calc = helpers::CalculatorLookup::extract_calculator(outer_features->first );
    for ( std::vector<std::string>::const_iterator inner_features = outer_features->second.begin();
          inner_features != outer_features->second.end();
          ++inner_features ) {
      FeatureExtractor extractor( calc, *inner_features );
      res[std::make_pair(outer_features->first, *inner_features)] = extractor.extract( trax1, trax2 );
    }
  }
  return res;
}


MultipleFeatureExtraction::CombinedFeatureMap MultipleFeatureExtraction::operator() ( const FeatureList& features,
                                                                                      const Traxel& trax1,
                                                                                      const Traxel& trax2,
                                                                                      const Traxel& trax3 ) const {
  CombinedFeatureMap res;
  for ( FeatureList::const_iterator outer_features = features.begin(); outer_features != features.end(); ++outer_features) {
    boost::shared_ptr<FeatureCalculator> calc = helpers::CalculatorLookup::extract_calculator(outer_features->first );
    for ( std::vector<std::string>::const_iterator inner_features = outer_features->second.begin();
          inner_features != outer_features->second.end();
          ++inner_features ) {
      FeatureExtractor extractor( calc, *inner_features );
      res[std::make_pair(outer_features->first, *inner_features)] = extractor.extract( trax1, trax2, trax3 );
    }
  }
  return res;
}



namespace helpers {


////
//// function convenience_feature_extraction
////
MultipleFeatureExtraction::CombinedFeatureMap convenience_feature_extraction( const MultipleFeatureExtraction::FeatureList& features,
                                                                              const Traxel& trax ) {
  return MultipleFeatureExtraction()( features, trax );
}


MultipleFeatureExtraction::CombinedFeatureMap convenience_feature_extraction( const MultipleFeatureExtraction::FeatureList& features,
                                                                              const Traxel& trax1,
                                                                              const Traxel& trax2 ) {
  return MultipleFeatureExtraction()( features, trax1, trax2 );
}


MultipleFeatureExtraction::CombinedFeatureMap convenience_feature_extraction( const MultipleFeatureExtraction::FeatureList& features,
                                                                              const Traxel& trax1,
                                                                              const Traxel& trax2,
                                                                              const Traxel& trax3 ) {
  return MultipleFeatureExtraction()( features, trax1, trax2, trax3 );
}


} 



} /* namespace feature_extraction */

} /* namespace pgmlink */


