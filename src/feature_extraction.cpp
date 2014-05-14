// stl
#include <string>
#include <stdexcept>
#include <map>

// boost
#include <boost/shared_ptr.hpp>

// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/feature_extraction.h"

namespace pgmlink {

namespace feature_extraction {

////
//// class FeatureCalculator
////
FeatureCalculator::~FeatureCalculator() {

}


feature_array FeatureCalculator::calculate(const feature_array& /* f1 */) const {
  throw std::runtime_error("Not implemented yet!");
}


feature_array FeatureCalculator::calculate(const feature_array& /* f1 */, const feature_array& /* f2 */) const {
  throw std::runtime_error("Not implemented yet!");
}


feature_array FeatureCalculator::calculate(const feature_array& /* f1 */, const feature_array& /* f2 */, const feature_array& /* f3 */) const {
  throw std::runtime_error("Not implemented yet!");
}


const std::string& FeatureCalculator::name() const {
  return name_;
}


bool FeatureCalculator::operator==(const FeatureCalculator& other) {
  return this->name_ == other.name();
}


bool FeatureCalculator::operator!=(const FeatureCalculator& other) {
  return !(*this == other);
}


const std::string FeatureCalculator::name_ = "";



////
//// class PairwiseSquaredDistanceCalculator
////
ElementWiseSquaredDistanceCalculator::~ElementWiseSquaredDistanceCalculator() {

}


feature_array ElementWiseSquaredDistanceCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  LOG(logDEBUG4) << "ElementWiseSquaredDistanceCalculator::calculate(const feature_array& f1, const feature_array& f2) -- entered";
  assert(f1.size() == f2.size() && "feature array sizes must agree.");
  feature_array result = feature_array(f1.size());
  {
    feature_array::iterator res_it = result.begin();
    for (
             feature_array::const_iterator it1 = f1.begin(),
                 it2 = f2.begin();
             it1 != f1.end();
             ++it1, ++it2, ++res_it) {
      feature_type diff = *it1 - *it2;
      *res_it = diff * diff;
    }
  }
  LOG(logDEBUG4) << "ElementWiseSquaredDistanceCalculator::calculate(const feature_array& f1, const feature_array& f2) -- exit";
  return result;
}


const std::string& ElementWiseSquaredDistanceCalculator::name() const {
  return name_;
}


const std::string ElementWiseSquaredDistanceCalculator::name_ = "squared distance";


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




namespace helpers {

namespace {
std::map<std::string, boost::shared_ptr<FeatureCalculator> > define_features() {
  // put here all the available features:
  std::map<std::string, boost::shared_ptr<FeatureCalculator> > feature_map;
  
  boost::shared_ptr<FeatureCalculator> calc = // boost::shared_ptr<FeatureCalculator>();
      boost::shared_ptr<FeatureCalculator>(new ElementWiseSquaredDistanceCalculator);
  feature_map.insert(std::make_pair("ElementWiseSquaredDistance", calc));

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


} 



} /* namespace feature_extraction */

} /* namespace pgmlink */


