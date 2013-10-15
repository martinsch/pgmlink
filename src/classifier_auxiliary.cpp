// stl
#include <vector>
#include <string>
#include <numeric>
#include <cassert>
#include <stdexcept>

// boost
#include <boost/shared_ptr.hpp>

// vigra

// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/traxels.h"
#include "pgmlink/classifier_auxiliary.h"


namespace pgmlink {

////
//// class FeatureCalculator
////
const std::string FeatureCalculator::name = "";

const unsigned FeatureCalculator::length = 0;


FeatureCalculator::~FeatureCalculator() {

}


feature_array FeatureCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  throw std::runtime_error("FeatureCalculator \"" + name + "\" does not take two feature arrays");
  return feature_array();
}


feature_array FeatureCalculator::calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const {
  throw std::runtime_error("FeatureCalculator \"" + name + "\" does not take three feature arrays");
  return feature_array();
}


bool FeatureCalculator::operator==(const FeatureCalculator& other) {
  return name == other.name;
}


////
//// class DistanceCalculator
////
const std::string DistanceCalculator::name = "distance";

const unsigned DistanceCalculator::length = 1;


DistanceCalculator::~DistanceCalculator() {

}


feature_array DistanceCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  assert(f1.size() == f2.size());
  feature_array ret(length, 0.);
  feature_array::const_iterator f1_it = f1.begin();
  feature_array::const_iterator f2_it = f2.begin();
  for (; f1_it != f1.end(); ++f1_it, ++f2_it) {
    ret[0] += (*f1_it - *f2_it)*(*f1_it - *f2_it);
  }
  return ret;
}


////
//// class FeatureExtractor
////
FeatureExtractor::FeatureExtractor(boost::shared_ptr<FeatureCalculator> calculator, const std::vector<std::string>& feature_names)
    : calculator_(calculator), feature_names_(feature_names) {
  
}


feature_array FeatureExtractor::extract(const Traxel& t1, const Traxel& t2) const {
  feature_array features1, features2;
  for (std::vector<std::string>::const_iterator feature_name = feature_names_.begin();
       feature_name != feature_names_.end();
       ++feature_name) {
    FeatureMap::const_iterator f1 = t1.features.find(*feature_name);
    FeatureMap::const_iterator f2 = t2.features.find(*feature_name);
    assert(f1 != t1.features.end());
    assert(f2 != t2.features.end());
    features1.insert(features1.end(), f1->second.begin(), f1->second.end());
    features2.insert(features2.end(), f2->second.begin(), f2->second.end());
  }
  return calculator_->calculate(features1, features2);
}


feature_array FeatureExtractor::extract(const Traxel& t1, const Traxel& t2, const Traxel& t3) const {
  feature_array features1, features2, features3;
  for (std::vector<std::string>::const_iterator feature_name = feature_names_.begin();
       feature_name != feature_names_.end();
       ++feature_name) {
    FeatureMap::const_iterator f1 = t1.features.find(*feature_name);
    FeatureMap::const_iterator f2 = t2.features.find(*feature_name);
    FeatureMap::const_iterator f3 = t3.features.find(*feature_name);
    assert(f1 != t1.features.end());
    assert(f2 != t2.features.end());
    assert(f3 != t3.features.end());
    features1.insert(features1.end(), f1->second.begin(), f1->second.end());
    features2.insert(features2.end(), f2->second.begin(), f2->second.end());
    features2.insert(features3.end(), f3->second.begin(), f3->second.end());
  }
  return calculator_->calculate(features1, features2, features3);
}


boost::shared_ptr<FeatureCalculator> FeatureExtractor::calculator() {
  return calculator_;
}



} /* namespace pgmlink */
