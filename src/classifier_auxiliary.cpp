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
const std::string FeatureCalculator::name_ = "";

const unsigned FeatureCalculator::length = 0;


FeatureCalculator::~FeatureCalculator() {

}


feature_array FeatureCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  throw std::runtime_error("FeatureCalculator \"" + name() + "\" does not take two feature arrays");
  return feature_array();
}


feature_array FeatureCalculator::calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const {
  throw std::runtime_error("FeatureCalculator \"" + name() + "\" does not take three feature arrays");
  return feature_array();
}


const std::string& FeatureCalculator::name() const {
  return FeatureCalculator::name_;
}


bool FeatureCalculator::operator==(const FeatureCalculator& other) {
  return name() == other.name();
}


////
//// class SquaredDifferenceCalculator
////
const std::string SquaredDifferenceCalculator::name_ = "SquaredDiff";

const unsigned SquaredDifferenceCalculator::length = 1;


SquaredDifferenceCalculator::~SquaredDifferenceCalculator() {

}


feature_array SquaredDifferenceCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  assert(f1.size() == f2.size());
  feature_array ret(length, 0.);
  feature_array::const_iterator f1_it = f1.begin();
  feature_array::const_iterator f2_it = f2.begin();
  for (; f1_it != f1.end(); ++f1_it, ++f2_it) {
    ret[0] += (*f1_it - *f2_it)*(*f1_it - *f2_it);
  }
  return ret;
}


const std::string& SquaredDifferenceCalculator::name() const {
  return SquaredDifferenceCalculator::name_;
}


////
//// class AbsoluteDifferenceCalculator
////
const std::string AbsoluteDifferenceCalculator::name_ = "AbsDiff";

const unsigned AbsoluteDifferenceCalculator::length = 1;


AbsoluteDifferenceCalculator::~AbsoluteDifferenceCalculator() {

}


feature_array AbsoluteDifferenceCalculator::calculate(const feature_array&f1, const feature_array& f2) const {
  assert(f1.size() == f2.size());
  feature_array ret(length, 0.);
  feature_array::const_iterator f1_it = f1.begin();
  feature_array::const_iterator f2_it = f2.begin();
  for (; f1_it != f1.end(); ++f1_it, ++f2_it) {
    float res = *f1_it - *f2_it;
    ret[0] += res > 0 ? res : -res;
  }
  return ret;
}


const std::string& AbsoluteDifferenceCalculator::name() const {
  return AbsoluteDifferenceCalculator::name_;
}


////
//// class RatioCalculator
////
const std::string RatioCalculator::name_ = "Ratio";

const unsigned RatioCalculator::length = 1;


RatioCalculator::~RatioCalculator() {

}


feature_array RatioCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  assert(f1.size() == f2.size());
  feature_array ret(length, 0.);
  // keep ratio <= 1
  // no zero check, as we do not have empty regions
  assert(f1[0] > 0);
  assert(f2[0] > 0);
  if (f1[0] < f2[0]) {
    ret[0] = f1[0]/f2[0];
  } else {
    ret[0] = f2[0]/f1[0];
  }
  return ret;
}


feature_array RatioCalculator::calculate(const feature_array&, const feature_array& f2, const feature_array& f3) const {
  return calculate(f2, f3);
}


const std::string& RatioCalculator::name() const {
  return RatioCalculator::name_;
}



////
//// class IntensityRatioCalculator
////
/* const std::string IntensityRatioCalculator::name_ = "intensity_ratio";

const unsigned IntensityRatioCalculator::length = 1;


IntensityRatioCalculator::~IntensityRatioCalculator() {

}


feature_array IntensityRatioCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  assert(f1.size() == f2.size());
  feature_array ret(length, 0.);
  // keep ratio <= 1
  // no zero check, as we do not have empty regions
  assert(f1[0] > 0);
  assert(f2[0] > 0);
  if (f1[0] > f2[0]) {
    ret[0] = f2[0]/f1[0];
  } else {
    ret[0] = f1[0]/f2[0];
  }
  return ret;
}


feature_array IntensityRatioCalculator::calculate(const feature_array&, const feature_array& f2, const feature_array& f3) const {
  return calculate(f2, f3);
}


const std::string& IntensityRatioCalculator::name() const {
  return IntensityRatioCalculator::name_;
} */


////
//// class ChildrenMeanParentIntensityRatioCalculator
////
const std::string ChildrenMeanParentIntensityRatioCalculator::name_ = "children_parent_intensity_ratio";

const unsigned ChildrenMeanParentIntensityRatioCalculator::length = 1;


ChildrenMeanParentIntensityRatioCalculator::~ChildrenMeanParentIntensityRatioCalculator() {

}


feature_array ChildrenMeanParentIntensityRatioCalculator::calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const {
  feature_array ret(length, 0.);
  // mean of children intensities divided by parent intensity
  ret[0] = 0.5*(f2[0] + f3[0])/f1[0];
  return ret;
}


const std::string& ChildrenMeanParentIntensityRatioCalculator::name() const {
  return ChildrenMeanParentIntensityRatioCalculator::name_;
}



////
//// class FeatureExtractor
////
FeatureExtractor::FeatureExtractor(boost::shared_ptr<FeatureCalculator> calculator, const std::string& feature_name)
    : calculator_(calculator), feature_name_(feature_name) {
  
}


feature_array FeatureExtractor::extract(const Traxel& t1, const Traxel& t2) const {
  FeatureMap::const_iterator f1 = t1.features.find(feature_name_);
  FeatureMap::const_iterator f2 = t2.features.find(feature_name_);
  assert(f1 != t1.features.end());
  assert(f2 != t2.features.end());
  feature_array features1 = f1->second;
  feature_array features2 = f2->second;
  return calculator_->calculate(features1, features2);
}


feature_array FeatureExtractor::extract(const Traxel& t1, const Traxel& t2, const Traxel& t3) const {
  FeatureMap::const_iterator f1 = t1.features.find(feature_name_);
  FeatureMap::const_iterator f2 = t2.features.find(feature_name_);
  FeatureMap::const_iterator f3 = t3.features.find(feature_name_);
  assert(f1 != t1.features.end());
  assert(f2 != t2.features.end());
  assert(f3 != t3.features.end());
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

namespace {
std::map<std::string, boost::shared_ptr<FeatureCalculator> > define_features() {
  // put here all the available features:
  std::map<std::string, boost::shared_ptr<FeatureCalculator> > feature_map;
  
  boost::shared_ptr<FeatureCalculator> calc =
      boost::shared_ptr<FeatureCalculator>(new SquaredDifferenceCalculator);
  feature_map.insert(std::make_pair("SquaredDiff", calc));

  calc = boost::shared_ptr<FeatureCalculator>(new RatioCalculator);
  feature_map.insert(std::make_pair("Ratio", calc));

  calc = boost::shared_ptr<FeatureCalculator>(new AbsoluteDifferenceCalculator);
  feature_map.insert(std::make_pair("AbsDiff", calc));

  return feature_map;
}
} /* namespace */

////
//// class AvailableCalculators
////
std::map<std::string, boost::shared_ptr<FeatureCalculator> > AvailableCalculators::features_ = define_features();


const std::map<std::string, boost::shared_ptr<FeatureCalculator> >& AvailableCalculators::get() {
  return features_;
}


} /* namespace pgmlink */
