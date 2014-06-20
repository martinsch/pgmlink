// stl
#include <vector>
#include <string>
#include <numeric>
#include <cassert>
#include <stdexcept>
#include <cmath>
#include <fstream>

// boost
#include <boost/shared_ptr.hpp>

// vigra
#include <vigra/random_forest.hxx>
#include <vigra/random_forest_hdf5_impex.hxx>
#include <vigra/error.hxx>

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


feature_array FeatureCalculator::calculate(const feature_array& ) const {
  throw std::runtime_error("FeatureCalculator \"" + name() + "\" does not take one feature array");
  return feature_array();
}


feature_array FeatureCalculator::calculate(const feature_array&, const feature_array& ) const {
  throw std::runtime_error("FeatureCalculator \"" + name() + "\" does not take two feature arrays");
  return feature_array();
}


feature_array FeatureCalculator::calculate(const feature_array&, const feature_array&, const feature_array& ) const {
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
//// class IdentityCalculator
////
const std::string IdentityCalculator::name_ = "Identity";

const unsigned IdentityCalculator::length = 0;


IdentityCalculator::~IdentityCalculator() {

}


feature_array IdentityCalculator::calculate(const feature_array& f1) const {
  return f1;
}


const std::string& IdentityCalculator::name() const {
  return IdentityCalculator::name_;
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
//// SquareRootSquaredDifferenceCalculator
////
const std::string SquareRootSquaredDifferenceCalculator::name_ = "SqrtSquaredDiff";

const unsigned SquareRootSquaredDifferenceCalculator::length = 1;


SquareRootSquaredDifferenceCalculator::~SquareRootSquaredDifferenceCalculator() {

}


feature_array SquareRootSquaredDifferenceCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  assert(f1.size() == f2.size());
  feature_array ret(length, 0.);
  feature_array::const_iterator f1_it = f1.begin();
  feature_array::const_iterator f2_it = f2.begin();
  for (; f1_it != f1.end(); ++f1_it, ++f2_it) {
    ret[0] += (*f1_it - *f2_it)*(*f1_it - *f2_it);
  }
  ret[0] = sqrt(ret[0]);
  return ret;
}


const std::string& SquareRootSquaredDifferenceCalculator::name() const {
  return SquareRootSquaredDifferenceCalculator::name_;
}


////
//// class VectorDifferenceCalculator
////
const std::string VectorDifferenceCalculator::name_ = "VectorDifference";

const unsigned VectorDifferenceCalculator::length = 0;

VectorDifferenceCalculator::~VectorDifferenceCalculator() {
}

feature_array VectorDifferenceCalculator::calculate(
  const feature_array& f1,
  const feature_array& f2
) const {
  assert(f1.size() == f2.size());
  feature_array ret(f1.size(), 0.0);

  feature_array::const_iterator f1_it = f1.begin();
  feature_array::const_iterator f2_it = f2.begin();
  feature_array::iterator ret_it = ret.begin();

  for (; f1_it != f1.end(); f1_it++, f2_it++, ret_it++) {
    (*ret_it) = *f2_it - *f1_it;
  }

  return ret;
}

const std::string& VectorDifferenceCalculator::name() const {
  return VectorDifferenceCalculator::name_;
}


////
//// Class CurvatureCalculator
////
const std::string CurvatureCalculator::name_ = "CurvatureCalculator";

const unsigned CurvatureCalculator::length = 0;

CurvatureCalculator::~CurvatureCalculator() {
}

feature_array CurvatureCalculator::calculate(
  const feature_array& f1,
  const feature_array& f2,
  const feature_array& f3
) const {
  assert(f1.size() == f2.size());
  assert(f1.size() == f3.size());
  feature_array ret(f1.size(), 0.0);

  feature_array::const_iterator f1_it = f1.begin();
  feature_array::const_iterator f2_it = f2.begin();
  feature_array::const_iterator f3_it = f3.begin();
  feature_array::iterator ret_it = ret.begin();

  for (; f1_it != f1.end(); f1_it++, f2_it++, f3_it++, ret_it++) {
    (*ret_it) = (*f1_it) + (*f3_it) - 2*(*f2_it);
  }

  return ret;
}

const std::string& CurvatureCalculator::name() const {
  return CurvatureCalculator::name_;
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
  feature_array ret(f1.size(), 0.);
  // keep ratio <= 1
  // no zero check, as we do not have empty regions
  for (size_t i = 0; i < ret.size(); ++i) {
    if (f1[i] == f2[i]) {
      ret[i] = 1.0;
    } else if (f1[i] < 0.0001 && f2[i]  < 0.0001) {
      ret[i] = 1.0;
    } else if (f1[i] < f2[i]) {
      ret[i] = f1[i]/f2[i];
    } else {
      ret[i] = f2[i]/f1[i];
    }
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
//// class AsymmetricRatioCalculator
////
const std::string AsymmetricRatioCalculator::name_ = "AsymmetricRatio";

const unsigned AsymmetricRatioCalculator::length = 1;


AsymmetricRatioCalculator::~AsymmetricRatioCalculator() {

}


feature_array AsymmetricRatioCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  assert(f1.size() == f2.size());
  feature_array ret(f1.size(), 0.);
  // keep ratio <= 1
  // no zero check, as we do not have empty regions
  for (size_t i = 0; i < ret.size(); ++i) {
    if (f1[i] == f2[i]) {
      ret[i] = 1.0;
    } else if (f1[i] < 0.0001 && f2[i] < 0.0001) {
      ret[i] = 1.0;
    } else {
      ret[i] = f1[i]/f2[i];
    }
  }
  return ret;
}


feature_array AsymmetricRatioCalculator::calculate(const feature_array&, const feature_array& f2, const feature_array& f3) const {
  return calculate(f2, f3);
}


const std::string& AsymmetricRatioCalculator::name() const {
  return AsymmetricRatioCalculator::name_;
}


////
//// class ParentRatios
////



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
/* const std::string ChildrenMeanParentIntensityRatioCalculator::name_ = "children_parent_intensity_ratio";

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
   } */



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
  assert(f1 != t1.features.end());
  feature_array features1 = f1->second;
  LOG(logDEBUG4) << "FeatureExtractor::extract: feature[0] =  " << features1[0];
  return calculator_->calculate(features1);
}


feature_array FeatureExtractor::extract(const Traxel& t1, const Traxel& t2) const {
  LOG(logDEBUG4) << "FeatureExtractor::extract: feature " << feature_name_;
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


////
//// class FeatureExtractorSelective
////
FeatureExtractorSelective::FeatureExtractorSelective(boost::shared_ptr<FeatureCalculator> calculator, const std::string& feature_name) :
    FeatureExtractor(calculator, feature_name) {

}


FeatureExtractorSelective::~FeatureExtractorSelective() {
  
}


feature_array FeatureExtractorSelective::extract(const Traxel& t1, const Traxel& t2) const {
  LOG(logDEBUG4) << "FeatureExtractorSelective::extract: feature " << feature_name_;
  FeatureMap::const_iterator indices = t1.features.find("intersect_indices");
  FeatureMap::const_iterator intersects = t1.features.find("intersects");
  FeatureMap::const_iterator count1 = t1.features.find("Count");
  FeatureMap::const_iterator count2 = t2.features.find("Count");

  assert(indices != t1.features.end());
  assert(intersects != t1.features.end());
  assert(count1 != t1.features.end());
  assert(count2 != t2.features.end());

  feature_array::const_iterator index = std::find(indices->second.begin(), indices->second.end(), t2.Id);
  size_t offset = index - indices->second.begin();
  float res = 0.0;
  if (index != indices->second.end()) {
    res = intersects->second[offset]/(count1->second[0] + count2->second[0] - intersects->second[offset]);
  }
  return feature_array(1, res);
}


std::string FeatureExtractorSelective::name() const {
  return feature_name_;
}


////
//// class FeatureExtractorDifferentFeatures
////
FeatureExtractorDifferentFeatures::FeatureExtractorDifferentFeatures(boost::shared_ptr<FeatureCalculator> calculator,
                                                                     const std::string& feature_name_1,
                                                                     const std::string& feature_name_2) :
    FeatureExtractor(calculator, feature_name_1),
    feature_name_1_(feature_name_1),
    feature_name_2_(feature_name_2) {

}


FeatureExtractorDifferentFeatures::~FeatureExtractorDifferentFeatures() {
  
}


feature_array FeatureExtractorDifferentFeatures::extract(const Traxel& t1) const {
  LOG(logDEBUG4) << "FeatureExtractorDifferentFeatures::extract: " << calculator_->name() <<  ": features " << name();
  FeatureMap::const_iterator feature_1 = t1.features.find(feature_name_1_);
  FeatureMap::const_iterator feature_2 = t1.features.find(feature_name_2_);

  assert(feature_1 != t1.features.end());
  assert(feature_2 != t1.features.end());

  if (feature_1->second.size() != feature_2->second.size()) {
    throw std::runtime_error("FeatureExtractorDifferentFeatures::extract: features " + name() +
                             " have different sizes");
  }

  return calculator_->calculate(feature_1->second, feature_2->second);
}


std::string FeatureExtractorDifferentFeatures::name() const {
  return feature_name_1_ + "," + feature_name_2_;
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
  feature_map.insert(std::make_pair("ChildrenRatio", calc));

  calc = boost::shared_ptr<FeatureCalculator>(new AsymmetricRatioCalculator);
  feature_map.insert(std::make_pair("AsymmetricRatio", calc));

  calc = boost::shared_ptr<FeatureCalculator>(new AbsoluteDifferenceCalculator);
  feature_map.insert(std::make_pair("AbsDiff", calc));

  calc = boost::shared_ptr<FeatureCalculator>(new SquareRootSquaredDifferenceCalculator);
  feature_map.insert(std::make_pair("SqrtSquaredDiff", calc));

  calc = boost::shared_ptr<FeatureCalculator>(new MaxParentRatio);
  feature_map.insert(std::make_pair("MaxParentRatio", calc));

  calc = boost::shared_ptr<FeatureCalculator>(new MinParentRatio);
  feature_map.insert(std::make_pair("MinParentRatio", calc));

  calc = boost::shared_ptr<FeatureCalculator>(new MeanParentRatio);
  feature_map.insert(std::make_pair("MeanParentRatio", calc));

  calc = boost::shared_ptr<FeatureCalculator>(new MaxParentSquaredDifference);
  feature_map.insert(std::make_pair("MaxParentSquaredDifference", calc));
  
  calc = boost::shared_ptr<FeatureCalculator>(new MinParentSquaredDifference);
  feature_map.insert(std::make_pair("MinParentSquaredDifference", calc));
  
  calc = boost::shared_ptr<FeatureCalculator>(new MeanParentSquaredDifference);
  feature_map.insert(std::make_pair("MeanParentSquaredDifference", calc));
  
  calc = boost::shared_ptr<FeatureCalculator>(new RatioParentSquaredDifference);
  feature_map.insert(std::make_pair("RatioParentSquaredDifference", calc));

  calc = boost::shared_ptr<FeatureCalculator>(new IdentityCalculator);
  feature_map.insert(std::make_pair("Identity", calc));

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


////
//// ClassifierStrategy
////
ClassifierStrategy::ClassifierStrategy(const std::string& name) :
    name_(name) {}


ClassifierStrategy::~ClassifierStrategy() {}


////
//// ClassifierLazy
////
ClassifierLazy::ClassifierLazy() :
    ClassifierStrategy("") {}


ClassifierLazy::~ClassifierLazy() {}


void ClassifierLazy::classify(Traxel&, bool) {

}


void ClassifierLazy::classify(Traxel& ,
                              const Traxel&, bool ) {

}


void ClassifierLazy::classify(const Traxel& ,
                              const Traxel& ,
                              std::map<unsigned, feature_array >& ,
                              bool ) {

}


void ClassifierLazy::classify(const Traxel&,
                              const Traxel&,
                              const Traxel&,
                              std::map<std::pair<unsigned, unsigned>, feature_array >&,
                              bool) {

}

namespace {

void append_to_file(const std::string& filename, const Traxel& trx1, double number) {
  std::ofstream log_file;
  log_file.open(filename.c_str(), std::ios::out | std::ios::app);
  log_file << trx1 << " " << number << '\n';
  log_file.close();
}

void append_to_file(const std::string& filename, const Traxel& trx1, const Traxel& trx2, double number) {
  std::ofstream log_file;
  log_file.open(filename.c_str(), std::ios::out | std::ios::app);
  log_file << trx1 << ',' << trx2 << " " << number << '\n';
  log_file.close();
}

void append_to_file(const std::string& filename, const Traxel& trx1, const Traxel& trx2, const Traxel& trx3, double number) {
  std::ofstream log_file;
  log_file.open(filename.c_str(), std::ios::out | std::ios::app);
  log_file << trx1 << ',' << trx2 << ',' << trx3 << " " << number << '\n';
  log_file.close();
}

}


////
//// ClassifierConstant
////
ClassifierConstant::ClassifierConstant(double probability, const std::string& name) :
    ClassifierStrategy(name),
    probability_(probability) {
  if (probability_ > 1.) {
    throw std::runtime_error("Probability > 1!");
  }
  if (probability_ < 0.) {
    throw std::runtime_error("Probability < 0!");
  }
}


ClassifierConstant::~ClassifierConstant() {}


void ClassifierConstant::classify(Traxel& trax, bool /*with_predict*/) {
  LOG(logDEBUG4) << "ClassifierConstant::classify() -- adding feature " << name_
                 << " to " << trax;
  trax.features[name_] = feature_array(1, 1-probability_);
  trax.features[name_].push_back(probability_);
}


void ClassifierConstant::classify(Traxel& trax_out,
                                  const Traxel&,
                                  bool /*with_predict*/) {
  LOG(logDEBUG3) << "ClassifierConstant::classify() -- entered";
  trax_out.features[name_].push_back(1-probability_);
  trax_out.features[name_].push_back(probability_);
}
                                  


void ClassifierConstant::classify(const Traxel& trax_out,
                                  const Traxel& trax_in,
                                  std::map<unsigned, feature_array >& feature_map,
                                  bool with_predict) {
  LOG(logDEBUG3) << "ClassifierConstant::classify() -- entered";
  LOG(logDEBUG4) << "ClassifierConstant::classify() -- " << trax_out << " -> " << trax_in;
  if (with_predict) {
    assert(feature_map[trax_in.Id].size() == 2);
    feature_map[trax_in.Id][0] = 1-probability_;
    feature_map[trax_in.Id][1] = probability_;
  } else {
    assert(feature_map[trax_in.Id].size() == 0);
    feature_map[trax_in.Id].resize(2);
  }
}

void ClassifierConstant::classify(const Traxel& trax_out,
                                  const Traxel& trax_in_first,
                                  const Traxel& trax_in_second,
                                  std::map<std::pair<unsigned, unsigned>, feature_array >& feature_map,
                                  bool with_predict) {
  LOG(logDEBUG3) << "ClassifierConstant::classify() -- entered";
  LOG(logDEBUG4) << "ClassifierConstant::classify() -- " << trax_out << " -> "
                 << trax_in_first << ',' << trax_in_second;
  feature_array& features = feature_map[std::make_pair(trax_in_first.Id, trax_in_second.Id)];
  if (with_predict) {
    assert(features.size() == 2);
    features[0] = 1-probability_;
    features[1] = probability_;
  } else {
    assert(features.size() == 0);
    features.resize(2);
  }
}


////
//// ClassifierRF
////
ClassifierRF::ClassifierRF(vigra::RandomForest<> rf,
                           const std::vector<boost::shared_ptr<FeatureExtractor> >& feature_extractors,
                           const std::string& name) :
    ClassifierStrategy(name),
    rf_(rf),
    feature_extractors_(feature_extractors),
    features_(vigra::MultiArray<2, feature_type>::difference_type(1, rf.feature_count())),
    probabilities_(vigra::MultiArray<2, feature_type>::difference_type(1, rf.class_count())) {
  assert(feature_extractors.size() == feature_extractors_.size());
  assert(feature_extractors_.size() > 0);
  if(rf_.class_count() < 2) {
    throw std::runtime_error("ClassifierRF -- Random Forest has less than two classes!");
  }
}


ClassifierRF::~ClassifierRF() {}


// ClassifierRF is not doing anything. Create derived class if you want to take action.
void ClassifierRF::classify(Traxel&, bool) {}


void ClassifierRF::classify(Traxel&,
                            const Traxel&, bool) {}


void ClassifierRF::classify(const Traxel&,
                            const Traxel&,
                            std::map<unsigned, feature_array >&,
                            bool) {}


void ClassifierRF::classify(const Traxel&,
                            const Traxel&,
                            const Traxel&,
                            std::map<std::pair<unsigned, unsigned>, feature_array >&,
                            bool) {}


void ClassifierRF::extract_features(const Traxel& t) {
  extract_features(t, features_, probabilities_);
}


void ClassifierRF::extract_features(const Traxel& t,
                                    vigra::MultiArrayView<2, feature_type> features,
                                    vigra::MultiArrayView<2, feature_type> /*probabilities*/) {
  size_t starting_index = 0;
  for (std::vector<boost::shared_ptr<FeatureExtractor> >::const_iterator it = feature_extractors_.begin();
       it != feature_extractors_.end();
       ++it) {
    LOG(logDEBUG4) << "ClassifierRF: extracting " << (*it)->name();
    feature_array feats = (*it)->extract(t);
    LOG(logDEBUG4) << "ClassifierRF: current feats.size(): " << feats.size()
                   << ", features.shape(): " << features.shape();
    assert(starting_index + feats.size() <= features.shape()[1]);
    std::copy(feats.begin(), feats.end(), features.begin() + starting_index);
    starting_index += feats.size();
  }
  if (starting_index != features.shape()[1]) {
    throw std::runtime_error("ClassifierRF -- extracted features size does not match random forest feature size");
  }
}


void ClassifierRF::extract_features(const Traxel& t1, const Traxel& t2) {
  extract_features(t1, t2, features_, probabilities_);
}


void ClassifierRF::extract_features(const Traxel& t1,
                                    const Traxel& t2,
                                    vigra::MultiArrayView<2, feature_type> features,
                                    vigra::MultiArrayView<2, feature_type> /*probabilities*/) {
  size_t starting_index = 0;
  for (std::vector<boost::shared_ptr<FeatureExtractor> >::const_iterator it = feature_extractors_.begin();
       it != feature_extractors_.end();
       ++it) {
    LOG(logDEBUG4) << "ClassifierRF: extracting " << (*it)->name();
    feature_array feats = (*it)->extract(t1, t2);
    assert(starting_index + feats.size() <= features.shape()[1]);
    std::copy(feats.begin(), feats.end(), features.begin() + starting_index);
    starting_index += feats.size();
  }
  if (starting_index != features.shape()[1]) {
    throw std::runtime_error("ClassifierRF -- extracted features size does not match random forest feature size");
  }
}


void ClassifierRF::extract_features(const Traxel& parent, const Traxel& child1, const Traxel& child2) {
  extract_features(parent, child1, child2, features_, probabilities_);
}


void ClassifierRF::extract_features(const Traxel& parent,
                                    const Traxel& child1,
                                    const Traxel& child2,
                                    vigra::MultiArrayView<2, feature_type> features,
                                    vigra::MultiArrayView<2, feature_type> /*probabilities*/) {
  size_t starting_index = 0;
  for (std::vector<boost::shared_ptr<FeatureExtractor> >::const_iterator it = feature_extractors_.begin();
       it != feature_extractors_.end();
       ++it) {
    LOG(logDEBUG4) << "ClassifierRF:extract_features() -- extracting " << (*it)->name();
    feature_array feats = (*it)->extract(parent, child1, child2);
    assert(starting_index + feats.size() <= features.shape()[1]);
    std::copy(feats.begin(), feats.end(), features.begin() + starting_index);
    starting_index += feats.size();
  }
  if (starting_index != features.shape()[1]) {
    throw std::runtime_error("ClassifierRF::extract_features() -- extracted features size does not match random forest feature size");
  }
}


ClassifierMoveRF::ClassifierMoveRF(vigra::RandomForest<> rf, 
                                   const std::vector<boost::shared_ptr<FeatureExtractor> >& feature_extractors,
                                   const std::string& name) :
    ClassifierRF(rf, feature_extractors, name) {}


ClassifierMoveRF::~ClassifierMoveRF() {}


void ClassifierMoveRF::classify(Traxel& trax_out,
                                const Traxel& trax_in, bool /*with_predict*/) {
  extract_features(trax_out, trax_in);
  feature_array& features = trax_out.features[name_];
  if (features.size() == 0) {
    features.push_back(1.0);
    features.push_back(0.0);
  }
  rf_.predictProbabilities(features_, probabilities_);
  if (probabilities_(0, 1) > features[1]) {
    features[0] = probabilities_(0, 0);
    features[1] = probabilities_(0, 1);
  }
  assert(features.size() == rf_.class_count());
}


void ClassifierMoveRF::classify(const Traxel& trax_out,
                                const Traxel& trax_in,
                                std::map<unsigned, feature_array >& feature_map,
                                bool with_predict) {
  feature_array& probabilities_fa = feature_map[trax_in.Id];
  if (with_predict) {
    vigra::MultiArray<2, feature_type> features(features_.shape());
    vigra::MultiArray<2, feature_type> probabilities(probabilities_.shape());
    extract_features(trax_out, trax_in, features, probabilities);
    rf_.predictProbabilities(features, probabilities);
    LOG(logDEBUG4) << "ClassifierMoveRF::classify() -- features[0] = " << features[0];
    if (probabilities_fa.size() != probabilities.shape()[1]) {
      LOG(logDEBUG4) << "ClassifierMoveRF::classify() -- the feature map has not been initialized yet, doing that now.";
      probabilities_fa.resize(probabilities.shape()[1]);
    }
    std::copy(probabilities.begin(),
              probabilities.end(),
              probabilities_fa.begin());
    LOG(logDEBUG4) << "ClassifierMoveRF::classify() -- " << trax_out
                   << " to " << trax_in << "probability: "
                   << probabilities[0] << ',' << probabilities[1];
    LOG(logDEBUG4) << "move_prior,t=" << trax_out.Timestep << ",id=" << trax_out.Id << ",prob=" << probabilities[1] << ",id=" << trax_in.Id;
  } else {
    assert(probabilities_fa.size() == 0);
    probabilities_fa.resize(probabilities_.shape()[1]);
  }
}


ClassifierDivisionRF::ClassifierDivisionRF(vigra::RandomForest<> rf, 
                                           const std::vector<boost::shared_ptr<FeatureExtractor> >& feature_extractors,
                                           const std::string& name) :
    ClassifierRF(rf, feature_extractors, name) {}


ClassifierDivisionRF::~ClassifierDivisionRF() {}


void ClassifierDivisionRF::classify(Traxel& trax_out,
                                    const Traxel& trax_in_first,
                                    const Traxel& trax_in_second,
                                    bool /*with_predict*/) {
  feature_array& probabilities_fa = trax_out.features[name_];
  if (probabilities_fa.size() != 2) {
    probabilities_fa.resize(2);
    probabilities_fa[0] = 1.0;
    probabilities_fa[1] = 0.0;
  }
  extract_features(trax_out, trax_in_first, trax_in_second);
  rf_.predictProbabilities(features_, probabilities_);
  if (probabilities_(0, 1) > probabilities_fa[1]) {
    probabilities_fa[0] = probabilities_(0, 0);
    probabilities_fa[1] = probabilities_(0, 1);
  }
  assert(probabilities_fa.size() == rf_.class_count());
}


void ClassifierDivisionRF::classify(const Traxel& trax_out,
                                    const Traxel& trax_in_first,
                                    const Traxel& trax_in_second,
                                    std::map<std::pair<unsigned, unsigned>, feature_array >& feature_map,
                                    bool with_predict) {
  feature_array& probabilities_fa = feature_map[std::make_pair(std::min(trax_in_first.Id, trax_in_second.Id),
                                                               std::max(trax_in_first.Id, trax_in_second.Id))];
  if (with_predict) {
    vigra::MultiArray<2, feature_type> features(features_.shape());
    vigra::MultiArray<2, feature_type> probabilities(probabilities_.shape());
    extract_features(trax_out, trax_in_first, trax_in_second, features, probabilities);
    rf_.predictProbabilities(features, probabilities);
    if (probabilities_fa.size() != probabilities.shape()[1]) {
      LOG(logDEBUG4) << "ClassifierDivisionRF::classify() -- the feature map has not been initialized yet, doing that now.";
      probabilities_fa.resize(probabilities.shape()[1]);
    }
    std::copy(probabilities.begin(),
              probabilities.end(),
              probabilities_fa.begin());
    LOG(logDEBUG4) << "division_prior,t=" << trax_out.Timestep << ",id=" << trax_out.Id << ",prob=" << probabilities[1] << ",id=" << trax_in_first.Id
                   <<",id=" << trax_in_second.Id;
  } else {
    assert(probabilities_fa.size() == 0);
    probabilities_fa.resize(probabilities_.shape()[1]);
  }
}


ClassifierCountRF::ClassifierCountRF(vigra::RandomForest<> rf,
                                     const std::vector<boost::shared_ptr<FeatureExtractor> >& feature_extractors,
                                     const std::string& name) :
    ClassifierRF(rf, feature_extractors, name) {}


ClassifierCountRF::~ClassifierCountRF() {}


void ClassifierCountRF::classify(Traxel& trax, bool /*with_predict*/) {
  extract_features(trax);
  rf_.predictProbabilities(features_, probabilities_);
  LOG(logDEBUG4) << "ClassifierCountRF::classify() -- adding feature "
                 << name_ << " to " << trax;
  assert(trax.features[name_].size() == 0);
  trax.features[name_].insert(trax.features[name_].end(),
                              probabilities_.begin(),
                              probabilities_.end()
                              );
}


ClassifierDetectionRF::ClassifierDetectionRF(vigra::RandomForest<> rf,
                                             const std::vector<boost::shared_ptr<FeatureExtractor> >& feature_extractors,
                                             const std::string& name) :
    ClassifierRF(rf, feature_extractors, name) {}


ClassifierDetectionRF::~ClassifierDetectionRF() {}


void ClassifierDetectionRF::classify(Traxel& trax, bool with_predict) {
  feature_array& probabilities_fa = trax.features[name_];
  if (with_predict) {
    vigra::MultiArray<2, feature_type> features(features_.shape());
    vigra::MultiArray<2, feature_type> probabilities(probabilities_.shape());
    extract_features(trax, features, probabilities);
    rf_.predictProbabilities(features, probabilities);
    LOG(logDEBUG4) << "ClassifierDetectionRF::classify() -- features[0] = "
                   << features[0];
    LOG(logDEBUG4) << "detection_prior,t=" << trax.Timestep << ",id=" << trax.Id << ",prob=" << probabilities[1];

    if (probabilities_fa.size() == 0) {
      LOG(logDEBUG4) << "ClassifierDetectionRF::classify() -- the feature map has not been initialized yet, doing that now.";
      probabilities_fa.insert(probabilities_fa.end(),
                              probabilities.begin(),
                              probabilities.end()
                              );
        assert(probabilities_fa.size() == 2);
    } else {
      std::copy(probabilities.begin(),
                probabilities.end(),
                probabilities_fa.begin());
      assert(probabilities_fa.size() == 2);
      LOG(logDEBUG4) << "ClassifierDetectionRF::classify() -- adding feature "
                     << name_ << " to " << trax << ':'
                     << *(probabilities.begin()) << ',' << *(probabilities.begin()+1) << '\n'
                     << probabilities_fa[0] << ',' << probabilities_fa[1];
    }
  } else {
    assert(probabilities_fa.size() == 0);
    probabilities_fa.resize(probabilities_.shape()[1], 0.);
    assert(probabilities_fa.size() == 2);
  }
}


////
//// class ClassifierStrategyBuilder
////
boost::shared_ptr<ClassifierStrategy> ClassifierStrategyBuilder::build(const Options& options) {
  boost::shared_ptr<ClassifierStrategy> classifier;
  
  if (options.type == CONSTANT) {
    classifier = boost::shared_ptr<ClassifierStrategy>(new ClassifierConstant(options.constant_probability,
                                                                              options.name)
                                                       );

  } else {
    vigra::RandomForest<> rf;
    try {
      LOG(logINFO) << "ClassifierStrategyBuilder::build() -- loading random forest from "
                   << options.rf_filename << "/" << options.rf_internal_path;
      vigra::rf_import_HDF5(rf, options.rf_filename, options.rf_internal_path);
    }
    catch (vigra::PreconditionViolation e) {
      if (options.constant_classifier_fallback) {
        LOG(logINFO) << "ClassifierStrategyBuilder::build() -- Could not load random forest from "
                     << options.rf_filename << "/" << options.rf_internal_path
                     << " -> falling back to constant classifier";
        Options options_fallback = options;
        options_fallback.type = CONSTANT;
        return build(options_fallback);
      } else {
        throw std::runtime_error("Could not load random forest from "
                                 + options.rf_filename + "/" + options.rf_internal_path);
      }
    }
    
    std::vector<boost::shared_ptr<FeatureExtractor> > extractors;
    const std::map<std::string, boost::shared_ptr<FeatureCalculator> >& cmap = AvailableCalculators::get();
    for (std::vector<std::pair<std::string, std::string> >::const_iterator feature = options.feature_list.begin();
         feature != options.feature_list.end();
         ++feature) {
      LOG(logDEBUG4) << "Creating Extractor for " << feature->first << ": " << feature->second;
      if (feature->first == "IntersectionUnionRatio") {
        extractors.push_back(boost::shared_ptr<FeatureExtractorSelective>(new FeatureExtractorSelective(cmap.find("Identity")->second, "IntersectionUnionRatio")));
      } else if (feature->first[0] == '_') {
        std::map<std::string, boost::shared_ptr<FeatureCalculator> >::const_iterator c = cmap.find(feature->first.substr(1));
        if (c == cmap.end()) {
          throw std::runtime_error("Calculator \"" + feature->first.substr(1) + "\" not available!");
        } else {
          size_t pos = feature->second.find(",");
          if (pos == feature->second.npos) {
            throw std::runtime_error("Features (" + feature->second + ") are not split by ,");
          }
          extractors.push_back(boost::shared_ptr<FeatureExtractor>(new FeatureExtractorDifferentFeatures(c->second,
                                                                                                         feature->second.substr(0, pos),
                                                                                                         feature->second.substr(pos+1)
                                                                                                         )
                                                                   )
                               );
        }        
      } else {
        std::map<std::string, boost::shared_ptr<FeatureCalculator> >::const_iterator c = cmap.find(feature->first);
        if (c == cmap.end()) {
          throw std::runtime_error("Calculator \"" + feature->first + "\" not available!");
        } else {
          extractors.push_back(boost::shared_ptr<FeatureExtractor>(new FeatureExtractor(c->second, feature->second)));
        }
      }
      LOG(logDEBUG4) << "ClassifierStrategyBuilder::builder() -- Created Feature Extractor " << (*extractors.rbegin())->name();
    }

    
    if (options.type == RF_MOVE) {
      LOG(logDEBUG) << "ClassifierStrategyBuilder::build: create MoveRF";
      classifier = boost::shared_ptr<ClassifierStrategy>(new ClassifierMoveRF(rf,
                                                                              extractors,
                                                                              options.name)
                                                         );
    }
    else if (options.type == RF_DIVISION) {
      classifier = boost::shared_ptr<ClassifierStrategy>(new ClassifierDivisionRF(rf,
                                                                                  extractors,
                                                                                  options.name)
                                                         );
    } else if (options.type == RF_COUNT) {
      classifier = boost::shared_ptr<ClassifierStrategy>(new ClassifierCountRF(rf,
                                                                               extractors,
                                                                               options.name)
                                                         );
    } else if (options.type == RF_DETECTION) {
      classifier = boost::shared_ptr<ClassifierStrategy>(new ClassifierDetectionRF(rf,
                                                                                   extractors,
                                                                                   options.name)
                                                         );
    }  else {
      throw std::runtime_error("Not a valid choice for classifier!");
    }
  }

  return classifier;
}


} /* namespace pgmlink */
