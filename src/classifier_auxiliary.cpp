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

// omp
#include <omp.h>

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


feature_array FeatureCalculator::calculate(const feature_array& f1) const {
  throw std::runtime_error("FeatureCalculator \"" + name() + "\" does not take one feature array");
  return feature_array();
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
  if (f1[0] == f2[0]) {
    ret[0] = 1;
  } else if (f1[0] < f2[0]) {
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


void ClassifierLazy::classify(std::vector<Traxel>& traxels, bool with_predict) {

}


void ClassifierLazy::classify(std::vector<Traxel>& traxels_out,
              const std::vector<Traxel>& traxels_in, bool with_predict) {

}


void ClassifierLazy::classify(const std::vector<Traxel>& traxels_out,
              const std::vector<Traxel>& traxels_in,
              std::map<Traxel, std::map<Traxel, feature_array> >& feature_map,
              bool with_predict) {

}


void ClassifierLazy::classify(const std::vector<Traxel>& traxels_out,
              const std::vector<Traxel>& traxels_in,
              std::map<Traxel, std::map<std::pair<Traxel, Traxel>, feature_array> >& feature_map,
              bool with_predict) {

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


void ClassifierConstant::classify(std::vector<Traxel>& traxels, bool with_predict) {
  for (std::vector<Traxel>::iterator t = traxels.begin();
       t != traxels.end();
       ++t) {
    LOG(logDEBUG4) << "ClassifierConstant::classify() -- adding feature " << name_
                   << " to " << *t;
    t->features[name_] = feature_array(1, 1-probability_);
    t->features[name_].push_back(probability_);
  }
}


void ClassifierConstant::classify(std::vector<Traxel>& traxels_out,
                                  const std::vector<Traxel>&, bool with_predict) {
  LOG(logDEBUG3) << "ClassifierConstant::classify() -- entered";
  for (std::vector<Traxel>::iterator out = traxels_out.begin(); out != traxels_out.end(); ++out) {
    out->features[name_].push_back(1-probability_);
    out->features[name_].push_back(probability_);
  }
}
                                  


void ClassifierConstant::classify(const std::vector<Traxel>& traxels_out,
                                  const std::vector<Traxel>& traxels_in,
                                  std::map<Traxel, std::map<Traxel, feature_array> >& feature_map,
                                  bool with_predict) {
  LOG(logDEBUG3) << "ClassifierConstant::classify() -- entered";
  for (std::vector<Traxel>::const_iterator out = traxels_out.begin(); out != traxels_out.end(); ++out) {
    for (std::vector<Traxel>::const_iterator in = traxels_in.begin(); in != traxels_in.end(); ++in) {
      LOG(logDEBUG4) << "ClassifierConstant::classify() -- " << *out << " -> " << *in;
      // assert(feature_map[*out][*in].size() == 0);
      feature_map[*out][*in].push_back(1-probability_);
      feature_map[*out][*in].push_back(probability_);
    }
  }
}

void ClassifierConstant::classify(const std::vector<Traxel>& traxels_out,
                                  const std::vector<Traxel>& traxels_in,
                                  std::map<Traxel, std::map<std::pair<Traxel, Traxel>, feature_array> >& feature_map,
                                  bool with_predict) {
  for (std::vector<Traxel>::const_iterator out = traxels_out.begin(); out != traxels_out.end(); ++out) {
    for (std::vector<Traxel>::const_iterator child1 = traxels_in.begin(); child1 != traxels_in.end(); ++child1) {
      for (std::vector<Traxel>::const_iterator child2 = child1 + 1; child2 != traxels_in.end(); ++child2) {
        assert(feature_map[*out][std::make_pair(*child1, *child2)].size() == 0);
        feature_map[*out][std::make_pair(*child1, *child2)].push_back(1-probability_);
        feature_map[*out][std::make_pair(*child1, *child2)].push_back(probability_);
      }
    }
  }
}


////
//// ClassifierRF
////
ClassifierRF::ClassifierRF(vigra::RandomForest<> rf,
                           const std::vector<FeatureExtractor>& feature_extractors,
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
void ClassifierRF::classify(std::vector<Traxel>& traxels, bool with_predict) {}


void ClassifierRF::classify(std::vector<Traxel>& traxels_out,
                            const std::vector<Traxel>& traxels_in, bool with_predict) {}


void ClassifierRF::classify(const std::vector<Traxel>& traxels_out,
                            const std::vector<Traxel>& traxels_in,
                            std::map<Traxel, std::map<Traxel, feature_array> >& feature_map,
                            bool with_predict) {}


void ClassifierRF::classify(const std::vector<Traxel>& traxels_out,
                            const std::vector<Traxel>& traxels_in,
                            std::map<Traxel, std::map<std::pair<Traxel, Traxel>, feature_array> >& feature_map,
                            bool with_predict) {}


void ClassifierRF::extract_features(const Traxel& t) {
  extract_features(t, features_, probabilities_);
}


void ClassifierRF::extract_features(const Traxel& t,
                                    vigra::MultiArrayView<2, feature_type> features,
                                    vigra::MultiArrayView<2, feature_type> probabilities) {
  size_t starting_index = 0;
  for (std::vector<FeatureExtractor>::const_iterator it = feature_extractors_.begin();
       it != feature_extractors_.end();
       ++it) {
    LOG(logDEBUG4) << "ClassifierRF: extracting " << it->name();
    feature_array feats = it->extract(t);
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
                                    vigra::MultiArrayView<2, feature_type> probabilities) {
  size_t starting_index = 0;
  for (std::vector<FeatureExtractor>::const_iterator it = feature_extractors_.begin();
       it != feature_extractors_.end();
       ++it) {
    LOG(logDEBUG4) << "ClassifierRF: extracting " << it->name();
    feature_array feats = it->extract(t1, t2);
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
                                    vigra::MultiArrayView<2, feature_type> probabilities) {
  size_t starting_index = 0;
  for (std::vector<FeatureExtractor>::const_iterator it = feature_extractors_.begin();
       it != feature_extractors_.end();
       ++it) {
    LOG(logDEBUG4) << "ClassifierRF:extract_features() -- extracting " << it->name();
    feature_array feats = it->extract(parent, child1, child2);
    assert(starting_index + feats.size() <= features.shape()[1]);
    std::copy(feats.begin(), feats.end(), features.begin() + starting_index);
    starting_index += feats.size();
  }
  if (starting_index != features.shape()[1]) {
    throw std::runtime_error("ClassifierRF::extract_features() -- extracted features size does not match random forest feature size");
  }
}


ClassifierMoveRF::ClassifierMoveRF(vigra::RandomForest<> rf, 
                                   const std::vector<FeatureExtractor>& feature_extractors,
                                   const std::string& name) :
    ClassifierRF(rf, feature_extractors, name) {}


ClassifierMoveRF::~ClassifierMoveRF() {}


void ClassifierMoveRF::classify(std::vector<Traxel>& traxels_out,
                                const std::vector<Traxel>& traxels_in, bool with_predict) {
  for (std::vector<Traxel>::iterator out = traxels_out.begin(); out != traxels_out.end(); ++out) {
    vigra::MultiArray<2, feature_type> max_prob(vigra::Shape2(1, rf_.class_count()), 0.);
    for (std::vector<Traxel>::const_iterator in = traxels_in.begin(); in != traxels_in.end(); ++in) {
      extract_features(*out, *in);
      rf_.predictProbabilities(features_, probabilities_);
      if (probabilities_(0, 1) > max_prob(0, 1)) {
        max_prob.swapData(probabilities_);
      }
    }
    assert(out->features[name_].size() == 0);
	out->features[name_].insert(out->features[name_].begin(), max_prob.begin(), max_prob.end());
    assert(out->features[name_].size() == rf_.class_count());
  }
}


void ClassifierMoveRF::classify(const std::vector<Traxel>& traxels_out,
                                const std::vector<Traxel>& traxels_in,
                                std::map<Traxel, std::map<Traxel, feature_array> >& feature_map,
                                bool with_predict) {
  vigra::MultiArray<2, feature_type> features(features_.shape());
  vigra::MultiArray<2, feature_type> probabilities(probabilities_.shape());
  for (std::vector<Traxel>::const_iterator out = traxels_out.begin(); out != traxels_out.end(); ++out) {
    for (std::vector<Traxel>::const_iterator in = traxels_in.begin(); in != traxels_in.end(); ++in) {
      if (with_predict) {
        extract_features(*out, *in, features, probabilities);
        rf_.predictProbabilities(features, probabilities);
        LOG(logDEBUG4) << "ClassifierMoveRF::classify() -- features[0] = " << features[0];

        if (feature_map[*out][*in].size() == 0) {
          LOG(logDEBUG4) << "ClassifierMoveRF::classify() -- the feature map has not been initialized yet, doing that now.";
          feature_map[*out][*in] = feature_array(probabilities.shape()[1]);
        }
        std::copy(probabilities.begin(),
                  probabilities.end(),
                  feature_map[*out][*in].begin());
        LOG(logDEBUG4) << "ClassifierMoveRF::classify() -- " << *out
                       << " to " << *in << "probability: "
                       << probabilities[0] << ',' << probabilities[1];
        LOG(logDEBUG4) << "move_prior,t=" << out->Timestep << ",id=" << out->Id << ",prob=" << probabilities[1] << ",id=" << in->Id;
        append_to_file("classifier_move.log", *out, *in, probabilities[1]);
      } else {
        assert(feature_map[*out][*in].size() == 0);
        feature_map[*out][*in] = feature_array(probabilities.shape()[1], 0.);
      }

    }
  }
}


ClassifierDivisionRF::ClassifierDivisionRF(vigra::RandomForest<> rf, 
                                           const std::vector<FeatureExtractor>& feature_extractors,
                                           const std::string& name) :
    ClassifierRF(rf, feature_extractors, name) {}


ClassifierDivisionRF::~ClassifierDivisionRF() {}


void ClassifierDivisionRF::classify(std::vector<Traxel>& traxels_out,
                                    const std::vector<Traxel>& traxels_in,
                                    bool with_predict) {
  for (std::vector<Traxel>::iterator out = traxels_out.begin(); out != traxels_out.end(); ++out) {
    vigra::MultiArray<2, feature_type> max_prob(vigra::Shape2(1, rf_.class_count()), 0.);
    for (std::vector<Traxel>::const_iterator child1 = traxels_in.begin(); child1 != traxels_in.end(); ++child1) {
      for (std::vector<Traxel>::const_iterator child2 = child1 + 1; child2 != traxels_in.end(); ++child2) {
        extract_features(*out, *child1, *child2);
        rf_.predictProbabilities(features_, probabilities_);
        if (probabilities_(0, 1) > max_prob(0, 1)) {
          max_prob.swapData(probabilities_);
        }
      }
    }
    assert(out->features[name_].size() == 0);
    out->features[name_].insert(out->features[name_].begin(), max_prob.begin(), max_prob.end());
    assert(out->features[name_].size() == rf_.class_count());
  }
}


void ClassifierDivisionRF::classify(const std::vector<Traxel>& traxels_out,
                                    const std::vector<Traxel>& traxels_in,
                                    std::map<Traxel, std::map<std::pair<Traxel, Traxel>, feature_array> >& feature_map,
                                    bool with_predict) {
  vigra::MultiArray<2, feature_type> features(features_.shape());
  vigra::MultiArray<2, feature_type> probabilities(probabilities_.shape());
  for (std::vector<Traxel>::const_iterator out = traxels_out.begin(); out != traxels_out.end(); ++out) {
    for (std::vector<Traxel>::const_iterator child1 = traxels_in.begin(); child1 != traxels_in.end(); ++child1) {
      for (std::vector<Traxel>::const_iterator child2 = child1 + 1; child2 != traxels_in.end(); ++child2) {
        if (with_predict) {
          extract_features(*out, *child1, *child2, features, probabilities);
        	rf_.predictProbabilities(features, probabilities);
        	if (feature_map[*out][std::make_pair(*child1, *child2)].size() == 0) {
                  LOG(logDEBUG4) << "ClassifierDivisionRF::classify() -- the feature map has not been initialized yet, doing that now.";
				  feature_map[*out][std::make_pair(*child1, *child2)] = feature_array(probabilities.shape()[1]);
                }
        	std::copy(probabilities.begin(),
                          probabilities.end(),
                          feature_map[*out][std::make_pair(*child1, *child2)].begin());
                LOG(logDEBUG4) << "division_prior,t=" << out->Timestep << ",id=" << out->Id << ",prob=" << probabilities[1] << ",id=" << child1->Id
                               <<",id=child2->Id";
                append_to_file("classifier_division.log", *out, *child1, *child2, probabilities[1]);
        } else {
        	assert(feature_map[*out][std::make_pair(*child1, *child2)].size() == 0);
        	feature_map[*out][std::make_pair(*child1, *child2)] = feature_array(probabilities.shape()[1], 0.);
        }
      }
    }
  }
}


ClassifierCountRF::ClassifierCountRF(vigra::RandomForest<> rf,
                                     const std::vector<FeatureExtractor>& feature_extractors,
                                     const std::string& name) :
    ClassifierRF(rf, feature_extractors, name) {}


ClassifierCountRF::~ClassifierCountRF() {}


void ClassifierCountRF::classify(std::vector<Traxel>& traxels, bool with_predict) {
  extract_features(traxels[0]);
  rf_.predictProbabilities(features_, probabilities_);
  LOG(logDEBUG4) << "ClassifierCountRF::classify() -- adding feature "
                 << name_ << " to " << traxels[0];
  assert(traxels[0].features[name_].size() == 0);
  traxels[0].features[name_].insert(traxels[0].features[name_].end(),
                                    probabilities_.begin(),
                                    probabilities_.end()
                                    );
}


ClassifierDetectionRF::ClassifierDetectionRF(vigra::RandomForest<> rf,
                                             const std::vector<FeatureExtractor>& feature_extractors,
                                             const std::string& name) :
    ClassifierRF(rf, feature_extractors, name) {}


ClassifierDetectionRF::~ClassifierDetectionRF() {}


void ClassifierDetectionRF::classify(std::vector<Traxel>& traxels, bool with_predict) {
  vigra::MultiArray<2, feature_type> features(features_.shape());
  vigra::MultiArray<2, feature_type> probabilities(probabilities_.shape());
  for (std::vector<Traxel>::iterator t = traxels.begin();
       t != traxels.end();
       ++t) {
    if (with_predict) {
      extract_features(*t, features, probabilities);
      rf_.predictProbabilities(features, probabilities);
      append_to_file("classifier_detection.log", *t, probabilities[1]);
      LOG(logDEBUG4) << "ClassifierDetectionRF::classify() -- features[0] = "
                     << features[0];
      LOG(logDEBUG4) << "detection_prior,t=" << t->Timestep << ",id=" << t->Id << ",prob=" << probabilities[1];

      if (t->features[name_].size() == 0) {
        LOG(logDEBUG4) << "ClassifierDetectionRF::classify() -- the feature map has not been initialized yet, doing that now.";
        t->features[name_].insert(t->features[name_].end(),
                                  probabilities.begin(),
                                  probabilities.end()
                                  );
        assert(t->features[name_].size() == 2);
      } else {

        std::copy(probabilities.begin(),
                  probabilities.end(),
                  t->features[name_].begin());
        assert(t->features[name_].size() == 2);
        LOG(logDEBUG4) << "ClassifierDetectionRF::classify() -- adding feature "
                       << name_ << " to " << *t << ':'
                       << *(probabilities.begin()) << ',' << *(probabilities.begin()+1) << '\n'
                   << t->features[name_][0] << ',' << t->features[name_][1];
      }
    } else {
    assert(t->features[name_].size() == 0);
    t->features[name_].resize(probabilities.shape()[1], 0.);
    assert(t->features[name_].size() == 2);
  }
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
    
    std::vector<FeatureExtractor> extractors;
    const std::map<std::string, boost::shared_ptr<FeatureCalculator> >& cmap = AvailableCalculators::get();
    for (std::vector<std::pair<std::string, std::string> >::const_iterator feature = options.feature_list.begin();
         feature != options.feature_list.end();
         ++feature) {
      std::map<std::string, boost::shared_ptr<FeatureCalculator> >::const_iterator c = cmap.find(feature->first);
      if (c == cmap.end()) {
        throw std::runtime_error("Calculator \"" + feature->first + "\" not available!");
      } else {
        extractors.push_back(FeatureExtractor(c->second, feature->second));
        LOG(logDEBUG4) << "ClassifierStrategyBuilder::builder() -- Created Feature Extractor " << extractors.rbegin()->name();
      }
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
