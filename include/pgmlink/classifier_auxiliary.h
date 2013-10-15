#ifndef CLASSIFIER_AUXILIARY_H
#define CLASSIFIER_AUXILIARY_H

// stl
#include <vector>
#include <string>

// boost
#include <boost/shared_ptr.hpp>

// vigra

// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/traxels.h"

namespace pgmlink {

class FeatureCalculator {
 public:
  static const std::string name;
  static const unsigned length;

  virtual ~FeatureCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const;

  bool operator==(const FeatureCalculator& other);
};


class DistanceCalculator : public FeatureCalculator {
 public:
  static const std::string name;
  static const unsigned length;

  virtual ~DistanceCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
};

  

class FeatureExtractor {
 public:
  FeatureExtractor(boost::shared_ptr<FeatureCalculator> calculator, const std::vector<std::string>& feature_names);
  feature_array extract(const Traxel& t1, const Traxel& t2) const;
  feature_array extract(const Traxel& t1, const Traxel& t2, const Traxel& t3) const;
  boost::shared_ptr<FeatureCalculator> calculator();

 private:
  boost::shared_ptr<FeatureCalculator> calculator_;
  std::vector<std::string> feature_names_;
};
  

} /* namesapce pgmlink */

#endif /* CLASSIFIER_H */
