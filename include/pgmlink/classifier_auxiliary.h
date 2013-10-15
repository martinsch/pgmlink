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
  static const std::string name_;
  static const unsigned length;

  virtual ~FeatureCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const;
  virtual const std::string& name() const;

  bool operator==(const FeatureCalculator& other);
};


class DistanceCalculator : public FeatureCalculator {
 public:
  static const std::string name_;
  static const unsigned length;

  virtual ~DistanceCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual const std::string& name() const;
};


class SizeRatioCalculator : public FeatureCalculator {
 public:
  static const std::string name_;
  static const unsigned length;

  virtual ~SizeRatioCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const;
  virtual const std::string& name() const;
};


class IntensityRatioCalculator : public FeatureCalculator {
 public:
  static const std::string name_;
  static const unsigned length;

  virtual ~IntensityRatioCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const;
  virtual const std::string& name() const;
};


class ChildrenMeanParentIntensityRatioCalculator : public FeatureCalculator {
  public:
  static const std::string name_;
  static const unsigned length;

  virtual ~ChildrenMeanParentIntensityRatioCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const;
  virtual const std::string& name() const;
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
