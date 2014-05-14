#ifndef FEATURE_EXTRACTION_STRATEGY_H
#define FEATURE_EXTRACTION_STRATEGY_H

// stl
#include <string>
#include <map>

// boost
#include <boost/shared_ptr.hpp>

// pgmlink
#include "feature.h"

namespace pgmlink {

namespace feature_extraction {


////
//// class FeatureCalculator
////
class FeatureCalculator {
 public:
  virtual ~FeatureCalculator();
  virtual feature_array calculate(const feature_array& /* f1 */ ) const;
  virtual feature_array calculate(const feature_array& /* f1 */, const feature_array& /* f2 */) const;
  virtual feature_array calculate(const feature_array& /* f1 */, const feature_array& /* f2 */, const feature_array& /* f3 */) const;
  virtual const std::string& name() const;

  bool operator==(const FeatureCalculator& other);
  bool operator!=(const FeatureCalculator& other);

 private:
  static const std::string name_;
};


////
//// class ElementWiseSquaredDistanceCalculator
////
class ElementWiseSquaredDistanceCalculator : public FeatureCalculator {
 public:
  virtual ~ElementWiseSquaredDistanceCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;

  virtual const std::string& name() const;

 private:
  static const std::string name_;
};


////
//// class FeatureExtractor
////
class FeatureExtractor {
 public:
  FeatureExtractor(boost::shared_ptr<FeatureCalculator> calculator, const std::string& feature_name);
  virtual ~FeatureExtractor();
  virtual feature_array extract(const Traxel& t1) const;
  virtual feature_array extract(const Traxel& t1, const Traxel& t2) const;
  virtual feature_array extract(const Traxel& t1, const Traxel& t2, const Traxel& t3) const;
  boost::shared_ptr<FeatureCalculator> calculator() const;
  virtual std::string name() const;

 protected:
  boost::shared_ptr<FeatureCalculator> calculator_;
  std::string feature_name_;
};


////
//// class MultipleFeaturesExtraction
////
class MultipleFeatureExtraction {
  
};



////
//// helpers
////
namespace helpers 
{

class CalculatorLookup {
 public:
  static boost::shared_ptr<FeatureCalculator> extract_calculator(const std::string& name);
  
 private:
  static const std::map<std::string, boost::shared_ptr<FeatureCalculator> > calculator_map_;
};

} /* namespace helpers */

} /* namespace feature_extraction */

} /* namespace pgmlink */

#endif /* FEATURE_EXTRACTION_STRATEGY_H */
