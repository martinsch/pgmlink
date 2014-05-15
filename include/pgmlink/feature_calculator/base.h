#ifndef FEATURE_CALCULATOR_SQUARE_ROOT_SQUARE_DIFFERENCE_H
#define FEATURE_CALCULATOR_SQUARE_ROOT_SQUARE_DIFFERENCE_H

// stl
#include <string>

// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator/base.h"
#include "pgmlink/feature_extraction.h"

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

namespace helpers {

////
//// class CalculatorLookup
////
class CalculatorLookup {
 public:
  static boost::shared_ptr<FeatureCalculator> extract_calculator(const std::string& name);
  
 private:
  static const std::map<std::string, boost::shared_ptr<FeatureCalculator> > calculator_map_;
};

} // namespace helpers


} // namespace feature_extraction

} // namespace pgmlink


#endif /* FEATURE_CALCULATOR_SQUARE_ROOT_SQUARE_DIFFERENCE_H */

