#ifndef FEATURE_CALCULATOR_ELEMENTWISE_SQUARED_DIFFERENCE_H
#define FEATURE_CALCULATOR_ELEMENTWISE_SQUARED_DIFFERENCE_H

// stl
#include <string>

// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator/base.h"
#include "pgmlink/feature_extraction.h"

namespace pgmlink {

namespace feature_extraction {

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




} // namespace feature_extraction

} // namespace pgmlink


#endif /* FEATURE_CALCULATOR_ELEMENTWISE_SQUARED_DIFFERENCE_H */

