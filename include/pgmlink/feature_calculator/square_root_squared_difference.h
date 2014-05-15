#ifndef FEATURE_CALCULATOR_SQUARE_ROOT_SQUARE_DIFFERENCE_H
#define FEATURE_CALCULATOR_SQUARE_ROOT_SQUARE_DIFFERENCE_H

// stl
#include <string>

// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator/base.h"


namespace pgmlink {

namespace feature_extraction {
////
//// SquareRootSquaredDifferenceCalculator
////
class SquareRootSquaredDifferenceCalculator : public FeatureCalculator {
 public:
  static const std::string name_;

  virtual ~SquareRootSquaredDifferenceCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual const std::string& name() const;
};

} // namespace feature_extraction

} // namespace pgmlink


#endif /* FEATURE_CALCULATOR_SQUARE_ROOT_SQUARE_DIFFERENCE_H */

