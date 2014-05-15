#ifndef FEATURE_CALCULATOR_BASE_H
#define FEATURE_CALCULATOR_BASE_H

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

} // namespace feature_extraction

} // namespace pgmlink


#endif /* FEATURE_CALCULATOR_BASE_H */

