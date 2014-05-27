#ifndef FEATURE_CALCULATOR_IDENTITY_H
#define FEATURE_CALCULATOR_IDENTITY_H

// stl
#include <string>

// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator/base.h"
#include "pgmlink/feature_extraction.h"

namespace pgmlink {

namespace feature_extraction {
////
//// IdentityCalculator
////
class IdentityCalculator : public FeatureCalculator {
 public:
  virtual ~IdentityCalculator();
  virtual feature_array calculate(const feature_array& f1) const;
  virtual const std::string& name() const;
private:
static const std::string name_;
};

} // namespace feature_extraction

} // namespace pgmlink


#endif /* FEATURE_CALCULATOR_IDENTITY_H */

