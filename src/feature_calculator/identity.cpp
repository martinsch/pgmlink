#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator/identity.h"

namespace pgmlink {

namespace feature_extraction {

////
//// class IdentityCalculator
////
const std::string IdentityCalculator::name_ = "Identity";

IdentityCalculator::~IdentityCalculator() {

}


feature_array IdentityCalculator::calculate(const feature_array& f1) const {
  return f1;
}


const std::string& IdentityCalculator::name() const {
  return IdentityCalculator::name_;
}

} // namespace feature_extraction

} // namespace pgmlink
