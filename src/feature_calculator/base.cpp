// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator/base.h"

namespace pgmlink {

namespace feature_extraction {

////
//// class FeatureCalculator
////
FeatureCalculator::~FeatureCalculator() {

}


feature_array FeatureCalculator::calculate(const feature_array& /* f1 */) const {
  throw std::runtime_error("Not implemented yet!");
}


feature_array FeatureCalculator::calculate(const feature_array& /* f1 */, const feature_array& /* f2 */) const {
  throw std::runtime_error("Not implemented yet!");
}


feature_array FeatureCalculator::calculate(const feature_array& /* f1 */, const feature_array& /* f2 */, const feature_array& /* f3 */) const {
  throw std::runtime_error("Not implemented yet!");
}


const std::string& FeatureCalculator::name() const {
  return name_;
}


bool FeatureCalculator::operator==(const FeatureCalculator& other) {
  return this->name_ == other.name();
}


bool FeatureCalculator::operator!=(const FeatureCalculator& other) {
  return !(*this == other);
}


const std::string FeatureCalculator::name_ = "";


} // namespace feature_extraction

} // namespace pgmlink
