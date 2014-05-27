#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator/square_root_squared_difference.h"

namespace pgmlink {

namespace feature_extraction {
////
//// SquareRootSquaredDifferenceCalculator
////
const std::string SquareRootSquaredDifferenceCalculator::name_ = "SqrtSquaredDiff";

SquareRootSquaredDifferenceCalculator::~SquareRootSquaredDifferenceCalculator() {

}


feature_array SquareRootSquaredDifferenceCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  assert(f1.size() == f2.size());
  feature_array ret(f1.size(), 0.);
  feature_array::const_iterator f1_it = f1.begin();
  feature_array::const_iterator f2_it = f2.begin();
  for (; f1_it != f1.end(); ++f1_it, ++f2_it) {
    ret[0] += (*f1_it - *f2_it)*(*f1_it - *f2_it);
  }
  ret[0] = sqrt(ret[0]);
  return ret;
}


const std::string& SquareRootSquaredDifferenceCalculator::name() const {
  return SquareRootSquaredDifferenceCalculator::name_;
}


} // namespace feature_extraction

} // namespace pgmlink
