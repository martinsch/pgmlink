#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator/squared_difference.h"

namespace pgmlink {

namespace feature_extraction {

////
//// class SquaredDifferenceCalculator
////
const std::string SquaredDifferenceCalculator::name_ = "SquaredDiff";


SquaredDifferenceCalculator::~SquaredDifferenceCalculator() {

}


feature_array SquaredDifferenceCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  assert(f1.size() == f2.size());
  feature_array ret(f1.size(), 0.);
  feature_array::const_iterator f1_it = f1.begin();
  feature_array::const_iterator f2_it = f2.begin();
  for (; f1_it != f1.end(); ++f1_it, ++f2_it) {
    ret[0] += (*f1_it - *f2_it)*(*f1_it - *f2_it);
  }
  return ret;
}


const std::string& SquaredDifferenceCalculator::name() const {
  return SquaredDifferenceCalculator::name_;
}


} // namespace feature_extraction

} // namespace pgmlink
