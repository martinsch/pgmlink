#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator/absolute_difference.h"

namespace pgmlink {

namespace feature_extraction {

////
//// class AbsoluteDifferenceCalculator
////
const std::string AbsoluteDifferenceCalculator::name_ = "AbsDiff";


AbsoluteDifferenceCalculator::~AbsoluteDifferenceCalculator() {

}


feature_array AbsoluteDifferenceCalculator::calculate(const feature_array&f1, const feature_array& f2) const {
  assert(f1.size() == f2.size());
  feature_array ret(f1.size(), 0.);
  feature_array::const_iterator f1_it = f1.begin();
  feature_array::const_iterator f2_it = f2.begin();
  for (; f1_it != f1.end(); ++f1_it, ++f2_it) {
    float res = *f1_it - *f2_it;
    ret[0] += res > 0 ? res : -res;
  }
  return ret;
}


const std::string& AbsoluteDifferenceCalculator::name() const {
  return AbsoluteDifferenceCalculator::name_;
}

} // namespace feature_extraction

} // namespace pgmlink
