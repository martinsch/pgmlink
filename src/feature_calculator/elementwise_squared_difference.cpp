#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator/elementwise_squared_difference.h"

namespace pgmlink {

namespace feature_extraction {

////
//// class ElementWiseSquaredDistanceCalculator
////
ElementWiseSquaredDistanceCalculator::~ElementWiseSquaredDistanceCalculator() {

}


feature_array ElementWiseSquaredDistanceCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  LOG(logDEBUG4) << "ElementWiseSquaredDistanceCalculator::calculate(const feature_array& f1, const feature_array& f2) -- entered";
  assert(f1.size() == f2.size() && "feature array sizes must agree.");
  feature_array result = feature_array(f1.size());
  {
    feature_array::iterator res_it = result.begin();
    for (
             feature_array::const_iterator it1 = f1.begin(),
                 it2 = f2.begin();
             it1 != f1.end();
             ++it1, ++it2, ++res_it) {
      feature_type diff = *it1 - *it2;
      *res_it = diff * diff;
    }
  }
  LOG(logDEBUG4) << "ElementWiseSquaredDistanceCalculator::calculate(const feature_array& f1, const feature_array& f2) -- exit";
  return result;
}


const std::string& ElementWiseSquaredDistanceCalculator::name() const {
  return name_;
}


const std::string ElementWiseSquaredDistanceCalculator::name_ = "squared distance";

} // namespace feature_extraction

} // namespace pgmlink
