#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator/asymmetric_ratio.h"

namespace pgmlink {

namespace feature_extraction {


////
//// class AsymmetricRatioCalculator
////
const std::string AsymmetricRatioCalculator::name_ = "AsymmetricRatio";


AsymmetricRatioCalculator::~AsymmetricRatioCalculator() {

}


feature_array AsymmetricRatioCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  assert(f1.size() == f2.size());
  feature_array ret(f1.size(), 0.);
  // keep ratio <= 1
  // no zero check, as we do not have empty regions
  for (size_t i = 0; i < ret.size(); ++i) {
    if (f1[i] == f2[i]) {
      ret[i] = 1.0;
    } else if (f1[i] < 0.0001 && f2[i] < 0.0001) {
      ret[i] = 1.0;
    } else {
      ret[i] = f1[i]/f2[i];
    }
  }
  return ret;
}


feature_array AsymmetricRatioCalculator::calculate(const feature_array&, const feature_array& f2, const feature_array& f3) const {
  return calculate(f2, f3);
}


const std::string& AsymmetricRatioCalculator::name() const {
  return AsymmetricRatioCalculator::name_;
}


} // namespace feature_extraction

} // namespace pgmlink
