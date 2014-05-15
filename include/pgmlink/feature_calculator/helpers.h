#ifndef FEATURE_CALCULATOR_HELPERS_H
#define FEATURE_CALCULATOR_HELPERS_H

// stl
#include <string>
#include <map>

// boost
#include <boost/shared_ptr.hpp>

namespace pgmlink {


namespace feature_extraction {

class FeatureCalculator;


namespace helpers {
////
//// class CalculatorLookup
////
class CalculatorLookup {
 public:
  static boost::shared_ptr<FeatureCalculator> extract_calculator(const std::string& name);
  
 private:
  static const std::map<std::string, boost::shared_ptr<FeatureCalculator> > calculator_map_;
};

} // namespace helpers

} // namespace feature_extraction

} // namespace pgmlink


#endif /* FEATURE_CALCULATOR_HELPERS_H */

