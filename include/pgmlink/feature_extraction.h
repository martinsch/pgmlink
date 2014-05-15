#ifndef FEATURE_EXTRACTION_STRATEGY_H
#define FEATURE_EXTRACTION_STRATEGY_H

// stl
#include <string>
#include <map>
#include <utility>
#include <numeric>

// boost
#include <boost/shared_ptr.hpp>

// pgmlink
#include "feature.h"

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


////
//// IdentityCalculator
////
class IdentityCalculator : public FeatureCalculator {
 public:
  static const std::string name_;

  virtual ~IdentityCalculator();
  virtual feature_array calculate(const feature_array& f1) const;
  virtual const std::string& name() const;

};


////
//// AbsoluteDifferenceCalculator
////
class AbsoluteDifferenceCalculator : public FeatureCalculator {
 public:
  static const std::string name_;

  virtual ~AbsoluteDifferenceCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual const std::string& name() const;
};


////
//// SquareRootSquaredDifferenceCalculator
////
class SquareRootSquaredDifferenceCalculator : public FeatureCalculator {
 public:
  static const std::string name_;

  virtual ~SquareRootSquaredDifferenceCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual const std::string& name() const;
};


////
//// RatioCalculator
////
class RatioCalculator : public FeatureCalculator {
 public:
  static const std::string name_;

  virtual ~RatioCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const;
  virtual const std::string& name() const;
};


////
//// AssymetricRatioCalculator
////
class AsymmetricRatioCalculator : public FeatureCalculator {
 public:
  static const std::string name_;

  virtual ~AsymmetricRatioCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const;
  virtual const std::string& name() const;
};


////
//// SquaredDifferenceCalculator
////
class SquaredDifferenceCalculator : public FeatureCalculator {
 public:
  static const std::string name_;

  virtual ~SquaredDifferenceCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual const std::string& name() const;
};


////
//// class ElementWiseSquaredDistanceCalculator
////
class ElementWiseSquaredDistanceCalculator : public FeatureCalculator {
 public:
  virtual ~ElementWiseSquaredDistanceCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;

  virtual const std::string& name() const;

 private:
  static const std::string name_;
};




////
//// class FeatureExtractor
////
class FeatureExtractor {
 public:
  FeatureExtractor(boost::shared_ptr<FeatureCalculator> calculator, const std::string& feature_name);
  virtual ~FeatureExtractor();
  virtual feature_array extract(const Traxel& t1) const;
  virtual feature_array extract(const Traxel& t1, const Traxel& t2) const;
  virtual feature_array extract(const Traxel& t1, const Traxel& t2, const Traxel& t3) const;
  boost::shared_ptr<FeatureCalculator> calculator() const;
  virtual std::string name() const;

 protected:
  boost::shared_ptr<FeatureCalculator> calculator_;
  std::string feature_name_;
};


////
//// class MultipleFeaturesExtraction
////
class MultipleFeatureExtraction {
 public:
  typedef std::pair<std::string, std::string> CombinedFeatureName;
  typedef std::map<CombinedFeatureName, feature_array> CombinedFeatureMap;
  typedef std::map< std::string, std::vector< std::string > > FeatureList;

  CombinedFeatureMap operator() ( const FeatureList& features, const Traxel& trax ) const;
  CombinedFeatureMap operator() ( const FeatureList& features, const Traxel& trax1, const Traxel& trax2 ) const;
  CombinedFeatureMap operator() ( const FeatureList& features, const Traxel& trax1, const Traxel& trax2, const Traxel& trax3 ) const;
};



////
//// helpers
////
namespace helpers 
{
////
//// class CalculatorLookup
////
class CalculatorLookup {
 public:
  static boost::shared_ptr<FeatureCalculator> extract_calculator(const std::string& name);
  
 private:
  static const std::map<std::string, boost::shared_ptr<FeatureCalculator> > calculator_map_;
};


////
//// function convenience_feature_extraction
////
MultipleFeatureExtraction::CombinedFeatureMap convenience_feature_extraction( const MultipleFeatureExtraction::FeatureList& features,
                                                                              const Traxel& trax );
MultipleFeatureExtraction::CombinedFeatureMap convenience_feature_extraction( const MultipleFeatureExtraction::FeatureList& features,
                                                                              const Traxel& trax1,
                                                                              const Traxel& trax2 );
MultipleFeatureExtraction::CombinedFeatureMap convenience_feature_extraction( const MultipleFeatureExtraction::FeatureList& features,
                                                                              const Traxel& trax1,
                                                                              const Traxel& trax2,
                                                                              const Traxel& trax3 );


} /* namespace helpers */

} /* namespace feature_extraction */

} /* namespace pgmlink */

#endif /* FEATURE_EXTRACTION_STRATEGY_H */
