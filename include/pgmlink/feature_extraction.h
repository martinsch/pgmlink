#ifndef FEATURE_EXTRACTION_STRATEGY_H
#define FEATURE_EXTRACTION_STRATEGY_H

// stl
#include <string>
#include <map>
#include <utility>

// boost
#include <boost/shared_ptr.hpp>

// pgmlink
#include "feature.h"
#include "pgmlink/pgmlink_export.h"

namespace pgmlink {

namespace feature_extraction {

// forward declaration of FeatureCalculator
class FeatureCalculator;

////
//// class FeatureExtractor
////
class FeatureExtractor 
{
 public:
  PGMLINK_EXPORT FeatureExtractor(boost::shared_ptr<FeatureCalculator> calculator, const std::string& feature_name);
  PGMLINK_EXPORT virtual ~FeatureExtractor();
  PGMLINK_EXPORT virtual feature_array extract(const Traxel& t1) const;
  PGMLINK_EXPORT virtual feature_array extract(const Traxel& t1, const Traxel& t2) const;
  PGMLINK_EXPORT virtual feature_array extract(const Traxel& t1, const Traxel& t2, const Traxel& t3) const;
  PGMLINK_EXPORT boost::shared_ptr<FeatureCalculator> calculator() const;
  PGMLINK_EXPORT virtual std::string name() const;

 protected:
  boost::shared_ptr<FeatureCalculator> calculator_;
  std::string feature_name_;
};


////
//// class MultipleFeaturesExtraction
////
class MultipleFeatureExtraction 
{
 public:
  typedef std::pair<std::string, std::string> CombinedFeatureName;
  typedef std::map<CombinedFeatureName, feature_array> CombinedFeatureMap;
  typedef std::map< std::string, std::vector< std::string > > FeatureList;

  PGMLINK_EXPORT CombinedFeatureMap operator() ( const FeatureList& features, const Traxel& trax ) const;
  PGMLINK_EXPORT CombinedFeatureMap operator() ( const FeatureList& features, const Traxel& trax1, const Traxel& trax2 ) const;
  PGMLINK_EXPORT CombinedFeatureMap operator() ( const FeatureList& features, const Traxel& trax1, const Traxel& trax2, const Traxel& trax3 ) const;
};



////
//// helpers
////
namespace helpers 
{

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
