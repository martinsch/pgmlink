/**
* @file
* @ingroup tracking
* @brief outlier detection
*/
#ifndef PGMLINK_OUTLIER_DETECTION_H
#define PGMLINK_OUTLIER_DETECTION_H

#include <vector>

// pgmlink
#include "pgmlink/tracks.h"
#include "pgmlink/traxels.h"
#include "pgmlink/classifier_auxiliary.h"

// boost
#include <boost/shared_ptr.hpp>

namespace pgmlink {

class Outlier {
  public:
    Outlier(
      boost::shared_ptr<TrackFeatureExtractor> extractor,
      boost::shared_ptr<OutlierCalculator> outlier_calculator
        = boost::shared_ptr<OutlierCalculator>(new MVNOutlierCalculator())
    );
    virtual const std::vector<size_t> calculate(const Track& track);
    virtual const feature_array& get_measures() const;
    boost::shared_ptr<TrackFeatureExtractor> extractor_;
    boost::shared_ptr<OutlierCalculator> outlier_calculator_;
  protected:
    feature_array measures_;
}; // class Outlier

} // Namespace pgmlink
#endif
