/**
* @file
* @ingroup tracking
* @brief outlier detection
*/
#ifndef PGMLINK_OUTLIER_DETECTION_H
#define PGMLINK_OUTLIER_DETECTION_H

#include <vector>
#include <stdexcept>
#include <cassert>

#include <armadillo>

// pgmlink
#include "pgmlink/tracks.h"
#include "pgmlink/traxels.h"
#include "pgmlink/classifier_auxiliary.h"

// boost
#include <boost/shared_ptr.hpp>

namespace pgmlink {

class OutlierCalculator {
  public:
    OutlierCalculator() {};
    ~OutlierCalculator() {};
    virtual const std::vector<size_t>& calculate(const feature_arrays& features) = 0;
    virtual const feature_array& get_measures() const;
    virtual const std::string& name() const;
  protected:
    static const std::string name_;
    feature_array measures_;
}; // class OutlierCalculator

class MVNOutlierCalculator : public OutlierCalculator {
  public:
    MVNOutlierCalculator() {};
    ~MVNOutlierCalculator() {};
    const std::vector<size_t>& calculate(const feature_arrays& features);
    const feature_array& get_measures() const;
    const arma::Mat<feature_type>& get_covariance() const;
    const arma::Mat<feature_type>& get_inverse_covariance() const;
    const arma::Col<feature_type>& get_mean() const;
    virtual const std::string& name() const;

  protected:
    static const std::string name_;
    static const feature_type sigma_threshold_;
    feature_array measures_;
    std::vector<size_t> outlier_ids_;
    arma::Col<feature_type> mean_;
    arma::Mat<feature_type> covariance_;
    arma::Mat<feature_type> inv_covariance_;
}; // class MVNOutlierCalculator

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
