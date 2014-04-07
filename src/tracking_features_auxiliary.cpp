#include <cmath>

#include "pgmlink/tracking_features_auxiliary.h"

namespace pgmlink{
//
// function to_arma_mat
//
arma::Mat<feature_type> to_arma_matrix(const feature_arrays& features) {
  assert(features.empty() == false);
  size_t cols = features.size();
  size_t rows = features[0].size();
  arma::Mat<feature_type> ret(rows, cols);

  typename std::vector<feature_array>::const_iterator feature_array_it;
  feature_array_it = features.begin();
  for (
    size_t j=0;
    feature_array_it != features.end();
    feature_array_it++, j++
  ) {
    assert(feature_array_it->size() == rows);
    arma::Col<feature_type> column(*feature_array_it);
    ret.col(j) = column;
  }
  return ret;
}

//
// Class OutlierCalculator
//
const std::string OutlierCalculator::name_ = "OutlierCalculator";

const std::string& OutlierCalculator::name() const {
  return name_;
}

const feature_array& OutlierCalculator::get_measures() const {
  throw std::runtime_error(
    "OutlierCalculator \"" + name() + "\"doesn't provide a measure"
  );
  return measures_;
}

//
// Class MVNOutlierCalculator
//
const std::string MVNOutlierCalculator::name_ = "MVNOutlierCalculator";

MVNOutlierCalculator::MVNOutlierCalculator(
  const feature_type sigma_threshold
) {
  MVNOutlierCalculator::sigma_threshold_ = sigma_threshold;
}

const std::string& MVNOutlierCalculator::name() const {
  return name_;
}

const feature_array& MVNOutlierCalculator::get_measures() const {
  return measures_;
}

const arma::Mat<feature_type>& MVNOutlierCalculator::get_covariance() const {
  return covariance_;
}

const arma::Mat<feature_type>& MVNOutlierCalculator::get_inverse_covariance() const {
  return inv_covariance_;
}

const arma::Col<feature_type>& MVNOutlierCalculator::get_mean() const {
  return mean_;
}

const std::vector<size_t>& MVNOutlierCalculator::calculate(
  const feature_arrays& features
) {
  if (features.size() <= features[0].size()) {
    measures_.clear();
    outlier_ids_.clear();
    mean_.clear();
    covariance_.clear();
    inv_covariance_.clear();
    std::cerr << "Too few data to calculate outliers" << std::endl;
  } else {
    // Get covariance and inverse covariance matrix
    arma::Mat<feature_type> features_mat(to_arma_matrix(features));
    arma::Mat<feature_type> features_mat_t(trans(features_mat));
    try {
      covariance_ = arma::cov(features_mat_t);
      inv_covariance_ = arma::inv_sympd(covariance_);
  
      // Get mean values
      mean_ = arma::mean(features_mat, 1);
  
      // Calculate the outliers
      outlier_ids_.clear();
      measures_.clear();
      feature_arrays::const_iterator features_it = features.begin();
      for(size_t id=0; features_it != features.end(); features_it++, id++) {
        arma::Col<feature_type> diff_vector(*features_it);
        diff_vector -= mean_;
        feature_type norm_residual = arma::dot(diff_vector, inv_covariance_*diff_vector);
        measures_.push_back(norm_residual);
        if (norm_residual >= sigma_threshold_) {
          outlier_ids_.push_back(id);
        }
      }
    } catch (std::exception& exception) {
      std::cout << "Too few data to calculate outliers" << std::endl;
    }
  } // else
  return outlier_ids_;
}




////
//// Class FeatureAggregator
////
const std::string FeatureAggregator::name_ = "";

feature_array FeatureAggregator::vector_valued(
  const feature_arrays& features
) {
  (void)features; // Casting to void to avoid the warning "unused parameter"
  throw std::runtime_error(
    "FeatureAggregator \"" + name() + "\" has no vector valued method"
  );
  feature_array ret;
  return ret;
}

feature_type FeatureAggregator::scalar_valued(
  const feature_arrays& features
) {
  (void)features; // Casting to void to avoid the warning "unused parameter"
  throw std::runtime_error(
    "FeatureAggregator \"" + name() + "\" has no scalar valued method"
  );
  feature_type ret;
  return ret;
}

feature_type FeatureAggregator::scalar_valued(
  const feature_array& features
) {
  feature_arrays f(features.size());
  feature_arrays::iterator f_it = f.begin();
  feature_array::const_iterator features_it = features.begin();
  for(; features_it != features.end(); features_it++, f_it++) {
    f_it->push_back(*features_it);
  }
  return scalar_valued(f);
}

const std::string& FeatureAggregator::name() const {
  return FeatureAggregator::name_;
}

////
//// OutlierBadnessAggregator
////
const std::string OutlierBadnessAggregator::name_ = "OutlierBadness";

const std::string& OutlierBadnessAggregator::name() const {
  return OutlierBadnessAggregator::name_;
}

feature_array OutlierBadnessAggregator::vector_valued(
  const feature_arrays& features
) {
  mvn_outlier_calculator_.calculate(features);
  return mvn_outlier_calculator_.get_measures();
}

feature_type OutlierBadnessAggregator::scalar_valued(
  const feature_arrays& features
) {
  feature_array vector = vector_valued(features);
  feature_type max = 0;
  for (feature_array::iterator it = vector.begin(); it!= vector.end(); it++) {
    max = max > *it ? max : *it;
  }
  return max;
}

////
//// OutlierCountAggregator
////
const std::string OutlierCountAggregator::name_ = "OutlierCount";

const std::string& OutlierCountAggregator::name() const {
  return OutlierCountAggregator::name_;
}

feature_type OutlierCountAggregator::scalar_valued(
  const feature_arrays& features
) {
  std::vector<size_t> out_ids = mvn_outlier_calculator_.calculate(features);
  return feature_type(out_ids.size());
}

////
//// SumAggregator
////
const std::string SumAggregator::name_ = "Sum";

const std::string& SumAggregator::name() const {
  return SumAggregator::name_;
}

feature_array SumAggregator::vector_valued(
  const feature_arrays& features
) {
  size_t n = features.size();
  (void)n; // Avoid variable unused warning
  assert(n > 0);
  size_t feature_dim = features.front().size();
  (void)feature_dim; // Avoid variable unused warning
  assert(feature_dim > 0);

  feature_array ret(feature_dim, 0.0);
  feature_arrays::const_iterator it_fs;
  for(it_fs = features.begin(); it_fs != features.end(); it_fs++) {
    assert(feature_dim == it_fs->size());
    feature_array::iterator it_ret = ret.begin();
    feature_array::const_iterator it_f = it_fs->begin();
    for(; it_f != it_fs->end(); it_f++, it_ret++) {
      (*it_ret) += (*it_f);
    }
  }

  return ret;
}

feature_type SumAggregator::scalar_valued(
  const feature_arrays& features
) {
  feature_type ret = 0;
  feature_array vector = vector_valued(features);
  feature_array::iterator v_it;
  for(v_it = vector.begin(); v_it != vector.end(); v_it++) {
    ret += (*v_it);
  }

  return ret;
}
////
//// TotalDiffAggregator
////
const std::string TotalDiffAggregator::name_ = "TotalDiff";

const std::string& TotalDiffAggregator::name() const {
  return TotalDiffAggregator::name_;
}

feature_array TotalDiffAggregator::vector_valued(
  const feature_arrays& features
) {
  assert(features.size() > 0);
  feature_array first = features.front();
  feature_array last = features.back();
  size_t feature_dim = first.size();
  assert(feature_dim > 0);
  assert(feature_dim == last.size());

  feature_array ret(feature_dim);
  feature_array::const_iterator it_first = first.begin();
  feature_array::const_iterator it_last = last.begin();
  feature_array::iterator it_ret = ret.begin();
  for(; it_first != first.end(); it_first++, it_last++, it_ret++) {
    *it_ret = *it_last - *it_first;
  }

  return ret;
}

feature_type TotalDiffAggregator::scalar_valued(
  const feature_arrays& features
) {
  feature_array vector = vector_valued(features);
  feature_type ret = 0;
  if (vector.size() == 1) {
    ret = vector.front();
  } else {
    for(feature_array::iterator it = vector.begin(); it != vector.end(); it++) {
      ret += (*it) * (*it);
    }
    ret = sqrt(ret);
  }
  return ret;
}

////
//// MinAggregator
////
const std::string MinAggregator::name_ = "Min";

const std::string& MinAggregator::name() const {
  return MinAggregator::name_;
}

feature_array MinAggregator::vector_valued(
  const feature_arrays& features
) {
  assert(features.size() > 0);
  size_t feature_dim = features.front().size();
  (void)feature_dim; // Avoid variable unused warning
  assert(feature_dim > 0);

  feature_array ret(features.front());
  feature_arrays::const_iterator it_fs;
  for(it_fs = features.begin()+1; it_fs != features.end(); it_fs++) {
    assert(feature_dim == it_fs->size());
    feature_array::iterator it_ret = ret.begin();
    feature_array::const_iterator it_f = it_fs->begin();
    for(; it_f != it_fs->end(); it_f++, it_ret++) {
      (*it_ret) = (*it_ret) < (*it_f) ? (*it_ret) : (*it_f);
    }
  }

  return ret;
}

feature_type MinAggregator::scalar_valued(
  const feature_arrays& features
) {
  assert(features.size() > 0);
  size_t feature_dim = features.front().size();
  (void)feature_dim; // Avoid variable unused warning
  assert(feature_dim > 0);

  feature_type ret = features[0][0];
  feature_arrays::const_iterator it_fs;
  for(it_fs = features.begin()+1; it_fs != features.end(); it_fs++) {
    assert(feature_dim == it_fs->size());
    feature_array::const_iterator it_f = it_fs->begin();
    for(; it_f != it_fs->end(); it_f++) {
      ret = ret < (*it_f) ? ret : (*it_f);
    }
  }

  return ret;
}

////
//// MaxAggregator
////
const std::string MaxAggregator::name_ = "Max";

const std::string& MaxAggregator::name() const {
  return MaxAggregator::name_;
}

feature_array MaxAggregator::vector_valued(
  const feature_arrays& features
) {
  assert(features.size() > 0);
  size_t feature_dim = features.front().size();
  (void)feature_dim; // Avoid variable unused warning
  assert(feature_dim > 0);

  feature_array ret(features.front());
  feature_arrays::const_iterator it_fs;
  for(it_fs = features.begin()+1; it_fs != features.end(); it_fs++) {
    assert(feature_dim == it_fs->size());
    feature_array::iterator it_ret = ret.begin();
    feature_array::const_iterator it_f = it_fs->begin();
    for(; it_f != it_fs->end(); it_f++, it_ret++) {
      (*it_ret) = (*it_ret) > (*it_f) ? (*it_ret) : (*it_f);
    }
  }

  return ret;
}

feature_type MaxAggregator::scalar_valued(
  const feature_arrays& features
) {
  assert(features.size() > 0);
  size_t feature_dim = features.front().size();
  (void)feature_dim; // Avoid variable unused warning
  assert(feature_dim > 0);

  feature_type ret = features[0][0];
  feature_arrays::const_iterator it_fs;
  for(it_fs = features.begin()+1; it_fs != features.end(); it_fs++) {
    assert(feature_dim == it_fs->size());
    feature_array::const_iterator it_f = it_fs->begin();
    for(; it_f != it_fs->end(); it_f++) {
      ret = ret > (*it_f) ? ret : (*it_f);
    }
  }

  return ret;
}

////
//// MeanAggregator
////
const std::string MeanAggregator::name_ = "Mean";

const std::string& MeanAggregator::name() const {
  return MeanAggregator::name_;
}

feature_array MeanAggregator::vector_valued(
  const feature_arrays& features
) {
  size_t n = features.size();
  assert(n > 0);
  size_t feature_dim = features.front().size();
  (void)feature_dim; // Avoid variable unused warning
  assert(feature_dim > 0);

  feature_array ret(feature_dim, 0.0);
  feature_arrays::const_iterator it_fs;
  for(it_fs = features.begin(); it_fs != features.end(); it_fs++) {
    assert(feature_dim == it_fs->size());
    feature_array::iterator it_ret = ret.begin();
    feature_array::const_iterator it_f = it_fs->begin();
    for(; it_f != it_fs->end(); it_f++, it_ret++) {
      (*it_ret) += (*it_f) / feature_type(n);
    }
  }

  return ret;
}

feature_type MeanAggregator::scalar_valued(
  const feature_arrays& features
) {
  feature_type ret = 0;
  feature_array vector = vector_valued(features);
  size_t n = vector.size();
  feature_array::iterator v_it;
  for(v_it = vector.begin(); v_it != vector.end(); v_it++) {
    ret += (*v_it) / feature_type(n);
  }

  return ret;
}

} // namespace pgmlink
