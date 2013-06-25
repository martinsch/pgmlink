// stl headers
#include <vector>

//boost
#include <boost/shared_ptr.hpp>

// armadillo
#include <armadillo>

// mlpack
#include <mlpack/core.hpp>
#include <mlpack/methods/kmeans/kmeans.hpp>
#include <mlpack/methods/gmm/gmm.hpp>

// pgmlink headers
#include "pgmlink/traxels.h"
#include "pgmlink/clustering.h"
#include "pgmlink/log.h"


namespace pgmlink {
  ////
  //// ClusteringMlpackBase
  ////
  void ClusteringMlpackBase::copy_centers_to_feature_array(const arma::mat& centers, feature_array& c) {
    int n = centers.n_cols;
    int stepSize = centers.n_rows;
    
    if (stepSize*n != (int)c.size()) {
      throw std::range_error("Source matrix dimensions and vector dimension do not agree!");
    }

    feature_array::iterator it = c.begin();
    for (int i = 0; i < n; ++i, it += stepSize) {
      arma::vec col = centers.col(i);
      std::copy(col.begin(), col.end(), it);
    }
  }


  ////
  //// KMeans
  ////
  feature_array KMeans::operator()() {
    mlpack::kmeans::KMeans<> kMeans;
    int n = data_.size()/3;
    arma::mat data(3, n);
    arma::mat centers(3, k_);
    arma::Col<size_t> labels;
    feature_array_to_arma_mat(data_, data);
    kMeans.Cluster(data, k_, labels);
    get_centers(data, labels, centers, k_);
    feature_array fa_centers(3*k_);
    copy_centers_to_feature_array(centers, fa_centers);
    return fa_centers;
  }


  ////
  //// GMM
  ////
  feature_array GMM::operator()() {
    mlpack::gmm::GMM<> gmm(k_, n_);
    int n_samples = data_.size()/3;
    arma::mat data(n_,n_samples);
    arma::Col<size_t> labels;
    LOG(logDEBUG1) << "GMM::operator(): n_=" << n_;
    if (n_ == 2) {
      feature_array_to_arma_mat_skip_last_dimension(data_, data, 3);
    } else if(n_ == 3) {
      feature_array_to_arma_mat(data_, data);
    } else {
      throw std::runtime_error("Number of spatial dimensions other than 2 or 3 would not make sense!");
    }
    score_ = gmm.Estimate(data, n_trials_);
    std::vector<arma::vec> centers = gmm.Means();
    feature_array fa_centers;
    for (std::vector<arma::vec>::iterator it = centers.begin(); it != centers.end(); ++it) {
      std::copy(it->begin(), it->end(), std::back_insert_iterator<feature_array >(fa_centers));
      if (n_ == 2) {
        fa_centers.push_back(0);
      }
    }
    return fa_centers;
  }


  /**
   * returns n*log_likelihood of the model
   */
  double GMM::score() const {
    return score_;
  }
  ////
  //// GMMInitalizeArma
  ////

  feature_array GMMInitializeArma::operator()() {
    feature_array ret;
    // std::back_insert_iterator<feature_array > bii(ret);
    // arma::vec d = this->operator()("this is so dirty!");
    // std::copy(d.begin(), d.end(), bii);
    return ret;
  }
  

  std::vector<arma::vec> GMMInitializeArma::operator()(const char*) {
    mlpack::gmm::GMM<> gmm(k_, data_.n_rows);
    score_ = gmm.Estimate(data_, n_trials_);
    std::vector<arma::vec> centers = gmm.Means();
    return centers;
  }

  
  double GMMInitializeArma::score() const {
    return score_;
  }





  /*void KMeans::copy_centers_to_feature_array(const arma::mat& centers, feature_array& c) {
    int n = centers.n_cols;
    int stepSize = centers.n_rows;
    
    if (stepSize*n != (int)c.size()) {
      throw std::range_error("Source matrix dimensions and vector dimension do not agree!");
    }

    feature_array::iterator it = c.begin();
    for (int i = 0; i < n; ++i, it += stepSize) {
      arma::vec col = centers.col(i);
      std::copy(col.begin(), col.end(), it);
    }
  }*/
}
