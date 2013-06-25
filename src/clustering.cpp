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


  const arma::mat& ClusteringMlpackBase::get_data_arma() const {
    return data_arma_;
  }


  ////
  //// KMeans
  ////
  feature_array KMeans::operator()() {
    mlpack::kmeans::KMeans<> kMeans;
    int n = data_.size()/3;
    if (data_arma_.is_empty()) {
      data_arma_ = arma::mat(3, n);
      feature_array_to_arma_mat(data_, data_arma_);
    }
    arma::mat centers(3, k_);
    arma::Col<size_t> labels;
    kMeans.Cluster(data_arma_, k_, labels);
    get_centers(data_arma_, labels, centers, k_);
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
    arma::Col<size_t> labels;
    LOG(logDEBUG1) << "GMM::operator(): n_=" << n_;
    if (data_arma_.is_empty()) {
      data_arma_ = arma::mat(n_,n_samples);
      if (n_ == 2) {
        feature_array_to_arma_mat_skip_last_dimension(data_, data_arma_, 3);
      } else if(n_ == 3) {
        feature_array_to_arma_mat(data_, data_arma_);
      } else {
        throw std::runtime_error("Number of spatial dimensions other than 2 or 3 would not make sense!");
      }
    }
    score_ = gmm.Estimate(data_arma_, n_trials_);
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


  unsigned GMM::get_cluster_assignment(const arma::vec& sample, const uint cluster_id) {
    unsigned res = 0;
    double prob = 0.0;
    for (int idx = 0; idx < k_; ++idx) {
    }
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
    mlpack::gmm::GMM<> gmm(k_, data_arma_.n_rows);
    score_ = gmm.Estimate(data_arma_, n_trials_);
    std::vector<arma::vec> centers = gmm.Means();
    return centers;
  }

  
  double GMMInitializeArma::score() const {
    return score_;
  }


  ////
  //// KMeansBuilder
  ////
  KMeansBuilder::KMeansBuilder() : k_(2) {}

  
  KMeansBuilder::KMeansBuilder(int k) : k_(k) {}


  ClusteringPtr KMeansBuilder::build(const feature_array& data) {
    return ClusteringPtr(new KMeans(k_, data));
  }


  ////
  //// GMMBuilder
  ////
  GMMBuilder::GMMBuilder() : k_(2), n_(3), n_trials_(1) {}

  
  GMMBuilder::GMMBuilder(int k, int n, int n_trials) : k_(k), n_(n), n_trials_(n_trials) {}


  ClusteringPtr GMMBuilder::build(const feature_array& data) {
    return ClusteringPtr(new GMM(k_, n_, data, n_trials_));
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
