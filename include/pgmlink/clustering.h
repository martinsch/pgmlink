#ifndef CLUSTERING_H
#define CLUSTERING_H


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



namespace pgmlink {
  class ClusteringMlpackBase;

  
  class ClusteringMlpackBuilderBase;

  
  typedef boost::shared_ptr<ClusteringMlpackBase> ClusteringPtr;

  typedef boost::shared_ptr<ClusteringMlpackBuilderBase> ClusteringBuilderPtr;
  ////
  //// ClusteringMlpackBase
  ////
  class ClusteringMlpackBase {
  private:
  protected:
    void copy_centers_to_feature_array(const arma::mat& centers, feature_array& c);
  public:
    virtual ~ClusteringMlpackBase() {}
    virtual feature_array operator()() = 0;
    virtual double score() const {return 0.0;}
    virtual const arma::mat& get_data_arma() const = 0;
    virtual unsigned get_cluster_assignment(const arma::vec& sample) {return 0;}
    virtual void set_k_clusters(unsigned k) {return;}
  };


  ////
  //// KMeans
  ////
  /**
   * @class KMeans 
   * @brief compatibility class for kMeans as an interface between feature_array and the mlpack library used for
   * kMeans
   *
   *
   * The library mlpack provides several clustering algorithms, one of which is kMeans. Instead of doing repetetive work
   * as rewriting kMeans would be, we are using the kMeans implementation provided in mlpack. To do so we need to convert
   * the data given in the form of a feature_array into an appropriate armadillo matrix (arma::mat), that can be used by
   * mlpack
   */
  class KMeans : public ClusteringMlpackBase {
  private:
    KMeans();
    int k_;
    const feature_array& data_;
    arma::mat data_arma_;
    // void copy_centers_to_feature_array(const arma::mat& centers, feature_array& c);
  public:
    // tested
    /**
     * @brief Constructor
     * @param [in] k number of clusters
     * @param [in] data feature_array storing data
     */
    KMeans(int k, const feature_array& data) :
      k_(k), data_(data), data_arma_(arma::mat()) {}

    // tested
    /**
     * @brief compute cluster centers and labels for datapoints
     * @returns feature_array that contains the coordinates of k clusters
     */
    virtual feature_array operator()();
    virtual const arma::mat& get_data_arma() const;
  };


  ////
  //// GMM
  ////
  class GMM : public ClusteringMlpackBase {
  private:
    GMM();
    int k_;
    int n_;
    const feature_array& data_;
    double score_;
    int n_trials_;
    arma::mat data_arma_;
    mlpack::gmm::GMM<> gmm_;
  public:
    // constructor needs to specify number of dimensions
    // for 2D data, ilastik provides coordinates with 3rd dimension 0
    // which will cause singular covariance matrix
    // therefore add option for dimensionality
    GMM(int k, int n, const feature_array& data, int n_trials=1);

    virtual feature_array operator()();
    virtual unsigned get_cluster_assignment(const arma::vec& sample);
    double score() const;
    virtual void set_k_clusters(unsigned k);
    virtual const arma::mat& get_data_arma() const;
  };

  ////
  //// GMMInitializeArma
  ////
  class GMMInitializeArma : public ClusteringMlpackBase {
  private:
    GMMInitializeArma();
    int k_;
    arma::mat data_arma_;
    double score_;
    int n_trials_;
  public:
    GMMInitializeArma(int k, const arma::mat& data, int n_trials=1) :
      k_(k), data_arma_(data), score_(0.0), n_trials_(n_trials) {}

    virtual feature_array operator()();
    std::vector<arma::vec> operator()(const char* dirty_hack);
    double score() const;
    virtual const arma::mat& get_data_arma() const;
  };

  ////
  //// ClusteringMlpackBuilderBase
  ////
  class ClusteringMlpackBuilderBase {
  private:
  public:
    virtual ClusteringPtr build(const feature_array& data) = 0;
  };

  ////
  //// KMeansBuilder
  ////
  class KMeansBuilder : public ClusteringMlpackBuilderBase {
  private:
    int k_;
  public:
    KMeansBuilder();
    KMeansBuilder(int k);
    virtual ClusteringPtr build(const feature_array& data);
  };


  ////
  //// GMMBuilder
  ////
  class GMMBuilder : public ClusteringMlpackBuilderBase {
  private:
    int k_;
    int n_;
    int n_trials_;
  public:
    GMMBuilder();
    GMMBuilder(int k, int n, int n_trials);
    virtual ClusteringPtr build(const feature_array& data);
  };
    

  ////
  //// helper functions
  ////
  template <typename T, typename U>
  // tested
  /**
   * @brief Helper function to convert feature_array to arma::Mat.
   * @param [in] in original data; specifying T=float will make in a feature_array
   * @param [in,out] out arma::Mat<U> that holds the converted data. For the use in
   * KMeans specify U=double
   */
  void feature_array_to_arma_mat(const std::vector<T>& in, arma::Mat<U>& out);


  template <typename T, typename U>
  void feature_array_to_arma_mat_skip_last_dimension(const std::vector<T>& in, arma::Mat<U>& out, unsigned int last_dimension);

  
  template <typename T>
  // tested
  /**
   * @brief Helper function to calculate center coordinates from data assignments.
   * @param [in] data data points (coordinates)
   * @param [in] labels assignments after running kMeans
   * @param [in,out] centers arma::Mat to hold the coordinates of the cluster centers
   * @param [in] k number of clusters used for kMeans
   *
   * The mlpack kMeans implementation does not return the coordinates of the cluster centers.
   * The centers can be computed using the original data and the assignments.
   */
  void get_centers(const arma::Mat<T>& data, const arma::Col<size_t> labels, arma::Mat<T>& centers, int k);


  ////
  //// IMPLEMENTATIONS ////
  ////



  template <typename T, typename U>
  void feature_array_to_arma_mat(const std::vector<T>& in, arma::Mat<U>& out) {
    int stepSize = out.n_rows;
    int n = out.n_cols;
    if (stepSize*n != (int)in.size()) {
      throw std::range_error("Source vector dimension and matrix dimensions do not agree!");
    }
    int count = 0;
    typename std::vector<T>::const_iterator srcIt = in.begin();
    while (count < n) {
      arma::Col<U> col(stepSize);
      std::copy(srcIt, srcIt+stepSize, col.begin());
      out.col(count) = col;
      ++count;
      srcIt += stepSize;
    }
  }

  
  template <typename T, typename U>
  void feature_array_to_arma_mat_skip_last_dimension(const std::vector<T>& in, arma::Mat<U>& out, unsigned int last_dimension) {
    unsigned int stepSize = out.n_rows;
    unsigned int n = out.n_cols;
    unsigned int count = 0;
    assert(last_dimension*n == in.size());
    assert(stepSize == last_dimension-1);
    typename std::vector<T>::const_iterator srcIt = in.begin();
    while (count < n) {
      arma::Col<U> col(stepSize);
      std::copy(srcIt, srcIt+stepSize, col.begin());
      out.col(count) = col;
      ++count;
      srcIt += last_dimension;
    }
  }

  
  template <typename T>
  void get_centers(const arma::Mat<T>& data, const arma::Col<size_t> labels, arma::Mat<T>& centers, int k) {
    arma::Col<size_t>::const_iterator labelIt = labels.begin();
    std::vector<int> clusterSize(k, 0);
    centers.zeros();
    for (unsigned int n = 0; n < data.n_cols; ++n, ++labelIt) {
      ++clusterSize[*labelIt];
      centers.col(*labelIt) = centers.col(*labelIt) + data.col(n);
    }
    for (int i = 0; i < k; ++i) {
      centers.col(i) /= clusterSize[i];
    }
  }
}


#endif /* CLUSTERING_H */
