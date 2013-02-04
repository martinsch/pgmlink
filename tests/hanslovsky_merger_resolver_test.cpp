#define BOOST_TEST_MODULE hanslovsky_merger_resolver_test

#include <stdexcept>
#include <cstring>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include "pgmlink/hypotheses.h"
#include "pgmlink/traxels.h"
#include "pgmlink/hanslovsky.h"

#include <armadillo>
#include <mlpack/core.hpp>
#include <mlpack/methods/kmeans/kmeans.hpp>

using namespace pgmlink;
using namespace std;
using namespace boost;


BOOST_AUTO_TEST_CASE( MergerResolver_constructor ) {
  HypothesesGraph g;
  MergerResolver m(&g);
  BOOST_CHECK_EQUAL(m.giveGraph(), &g);
}


BOOST_AUTO_TEST_CASE( MergerResolver_resolve_mergers ) {
  HypothesesGraph g;
  MergerResolver m(&g);
}


BOOST_AUTO_TEST_CASE( MergerResolver_kmeans) {
  // KMeans(int k, const feature_array& data, feature_array& centers
  // KMeans::predict()
  float arr[] = {-6, 0, 0, -5, 0, 0, -4, 0, 0, 6, 0, 0, 5, 0, 0, 4, 0, 0};
  float arr_res[] = {5, 0, 0, -5, 0, 0};
  feature_array data(arr, arr + sizeof(arr)/sizeof(arr[0]));
  KMeans kMeans(2, data);
  feature_array centers = kMeans();
  BOOST_CHECK_EQUAL_COLLECTIONS(centers.begin(), centers.end(), arr_res, arr_res+sizeof(arr_res)/sizeof(arr_res[0]));
}


BOOST_AUTO_TEST_CASE( MergerResolver_helper_functions_vector_to_mat ) {
  // void feature_array_to_arma_mat(const feature_array& in, arma::mat& out);
  
  float arr[] = {1, 2, 0, 3, 1, 2};
  feature_array ft(arr, arr + sizeof(arr)/sizeof(arr[0]));
  arma::fmat m1(6,1);
  arma::fmat m2(3,2);
  arma::fmat m3(6,2);
  feature_array::iterator ftIt;

  feature_array_to_arma_mat(ft, m1);
  ftIt = ft.begin();
  for (unsigned c = 0; c < m1.n_cols; ++c, ftIt += m1.n_rows) {
    arma::Col<float> vc = m1.col(c);
    BOOST_CHECK_EQUAL_COLLECTIONS(ftIt, ftIt + m1.n_rows, vc.begin(), vc.end());
  }
				
  feature_array_to_arma_mat(ft, m2);
  ftIt = ft.begin();
  for (unsigned c = 0; c < m2.n_cols; ++c, ftIt += m2.n_rows) {
    arma::Col<float> vc = m2.col(c);
    BOOST_CHECK_EQUAL_COLLECTIONS(ftIt, ftIt + m2.n_rows, vc.begin(), vc.end());
  }

  BOOST_CHECK_THROW(feature_array_to_arma_mat(ft, m3), std::range_error);
}


BOOST_AUTO_TEST_CASE( MergerResolver_helper_functions_get_centers ) {
  // void get_centers(const arma::mat& data, const arma::Col<size_t> labels, arma::mat& centers, int k);
  
  float arr[] = {1, 2, 0, 3, 1, 2, -1, -2, 0, 1, 1, 1};
  size_t ass[] = {0, 2, 0, 1};
  std::vector<std::vector<float> > c;
  for (int i = 0; i < 3; ++i) c.push_back(std::vector<float>(3));
  c[0][0] = 0; c[0][1] = 0; c[0][2] = 0;
  c[1][0] = 1; c[1][1] = 1; c[1][2] = 1;
  c[2][0] = 3; c[2][1] = 1; c[2][2] = 2;
  int k = 3;
  feature_array ft(arr, arr + sizeof(arr)/sizeof(arr[0]));
  std::vector<size_t> asgn(ass, ass + sizeof(ass)/sizeof(ass[0]));
  arma::mat data(3,4);
  arma::mat centers(3,3);
  feature_array_to_arma_mat(ft, data);
  arma::Col<size_t> assign(asgn);
  get_centers(data, assign, centers, k);
  for (unsigned n = 0; n < centers.n_cols; ++n) {
    arma::Col<double> fvec = centers.col(n);
    BOOST_CHECK_EQUAL_COLLECTIONS(c[n].begin(), c[n].end(), fvec.begin(), fvec.end());
  }
}


// EOF
