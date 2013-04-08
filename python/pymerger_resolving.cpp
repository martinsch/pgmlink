#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

#include <boost/python.hpp>
#include "pgmlink/merger_resolving.h"

#include <vigra/multi_array.hxx>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>


using namespace std;
using namespace pgmlink;
using namespace boost::python;

// feature dimension == number of rows in arma::mat!!!
// -> shape = [n_features, n_samples]

template<typename T>
void gmm_priors_and_centers_numpy_to_arma(const vigra::NumpyArray<2, T>& data, feature_array& priors, feature_array& centers, int k_max, int ndim, double regularization_weight) {
  arma::mat d(data.shape()[0], data.shape()[1]);
  assert((unsigned)ndim == d.n_rows);
  arma::mat::iterator dest_it = d.begin();
  typename vigra::NumpyArray<2, T>::const_iterator src_it = data.begin();
  for (; src_it != data.end(); ++src_it, ++dest_it) {
    *dest_it = *src_it;
  }
  gmm_priors_and_centers_arma(d, priors, centers, k_max, ndim, regularization_weight);
}

void export_gmm() {
  def("gmm_priors_and_centers", gmm_priors_and_centers);
  def("gmm_priors_and_centers", gmm_priors_and_centers_numpy_to_arma<double>);
  def("gmm_priors_and_centers", gmm_priors_and_centers_numpy_to_arma<unsigned char>);
  def("gmm_priors_and_centers", gmm_priors_and_centers_numpy_to_arma<float>);
  def("gmm_priors_and_centers", gmm_priors_and_centers_numpy_to_arma<unsigned>);
  def("gmm_priors_and_centers", gmm_priors_and_centers_numpy_to_arma<int>);
}
