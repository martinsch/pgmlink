#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

#include <boost/python.hpp>
#include "pgmlink/merger_resolving.h"

#include <vigra/multi_array.hxx>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/tinyvector.hxx>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

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



class PyTimestepIdCoordinateMap
{
 public:
  PyTimestepIdCoordinateMap() : map_(TimestepIdCoordinateMapPtr(new TimestepIdCoordinateMap)) {}
  void initialize() {if (!map_) map_ = TimestepIdCoordinateMapPtr(new TimestepIdCoordinateMap);}
  TimestepIdCoordinateMapPtr get() {return map_;}
  void set(TimestepIdCoordinateMapPtr ptr){ map_ = ptr; }
  size_t size() const
  {
      if(!map_)
      {
          return 0;
      }
      return map_->size();
  }
 private:
  TimestepIdCoordinateMapPtr map_;
};

template <int N, typename T>
void py_extract_coordinates(PyTimestepIdCoordinateMap coordinates,
                            const vigra::NumpyArray<N, T>& image,
                            const vigra::NumpyArray<1, vigra::Int64>& offsets,
                            const Traxel& trax) {
  if (offsets.shape()[0] != N) {
    throw std::runtime_error("py_extract_coordinates() -- Number of offsets and image dimensions disagree!");
  }
  vigra::TinyVector<long int, N> offsets_tv;
  for (size_t idx = 0; idx < N; ++idx) {
    offsets_tv[idx] = offsets[idx];
  }
  extract_coordinates<N, T>(coordinates.get(), image, offsets_tv, trax);
}

template <int N, typename T>
void py_extract_coord_by_timestep_id(PyTimestepIdCoordinateMap coordinates,
                            const vigra::NumpyArray<N, T>& image,
                            const vigra::NumpyArray<1, vigra::Int64>& offsets,
                            const size_t timestep,
                            const size_t traxel_id,
                            const size_t traxel_size) {
  if (offsets.shape()[0] != N) {
    throw std::runtime_error("py_extract_coord_by_timestep_id() -- Number of offsets and image dimensions disagree!");
  }
  vigra::TinyVector<long int, N> offsets_tv;
  for (size_t idx = 0; idx < N; ++idx) {
    offsets_tv[idx] = offsets[idx];
  }
  extract_coord_by_timestep_id<N, T>(coordinates.get(),
                                     image,
                                     offsets_tv,
                                     timestep,
                                     traxel_id,
                                     traxel_size);
}

template <int N, typename T>
void py_update_labelimage(PyTimestepIdCoordinateMap coordinates,
                          vigra::NumpyArray<N, T> image,
                          vigra::NumpyArray<1, vigra::Int64> offsets,
                          const size_t timestep,
                          const size_t traxel_id) {
  vigra::TinyVector<long int, N> offsets_tv;
  for (size_t idx = 0; idx < N; ++idx) {
    offsets_tv[idx] = offsets[idx];
  }
  update_labelimage<N, T>(coordinates.get(), image, offsets_tv, timestep, traxel_id);
}

class CoordinateMapPickleSuite : public boost::python::pickle_suite
{
public:
    static std::string getstate(PyTimestepIdCoordinateMap& coordinates)
    {
        std::stringstream ss;
        boost::archive::binary_oarchive oa(ss);
        TimestepIdCoordinateMapPtr map_ptr = coordinates.get();
        oa & map_ptr;
        return ss.str();
    }

    static void setstate(PyTimestepIdCoordinateMap& coordinates, const std::string& state)
    {
        std::stringstream ss(state);
        boost::archive::binary_iarchive ia(ss);
        TimestepIdCoordinateMapPtr map_ptr;
        ia & map_ptr;
        coordinates.set(map_ptr);
    }
};


void export_gmm() {
  def("gmm_priors_and_centers", gmm_priors_and_centers);
  def("gmm_priors_and_centers", gmm_priors_and_centers_numpy_to_arma<double>);
  def("gmm_priors_and_centers", gmm_priors_and_centers_numpy_to_arma<unsigned char>);
  def("gmm_priors_and_centers", gmm_priors_and_centers_numpy_to_arma<float>);
  def("gmm_priors_and_centers", gmm_priors_and_centers_numpy_to_arma<unsigned>);
  def("gmm_priors_and_centers", gmm_priors_and_centers_numpy_to_arma<int>);

  class_<PyTimestepIdCoordinateMap>("TimestepIdCoordinateMap")
      .def("initialize", &PyTimestepIdCoordinateMap::initialize)
      .def("get", &PyTimestepIdCoordinateMap::get)
      .def("size", &PyTimestepIdCoordinateMap::size)
      .def_pickle(CoordinateMapPickleSuite())
      ;

  class_<TimestepIdCoordinateMapPtr>("TimestepIdCoordinateMapPtr");

  def("extract_coordinates", vigra::registerConverters(&py_extract_coordinates<2, vigra::UInt8>));
  def("extract_coordinates", vigra::registerConverters(&py_extract_coordinates<3, vigra::UInt8>));

  def("extract_coordinates", vigra::registerConverters(&py_extract_coordinates<2, vigra::UInt16>));
  def("extract_coordinates", vigra::registerConverters(&py_extract_coordinates<3, vigra::UInt16>));

  def("extract_coordinates", vigra::registerConverters(&py_extract_coordinates<2, vigra::UInt32>));
  def("extract_coordinates", vigra::registerConverters(&py_extract_coordinates<3, vigra::UInt32>));

  def("extract_coord_by_timestep_id", vigra::registerConverters(&py_extract_coord_by_timestep_id<2, vigra::UInt8>));
  def("extract_coord_by_timestep_id", vigra::registerConverters(&py_extract_coord_by_timestep_id<3, vigra::UInt8>));

  def("extract_coord_by_timestep_id", vigra::registerConverters(&py_extract_coord_by_timestep_id<2, vigra::UInt16>));
  def("extract_coord_by_timestep_id", vigra::registerConverters(&py_extract_coord_by_timestep_id<3, vigra::UInt16>));

  def("extract_coord_by_timestep_id", vigra::registerConverters(&py_extract_coord_by_timestep_id<2, vigra::UInt32>));
  def("extract_coord_by_timestep_id", vigra::registerConverters(&py_extract_coord_by_timestep_id<3, vigra::UInt32>));

  def("update_labelimage", vigra::registerConverters(&py_update_labelimage<2, vigra::UInt8>));
  def("update_labelimage", vigra::registerConverters(&py_update_labelimage<3, vigra::UInt8>));

  def("update_labelimage", vigra::registerConverters(&py_update_labelimage<2, vigra::UInt16>));
  def("update_labelimage", vigra::registerConverters(&py_update_labelimage<3, vigra::UInt16>));

  def("update_labelimage", vigra::registerConverters(&py_update_labelimage<2, vigra::UInt32>));
  def("update_labelimage", vigra::registerConverters(&py_update_labelimage<3, vigra::UInt32>));
}
