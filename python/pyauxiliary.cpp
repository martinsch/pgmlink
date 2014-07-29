#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

// stl
#include <stdexcept>
#include <string>
#include <iostream>

// boost
#include <boost/python.hpp>
#include <boost/python/return_internal_reference.hpp>

// vigra
#undef tolower
#include <vigra/multi_array.hxx>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

// pgmlink
#include "pgmlink/auxiliary.h"
#include "pgmlink/traxels.h"

using namespace pgmlink;
using namespace boost::python;

struct IntersectIndexPairs {
  std::vector<unsigned> indices;
  std::vector<unsigned> intersects;
};

class IntersectsDictionary {
 public:
  void get_intersects_2d(const vigra::NumpyArray<2, unsigned> arr1,
                      const vigra::NumpyArray<2, unsigned> arr2) {
    extract_intersects<2, unsigned>(arr1, arr2, intersects_);
  }
  
  void get_intersects_3d(const vigra::NumpyArray<3, unsigned> arr1,
                      const vigra::NumpyArray<3, unsigned> arr2) {
    extract_intersects<3, unsigned>(arr1, arr2, intersects_);
  }

  void get_pairs_for(unsigned label, IntersectIndexPairs& res) {
    const std::map<unsigned, unsigned> intersects = intersects_[label];
    for (std::map<unsigned, unsigned>::const_iterator it = intersects.begin();
         it != intersects.end();
         ++it) {
      res.indices.push_back(it->first);
      res.intersects.push_back(it->second);
    }
  }

  void add_to_traxel(Traxel& trax, unsigned label) {
    feature_array& indices = trax.features["intersect_indices"];
    feature_array& intersects = trax.features["intersects"];
    if (indices.size() != intersects.size()) {
      throw std::runtime_error("add_to_traxe(): indices.size() != intersects.size()");
    }
    const std::map<unsigned, unsigned>& intersects_map = intersects_[label];
    for (std::map<unsigned, unsigned>::const_iterator it = intersects_map.begin();
         it != intersects_map.end();
         ++it) {
      indices.push_back(it->first);
      intersects.push_back(it->second);
    }
  }
  
 private:
  std::map<unsigned, std::map<unsigned, unsigned> > intersects_;
};

void print_traxel_feature(const Traxel& trax, const std::string& feature) {
  std::cout << feature << '\n';
  const feature_array& fa = trax.features.find(feature)->second;
  for (feature_array::const_iterator it = fa.begin(); it != fa.end(); ++it) {
    std::cout << *it << ',';
  }
  std::cout << "\b \n";
}
  

void export_auxiliary() {
  class_<IntersectsDictionary, boost::noncopyable>("IntersectsDictionary")
      .def("getIntersects", &IntersectsDictionary::get_intersects_2d)
      .def("getIntersects", &IntersectsDictionary::get_intersects_3d)
      .def("addToTraxel", &IntersectsDictionary::add_to_traxel)
      ;

  def("printTraxelFeature", &print_traxel_feature);
}



