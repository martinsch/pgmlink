#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

#include <vector>
#include <string>
#include <sstream>

#include "../include/pgmlink/traxels.h"
#include <vigra/multi_array.hxx>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/utility.hpp>



using namespace pgmlink;
using namespace boost::python;

namespace pgmlink {
  using namespace std;
  using namespace vigra;

  // extending Traxel
  void set_intmaxpos_locator(Traxel& t) {
    Locator* l = new IntmaxposLocator();
    t.set_locator(l); // takes ownership of pointer
  }

  void set_x_scale(Traxel& t, double s) {
    t.locator()->x_scale = s;
  }
  void set_y_scale(Traxel& t, double s) {
    t.locator()->y_scale = s;
  }
  void set_z_scale(Traxel& t, double s) {
    t.locator()->z_scale = s;
  }

  void add_feature_array(Traxel& t, string key, size_t size) {
    t.features[key] = feature_array(size, 0);
  }

  float get_feature_value(Traxel& t, string key, MultiArrayIndex i) {
    FeatureMap::const_iterator it = t.features.find(key);
    if(it == t.features.end()) {
      throw std::runtime_error("key not present in feature map");
    }
    if( !(static_cast<size_t>(i) < it->second.size())) {
      throw std::runtime_error("index out of range");
    }

    return it->second[i];
  }

  void set_feature_value(Traxel& t, string key, MultiArrayIndex i, float value) {
    FeatureMap::iterator it = t.features.find(key);
    if(it == t.features.end()) {
      throw std::runtime_error("key not present in feature map");
    }
    if( !(static_cast<size_t>(i) < it->second.size())) {
      throw std::runtime_error("index out of range");
    }
    it->second[i] = value;
  }

  // extending Traxels
  void add_traxel(Traxels& ts, const Traxel& t) {
    ts[t.Id] = t;
  }

  // extending TraxelStore
  struct TraxelStore_pickle_suite : pickle_suite {
    static std::string getstate( const TraxelStore& g ) {
      std::stringstream ss;
      boost::archive::text_oarchive oa(ss);
      oa & g;
      return ss.str();
    }
    
    static void setstate( TraxelStore& g, const std::string& state ) {
      std::stringstream ss(state);
      boost::archive::text_iarchive ia(ss);
      ia & g;
    }
  };

  void add_traxel_to_traxelstore(TraxelStore& ts, const Traxel& t) {
    add(ts, t);
  }

  void add_Traxels_to_traxelstore(TraxelStore& ts, const Traxels& traxels) {
    for(Traxels::const_iterator it = traxels.begin(); it!= traxels.end(); ++it){
      add(ts, it->second);
    }
  }

} /* namespace pgmlink */

void export_traxels() {
    class_< feature_array >("feature_array")
      .def(vector_indexing_suite< feature_array >() )
      .def("push_back", (void (feature_array::*)(feature_array::const_reference))&feature_array::push_back)
      .def("size", &feature_array::size)
      ;

    class_< ComLocator >("ComLocator")
      .def_readwrite("x_scale", &ComLocator::x_scale)
      .def_readwrite("y_scale", &ComLocator::y_scale)
      .def_readwrite("z_scale", &ComLocator::z_scale)
    ;

    class_< IntmaxposLocator >("IntmaxposLocator")
      .def_readwrite("x_scale", &ComLocator::x_scale)
      .def_readwrite("y_scale", &ComLocator::y_scale)
      .def_readwrite("z_scale", &ComLocator::z_scale)
    ;

    class_< std::map<std::string,feature_array> >("FeatureMap")
	.def(map_indexing_suite<std::map<std::string,feature_array> >())
    ;
    class_<Traxel>("Traxel")
	.def_readwrite("Id", &Traxel::Id)
	.def_readwrite("Timestep", &Traxel::Timestep)
        //.def("set_locator", &Traxel::set_locator, return_self<>())
        .def("set_x_scale", &set_x_scale)
        .def("set_y_scale", &set_y_scale)
        .def("set_z_scale", &set_z_scale)
        .def("set_intmaxpos_locator", &set_intmaxpos_locator, args("self"))
        .def("X", &Traxel::X)
        .def("Y", &Traxel::Y)
        .def("Z", &Traxel::Z)
	.def_readwrite("features", &Traxel::features)
        .def("add_feature_array", &add_feature_array, args("self","name", "size"), "Add a new feature array to the features map; initialize with zeros. If the name is already present, the old feature array will be replaced.")
	.def("get_feature_value", &get_feature_value, args("self", "name", "index"))
	.def("set_feature_value", &set_feature_value, args("self", "name", "index", "value"))
    ;

    class_<map<unsigned int, Traxel> >("Traxels")
	.def("add_traxel", &add_traxel)
	.def("__len__", &Traxels::size)
    ;

    class_< std::vector<double> >("VectorOfDouble")
      .def(vector_indexing_suite< std::vector<double> >() )
    ;

    class_< std::vector<int> >("VectorOfInt")
      .def(vector_indexing_suite< std::vector<int> >() )
    ;

    TraxelStoreByTimeid& (TraxelStore::*get_by_timeid)() = &TraxelStore::get<by_timeid>; 
    class_<TraxelStore>("TraxelStore")
      .def("add", &add_traxel_to_traxelstore)
      .def("add_from_Traxels", &add_Traxels_to_traxelstore)
      .def("bounding_box", &bounding_box)
      .def("get_by_timeid", get_by_timeid, return_internal_reference<>())
      .def("size", &TraxelStore::size)
      .def_pickle(TraxelStore_pickle_suite())
      ;
}
