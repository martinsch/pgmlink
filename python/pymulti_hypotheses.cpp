#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

// stl

// boost
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/python.hpp>
#include <boost/python/return_internal_reference.hpp>

// vigra
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

// pgmlink
#include "pgmlink/multi_hypotheses_graph.h"


using namespace pgmlink;
using namespace boost::python;

MultiHypothesesGraph::ContainedRegionsMap& addContainedRegionsMap(MultiHypothesesGraph* g) {
  g->add(node_regions_in_component());
  return g->get(node_regions_in_component());
}

MultiHypothesesGraph::ContainedRegionsMap& getContainedRegionsMap(MultiHypothesesGraph* g) {
  return g->get(node_regions_in_component());
}

inline object pass_through(object const& o) { return o; }

template <typename ITERABLE_VALUE_MAP>
struct IterableValueMap_ValueIterator {
  typedef ITERABLE_VALUE_MAP map_type;
  IterableValueMap_ValueIterator( const map_type& m ) {
    it_ = m.beginValue();
    end_ = m.endValue();
  }
  typename map_type::Value next() {
    if (it_ == end_) {
      PyErr_SetString(PyExc_StopIteration, "No more data.");
      boost::python::throw_error_already_set();
    }
    const typename map_type::Value& result = *it_;
    ++it_;
    return result;
  }

  static IterableValueMap_ValueIterator<map_type>
  values( const map_type& m ) {
    return IterableValueMap_ValueIterator<map_type>( m );
  }

  static void
  wrap( const char* python_name) {
    class_<IterableValueMap_ValueIterator<map_type> >(python_name, init<const map_type&>(args("iterable_value_map")) )
        .def("next", &IterableValueMap_ValueIterator::next)
        .def("__iter__", &pass_through)
        ;
  }

  typename map_type::ValueIt it_;
  typename map_type::ValueIt end_;
};



class PyMultiHypothesesTraxelStoreBuilder {
 public:
  PyMultiHypothesesTraxelStoreBuilder(MultiHypothesesTraxelStore& ts) : ts_(ts) {
    
  }
  template<int N, typename LABEL_TYPE>
  void build(vigra::NumpyArray<N+1, LABEL_TYPE> arr,
             vigra::NumpyArray<N, LABEL_TYPE> components,
             unsigned object_layer,
             int timestep,
             LABEL_TYPE object_label
             ) {
    MultiHypothesesTraxelStoreBuilder builder;
    builder.build<N, LABEL_TYPE>(ts_, arr, components, object_layer, timestep, object_label);
  }
 private:
  MultiHypothesesTraxelStore& ts_;
};


void export_multi_hypotheses() {
  class_<typename MultiHypothesesGraph::Arc>("Arc");
  class_<typename MultiHypothesesGraph::ArcIt>("ArcIt");
  class_<typename MultiHypothesesGraph::Node>("Node");
  class_<typename MultiHypothesesGraph::NodeIt>("NodeIt");
  class_<typename MultiHypothesesGraph::InArcIt>("InArcIt");
  class_<typename MultiHypothesesGraph::OutArcIt>("OutArcIt");

  IterableValueMap_ValueIterator<MultiHypothesesGraph::ContainedRegionsMap>::wrap("ContainedRegionsMap_ValueIt");
  // IterableValueMap_ValueIterator<node_traxel>::wrap("NodeTraxelMap_ValueIt");

  class_<MultiHypothesesGraph::ContainedRegionsMap, boost::noncopyable>("ContainedRegionsMap",
                                                                        init<const MultiHypothesesGraph&>(args("multi_hypotheses_graph"))
                                                                        )
      .def("__getitem__",
           &MultiHypothesesGraph::ContainedRegionsMap::operator[],
           return_internal_reference<>())
      .def("get_value",
           &MultiHypothesesGraph::ContainedRegionsMap::get_value,
           return_internal_reference<>())
      .def("values", &IterableValueMap_ValueIterator<MultiHypothesesGraph::ContainedRegionsMap>::values)
      ;

  // handle function overloading
  MultiHypothesesGraph::Node (HypothesesGraph::*addnode1)(const int)
      = &MultiHypothesesGraph::add_node;
  MultiHypothesesGraph::Node (HypothesesGraph::*addnode2)(const std::vector<int>)
      = &MultiHypothesesGraph::add_node;

  class_<MultiHypothesesGraph, boost::noncopyable>("MultiHypothesesGraph")
      .def("addNode", addnode1)
      .def("addNode", addnode2)
      ;


  ////
  //// MultiHypothesesTraxelStore
  ////
  class_<MultiHypothesesTraxelStore>("MultiHypothesesTraxelStore")
      .def("add", &MultiHypothesesTraxelStore::add,
           return_internal_reference<>())
      .def("startComponent", &MultiHypothesesTraxelStore::start_component,
           return_internal_reference<>())
      .def("make_string", &MultiHypothesesTraxelStore::print)
      .def("__str__", &MultiHypothesesTraxelStore::print)
      ;


  ////
  //// class MultiHypothesesTraxelStoreBuilder
  ////
  class_<PyMultiHypothesesTraxelStoreBuilder, boost::noncopyable>("MultiHypothesesTraxelStoreBuilder",
                                                                  init<MultiHypothesesTraxelStore&>()[with_custodian_and_ward<1, 2>()]
                                                                  )
      .def("build", vigra::registerConverters(&PyMultiHypothesesTraxelStoreBuilder::build<2, unsigned>), return_internal_reference<>())
      .def("build", vigra::registerConverters(&PyMultiHypothesesTraxelStoreBuilder::build<3, unsigned>), return_internal_reference<>())
      .def("build", vigra::registerConverters(&PyMultiHypothesesTraxelStoreBuilder::build<2, unsigned long>), return_internal_reference<>())
      .def("build", vigra::registerConverters(&PyMultiHypothesesTraxelStoreBuilder::build<3, unsigned long>), return_internal_reference<>())
      ;


  
}


