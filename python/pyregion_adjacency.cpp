#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

#include <string>
#include <sstream>

#include <boost/python.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/utility.hpp>

#include <lemon/core.h>

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

#include "../include/pgmlink/region_adjacency.h"

using namespace pgmlink;
using namespace boost::python;


typedef property_map<edge_weight, typename RegionAdjacencyGraph::base_graph>::type edge_weight_m;


template<unsigned int N, class T>
void buildGraphPython(RegionAdjacencyGraph& g, const vigra::NumpyArray<N,T> segmentImage) {
	g.buildRegionAdjacencyGraph<N,T>(segmentImage);
}



inline object pass_through(object const& o) { return o; }
template < typename ITERABLE_VALUE_MAP >
struct IterableValueMap_ValueIterator {
  typedef ITERABLE_VALUE_MAP map_type;
  IterableValueMap_ValueIterator( const map_type& m ) {
    it_ = m.beginValue();
    end_ = m.endValue();
  }
  typename map_type::Value next() {
    if( it_ == end_) {
        PyErr_SetString(PyExc_StopIteration, "No more data.");
        boost::python::throw_error_already_set();
    }
    const typename map_type::Value& result = *it_;
    ++it_;
    return result;
  }

  static IterableValueMap_ValueIterator<map_type>
  values ( const map_type& m ) {
    return IterableValueMap_ValueIterator<map_type>( m );
  }

  static void
  wrap( const char* python_name) {
    class_<IterableValueMap_ValueIterator<map_type> >( python_name, init<const map_type&>(args("iterable_value_map")) )
      .def("next", &IterableValueMap_ValueIterator::next)
      .def("__iter__", &pass_through)
      ;
  }

  typename map_type::ValueIt it_;
  typename map_type::ValueIt end_;
};


void export_region_adjacency() {
  class_<Region>("Region")
	.def_readwrite("id", &Region::id)
	.def_readwrite("contains_labels", &Region::contains_labels)
	;

  class_<typename RegionAdjacencyGraph::Edge>("Edge");
  class_<typename RegionAdjacencyGraph::EdgeIt>("EdgeIt");
  class_<typename RegionAdjacencyGraph::Node>("Node");
  class_<typename RegionAdjacencyGraph::NodeIt>("NodeIt");
  class_<typename RegionAdjacencyGraph::IncEdgeIt>("IncEdgeIt");

  IterableValueMap_ValueIterator<edge_weight_m>::wrap("EdgeWeightMap_ValueIt");
  class_<edge_weight_m, boost::noncopyable>("EdgeWeightMap", init<const RegionAdjacencyGraph&>(args("region_adjacency_graph")))
      .def("__getitem__", &edge_weight_m::operator[], return_value_policy<copy_const_reference>() )
      .def("__setitem__", &edge_weight_m::set)
      .def("values", &IterableValueMap_ValueIterator<edge_weight_m>::values)
      ;

  class_<std::vector<std::vector<int> > >("VectorOfVectorOfInt")
        .def(vector_indexing_suite<std::vector<std::vector<int> > >())
      ;

  class_<std::vector<Region> >("RegionVector")
		 .def(vector_indexing_suite<std::vector<Region> >())
		 ;

  class_<std::vector<int> >("IntVector")
		  .def(vector_indexing_suite<std::vector<int> >())
		  ;

  class_<std::map<int, std::vector<std::vector<int> > > >("ConflictSetsMap")
		  .def(map_indexing_suite<std::map<int, std::vector<std::vector<int> > > >())
		  ;


  // handle function overloading
  void (RegionAdjacencyGraph::*erase1)(const RegionAdjacencyGraph::Node) = &RegionAdjacencyGraph::erase;
  void (RegionAdjacencyGraph::*erase2)(const RegionAdjacencyGraph::Edge) = &RegionAdjacencyGraph::erase;
  bool (RegionAdjacencyGraph::*valid1)(const RegionAdjacencyGraph::Node) const = &RegionAdjacencyGraph::valid;
  bool (RegionAdjacencyGraph::*valid2)(const RegionAdjacencyGraph::Edge) const = &RegionAdjacencyGraph::valid;
  int (*id1)(RegionAdjacencyGraph::Node) = &RegionAdjacencyGraph::id;
  int (*id2)(RegionAdjacencyGraph::Edge) = &RegionAdjacencyGraph::id;

  class_<RegionAdjacencyGraph, boost::noncopyable>("RegionAdjacencyGraph")
    .def("addNode", &RegionAdjacencyGraph::addNode)

    // lemon graph interface
    .def("addEdge", &RegionAdjacencyGraph::addEdge)
    .def("erase", erase1)
    .def("erase", erase2)
    .def("valid", valid1)
    .def("valid", valid2)
    .def("target", &RegionAdjacencyGraph::target)
    .def("source", &RegionAdjacencyGraph::source)
    .def("id", id1)
    .def("id", id2)
    .def("nodeFromId", &RegionAdjacencyGraph::nodeFromId)
    .def("edgeFromId", &RegionAdjacencyGraph::edgeFromId)
    .def("maxNodeId", &RegionAdjacencyGraph::maxNodeId)
    .def("maxEdgeId", &RegionAdjacencyGraph::maxEdgeId)

    .def("buildGraph", vigra::registerConverters(&buildGraphPython<2u,double>))
    .def("buildGraph", vigra::registerConverters(&buildGraphPython<3u,double>))
    .def("buildGraph", vigra::registerConverters(&buildGraphPython<2u,unsigned int>))
    .def("buildGraph", vigra::registerConverters(&buildGraphPython<3u,unsigned int>))

    .def("getEdgeWeightMap", &RegionAdjacencyGraph::get_edge_weight_map, return_internal_reference<>())
    .def("getLabelToNodeMap", &RegionAdjacencyGraph::get_label_to_node_map, return_internal_reference<>())
    .def("mergeNodesThreshold", &RegionAdjacencyGraph::merge_nodes_threshold)

    .def("getLabelsVector", &RegionAdjacencyGraph::get_labels_vector)
    .def("getConnectedComponentIds", &RegionAdjacencyGraph::get_connected_component_ids)
    .def("getRegions", &RegionAdjacencyGraph::get_regions)
    .def("getConflictSets", &RegionAdjacencyGraph::get_conflict_sets)

    ;

  //
  // lemon
  //
  int (*countNodes)(const RegionAdjacencyGraph&) = lemon::countNodes;
  def("countNodes", countNodes);

  int (*countEdges)(const RegionAdjacencyGraph&) = lemon::countEdges;
  def("countEdges", countEdges);
}
