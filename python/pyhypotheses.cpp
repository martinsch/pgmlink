#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

#include <string>
#include <sstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/python.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/utility.hpp>

#include <lemon/core.h>

#include "../include/pgmlink/hypotheses.h"

using namespace pgmlink;
using namespace boost::python;

typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_m;


node_traxel_m& addNodeTraxelMap(HypothesesGraph* g) {
  g->add(node_traxel());
  return g->get(node_traxel());
}
node_traxel_m& getNodeTraxelMap(HypothesesGraph* g) {
  return g->get(node_traxel());
}

struct HypothesesGraph_pickle_suite : pickle_suite {
  static std::string getstate( const HypothesesGraph& g ) {
    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa & g;
    return ss.str();
  }

  static void setstate( HypothesesGraph& g, const std::string& state ) {
    std::stringstream ss(state);
    boost::archive::text_iarchive ia(ss);
    ia & g;
  }
};

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



void export_hypotheses() {
  class_<HypothesesGraph::Arc>("Arc");
  class_<HypothesesGraph::ArcIt>("ArcIt");
  class_<HypothesesGraph::Node>("Node");
  class_<HypothesesGraph::NodeIt>("NodeIt");
  class_<HypothesesGraph::InArcIt>("InArcIt");
  class_<HypothesesGraph::OutArcIt>("OutArcIt");


  IterableValueMap_ValueIterator<node_traxel_m>::wrap("NodeTraxelMap_ValueIt");

  class_<node_traxel_m, boost::noncopyable>("NodeTraxelMap", init<const HypothesesGraph&>(args("hypotheses_graph")))
    .def("__getitem__", &node_traxel_m::operator[],
	 return_internal_reference<>())
    .def("__setitem__", &node_traxel_m::set)
    .def("values", &IterableValueMap_ValueIterator<node_traxel_m>::values)
    ;


  // handle function overloading
  HypothesesGraph::Node (HypothesesGraph::*addnode1)(const int) 
    = &HypothesesGraph::add_node;
  HypothesesGraph::Node (HypothesesGraph::*addnode2)(const std::vector<int>) 
    = &HypothesesGraph::add_node;
  void (HypothesesGraph::*erase1)(const HypothesesGraph::Node) = &HypothesesGraph::erase;
  void (HypothesesGraph::*erase2)(const HypothesesGraph::Arc) = &HypothesesGraph::erase;
  bool (HypothesesGraph::*valid1)(const HypothesesGraph::Node) const = &HypothesesGraph::valid;
  bool (HypothesesGraph::*valid2)(const HypothesesGraph::Arc) const = &HypothesesGraph::valid;
  int (*id1)(HypothesesGraph::Node) = &HypothesesGraph::id;
  int (*id2)(HypothesesGraph::Arc) = &HypothesesGraph::id;
  HypothesesGraph::Node (HypothesesGraph::*baseNode1)(const HypothesesGraph::InArcIt&) const =
    &HypothesesGraph::baseNode;
  HypothesesGraph::Node (HypothesesGraph::*baseNode2)(const HypothesesGraph::OutArcIt&) const =
    &HypothesesGraph::baseNode;
  HypothesesGraph::Node (HypothesesGraph::*runningNode1)(const HypothesesGraph::InArcIt&) const =
    &HypothesesGraph::runningNode;
  HypothesesGraph::Node (HypothesesGraph::*runningNode2)(const HypothesesGraph::OutArcIt&) const =
    &HypothesesGraph::runningNode;

  class_<HypothesesGraph,shared_ptr<HypothesesGraph>, boost::noncopyable>("HypothesesGraph")
    .def("addNode", addnode1)
    .def("addNode", addnode2)
    //.def("timesteps", &HypothesesGraph::timesteps, 
    //   return_internal_reference<>())
    .def("earliest_timestep", &HypothesesGraph::earliest_timestep)
    .def("latest_timestep", &HypothesesGraph::latest_timestep)

    // lemon graph interface
    .def("addArc", &HypothesesGraph::addArc)
    .def("erase", erase1)
    .def("erase", erase2)
    .def("valid", valid1)
    .def("valid", valid2)
    .def("changeTarget", &HypothesesGraph::changeTarget)
    .def("changeSource", &HypothesesGraph::changeSource)
    .def("target", &HypothesesGraph::target)
    .def("source", &HypothesesGraph::source)
    .def("id", id1)
    .def("id", id2)
    .def("nodeFromId", &HypothesesGraph::nodeFromId)
    .def("arcFromId", &HypothesesGraph::arcFromId)
    .def("maxNodeId", &HypothesesGraph::maxNodeId)
    .def("maxArcId", &HypothesesGraph::maxArcId)
    .def("baseNode", baseNode1)
    .def("runningNode", runningNode1)
    .def("baseNode", baseNode2)
    .def("runningNode", runningNode2)
    .def("oppositeNode", &HypothesesGraph::oppositeNode)

    // extensions
    .def("addNodeTraxelMap", &addNodeTraxelMap,
	 return_internal_reference<>())
    .def("getNodeTraxelMap", &getNodeTraxelMap,
	 return_internal_reference<>())
    .def_pickle(HypothesesGraph_pickle_suite())
    ;

  //
  // lemon
  //
  int (*countNodes)(const HypothesesGraph&) = lemon::countNodes;
  def("countNodes", countNodes);

  int (*countArcs)(const HypothesesGraph&) = lemon::countArcs;
  def("countArcs", countArcs);
}
