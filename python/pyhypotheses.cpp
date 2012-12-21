#include <boost/python.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/utility.hpp>

#include <lemon/core.h>

#include "../include/pgmlink/hypotheses.h"

using namespace pgmlink;
using namespace boost::python;

typedef typename property_map<node_traxel, typename HypothesesGraph::base_graph>::type node_traxel_m;


node_traxel_m& addNodeTraxelMap(HypothesesGraph* g) {
  g->add(node_traxel());
  return g->get(node_traxel());
}
node_traxel_m& getNodeTraxelMap(HypothesesGraph* g) {
  return g->get(node_traxel());
}


inline object pass_through(object const& o) { return o; }

struct lemon_iterator_wrapper {
  lemon_iterator_wrapper( node_traxel_m::ValueIt endValue )
    : endValue_(endValue) {}

  static node_traxel_m::Value
  next( node_traxel_m::ValueIt& it ){
    if( it == endValue_) {
        PyErr_SetString(PyExc_StopIteration, "No more data.");
        boost::python::throw_error_already_set();
    }
    const node_traxel_m::Value& result = *it;
    ++it;
    return result;
  }

  static void
  wrap( const char* python_name) {
    class_<typename node_traxel_m::ValueIt>( python_name )
      .def("next", next)
      .def("__iter__", pass_through) 
      ;
  }

  node_traxel_m::ValueIt endValue_;
};

void export_hypotheses() {
  class_<typename HypothesesGraph::Arc>("Arc");
  class_<typename HypothesesGraph::ArcIt>("ArcIt");
  class_<typename HypothesesGraph::Node>("Node");
  class_<typename HypothesesGraph::NodeIt>("NodeIt");
  class_<typename HypothesesGraph::InArcIt>("InArcIt");
  class_<typename HypothesesGraph::OutArcIt>("OutArcIt");

  // NodeTraxelMap
  // class_<typename node_traxel_m::ValueIt>("NodeTraxelMap_ValueIt")
  //   .def("next", &node_traxel_m::ValueIt::operator++)
  //   ;

  class_<node_traxel_m, boost::noncopyable>("NodeTraxelMap", init<const HypothesesGraph&>(args("hypotheses_graph")))
    .def("__getitem__", &node_traxel_m::operator[],
	 return_internal_reference<>())
    .def("__setitem__", &node_traxel_m::set)
    ;


  // handle function overloading
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

  class_<HypothesesGraph, boost::noncopyable>("HypothesesGraph")
    .def("addNode", &HypothesesGraph::add_node)
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
    ;

  //
  // lemon
  //
  int (*countNodes)(const HypothesesGraph&) = lemon::countNodes;
  def("countNodes", countNodes);

  int (*countArcs)(const HypothesesGraph&) = lemon::countArcs;
  def("countArcs", countArcs);
}
