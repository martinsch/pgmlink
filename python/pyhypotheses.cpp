#include <boost/python.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/utility.hpp>

#include <lemon/core.h>

#include "../include/pgmlink/hypotheses.h"

using namespace pgmlink;
using namespace boost::python;

void export_hypotheses() {
  class_<typename HypothesesGraph::Arc>("Arc");
  class_<typename HypothesesGraph::ArcIt>("ArcIt");
  class_<typename HypothesesGraph::Node>("Node");
  class_<typename HypothesesGraph::NodeIt>("NodeIt");
  class_<typename HypothesesGraph::InArcIt>("InArcIt");
  class_<typename HypothesesGraph::OutArcIt>("OutArcIt");

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
    ;

  //
  // lemon
  //
  int (*countNodes)(const HypothesesGraph&) = lemon::countNodes;
  def("countNodes", countNodes);

  int (*countArcs)(const HypothesesGraph&) = lemon::countArcs;
  def("countArcs", countArcs);
}
