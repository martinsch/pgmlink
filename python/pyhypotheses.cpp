#include <boost/python.hpp>
#include <boost/utility.hpp>

#include "../include/pgmlink/hypotheses.h"

using namespace pgmlink;
using namespace boost::python;

void export_hypotheses() {
      class_<HypothesesGraph, boost::noncopyable>("HypothesesGraph")
      //.def("addNode", &HypothesesGraph::add_node)
      //.def("timesteps", &HypothesesGraph::timesteps)
      //.def("earliest_timestep", &HypothesesGraph::earliest_timestep)
      //.def("latest_timestep", &HypothesesGraph::latest_timestep)
    ;
}
