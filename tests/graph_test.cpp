#define BOOST_TEST_MODULE graph_test

#include <vector>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <lemon/core.h>
#include <lemon/concepts/digraph.h>
#include <lemon/list_graph.h>
#include <lemon/maps.h>

#include "graph.h"

using namespace Tracking;
using namespace std;
using namespace boost;


namespace Tracking {
// graph property for testing purposes
struct node_testprop {};
template <typename Graph>
struct property_map<node_testprop, Graph> {
  typedef lemon::IterableValueMap< Graph, typename Graph::Node, int > type;
  static const std::string name;
};
template <typename Graph>
const std::string property_map<node_testprop,Graph>::name = "node_testprop";
}

BOOST_AUTO_TEST_CASE( PropertyGraph_add_get ) {
  typedef PropertyGraph<lemon::ListDigraph> prop_graph;
  prop_graph g;
  prop_graph::base_graph::Node n = g.addNode();
  g.add(node_testprop());
  property_map<node_testprop,prop_graph::base_graph>::type& m 
    = g.get(node_testprop());
  m.set(n, 72);
  property_map<node_testprop,prop_graph::base_graph>::type& other_m 
    = g.get(node_testprop());
  BOOST_CHECK_EQUAL(other_m[n], 72);
}

BOOST_AUTO_TEST_CASE( PropertyGraph_addNode ) {
    PropertyGraph<lemon::ListDigraph> graph;
    graph.addNode();
}

// EOF

