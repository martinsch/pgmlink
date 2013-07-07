#define BOOST_TEST_MODULE region_graph_test


// boost
#include <boost/test/unit_test.hpp>

// pgmlink
#include <pgmlink/region_graph.h>

using namespace pgmlink;

BOOST_AUTO_TEST_CASE( RegionGraph_constructor ) {
  RegionGraph g;
  BOOST_CHECK(g.has_property(node_neighbors()));
  BOOST_CHECK(g.has_property(arc_dissimilarity()));
  BOOST_CHECK(g.has_property(node_label()));
  BOOST_CHECK(g.has_property(node_contains()));
  BOOST_CHECK(g.has_property(node_conflicts()));
  BOOST_CHECK(g.has_property(node_level()));
  BOOST_CHECK_EQUAL(g.get_maximum_label(), 0u);
}


BOOST_AUTO_TEST_CASE( RegionGraph_functions ) {
  pgmlink::RegionGraph g;
  RegionGraph::Region r1 = g.add_region(1u);
  BOOST_CHECK_EQUAL(g.get_maximum_label(), 1u);
  RegionGraph::Region r2 = g.add_region(2u);
  BOOST_CHECK_EQUAL(g.get_maximum_label(), 2u);
  RegionGraph::Region r3 = g.add_region(3u);
  BOOST_CHECK_EQUAL(g.get_maximum_label(), 3u);
  g.merge_regions(r1,r2);
  g.merge_regions(r2,r3);
}
