#define BOOST_TEST_MODULE region_graph_test


// boost
#include <boost/test/unit_test.hpp>

// pgmlink
#include <pgmlink/region_graph.h>


BOOST_AUTO_TEST_CASE( RegionGraph_constructor ) {
  pgmlink::RegionGraph g;
  BOOST_CHECK(g.has_property(pgmlink::node_neighbors()));
  BOOST_CHECK(g.has_property(pgmlink::arc_dissimilarity()));
  BOOST_CHECK(g.has_property(pgmlink::node_label()));
  BOOST_CHECK(g.has_property(pgmlink::node_contains()));
  BOOST_CHECK(g.has_property(pgmlink::node_conflicts()));
  BOOST_CHECK(g.has_property(pgmlink::node_level()));
}
