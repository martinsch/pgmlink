#define BOOST_TEST_MODULE region_graph_test

// stdl
#include <set>
#include <iostream>

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
  std::set<RegionGraph::Region> comparison;
  
  pgmlink::RegionGraph g;
  
  RegionGraph::ContainingMap& containing_map = g.get(node_contains());
  RegionGraph::ConflictMap& conflict_map = g.get(node_conflicts());
  RegionGraph::NeighborMap& neighbor_map = g.get(node_neighbors());

  RegionGraph::Region r1 = g.add_region(1u);
  BOOST_CHECK_EQUAL(g.get_maximum_label(), 1u);
  RegionGraph::Region r2 = g.add_region(2u);
  BOOST_CHECK_EQUAL(g.get_maximum_label(), 2u);
  RegionGraph::Region r3 = g.add_region(3u);
  BOOST_CHECK_EQUAL(g.get_maximum_label(), 3u);

  neighbor_map.get_value(r1).insert(r2);
  neighbor_map.get_value(r1).insert(r3);

  neighbor_map.get_value(r2).insert(r3);
  neighbor_map.get_value(r2).insert(r1);

  neighbor_map.get_value(r3).insert(r1);
  neighbor_map.get_value(r3).insert(r2);
  
  
  RegionGraph::Region r4 = g.merge_regions(r1,r2);
  comparison.insert(r1);
  comparison.insert(r2);
  BOOST_CHECK(conflict_map[r4] == comparison);
  
  RegionGraph::Region r5 = g.merge_regions(2u, 3u);

  RegionGraph::LabelMap& label_map = g.get(node_label());

  BOOST_CHECK_EQUAL(label_map[r1], 1u);
  BOOST_CHECK_EQUAL(label_map[r2], 2u);
  BOOST_CHECK_EQUAL(label_map[r3], 3u);
  BOOST_CHECK_EQUAL(label_map[r4], 4u);
  BOOST_CHECK_EQUAL(label_map[r5], 5u);

  BOOST_CHECK_EQUAL(g.get_maximum_label(), 5u);

  BOOST_CHECK_EQUAL(containing_map[r1].size(), 0u);
  BOOST_CHECK_EQUAL(containing_map[r2].size(), 0u);
  BOOST_CHECK_EQUAL(containing_map[r3].size(), 0u);
  BOOST_CHECK_EQUAL(containing_map[r4].size(), 2u);
  BOOST_CHECK_EQUAL(containing_map[r5].size(), 2u);

  comparison.clear();
  comparison.insert(r1);
  comparison.insert(r2);
  BOOST_CHECK(containing_map[r4] == comparison);

  comparison.insert(r5);
  BOOST_CHECK(conflict_map[r4] == comparison);

  
  comparison.clear();
  comparison.insert(r2);
  comparison.insert(r3);
  BOOST_CHECK(containing_map[r5] == comparison);

  comparison.insert(r4);
  BOOST_CHECK(conflict_map[r5] == comparison);

  comparison.clear();
  comparison.insert(r1);
  comparison.insert(r2);
  comparison.insert(r4);
  BOOST_CHECK(neighbor_map[r3] == comparison);

  comparison.clear();
  comparison.insert(r2);
  comparison.insert(r3);
  comparison.insert(r5);
  BOOST_CHECK(neighbor_map[r1] == comparison);

  comparison.clear();
  comparison.insert(r4);
  comparison.insert(r5);
  BOOST_CHECK(conflict_map[r2] == comparison);

  BOOST_CHECK(g.merge_regions(10u, 11u) == lemon::INVALID);

  RegionGraph::Region r6 = g.merge_regions(r4, r5);
  comparison.clear();
  comparison.insert(r1);
  comparison.insert(r2);
  comparison.insert(r3);
  comparison.insert(r4);
  comparison.insert(r5);
  BOOST_CHECK(conflict_map[r6] == comparison);

  comparison.clear();
  comparison.insert(r4);
  comparison.insert(r5);
  BOOST_CHECK(containing_map[r6] == comparison);
  
  BOOST_CHECK_EQUAL(neighbor_map[r6].size(), 0u);
  


}
