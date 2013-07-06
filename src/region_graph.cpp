// stl
#include <set> // std::set
#include <algorithm> // std::copy
#include <functional> //std::not_equal_to
#include <iterator> // std::insert_iterator

// boost
#include <boost/range/adaptors.hpp> // boost::adaptors::filter
#include <boost/range/algorithm.hpp> // boost::copy
#include <boost/mpl/copy_if.hpp> // boost::copy_if

// lemon
#include <lemon/list_graph.h> // lemon::ListGraph
#include <lemon/maps.h> // lemon::IterableValueMap

//pgmlink
#include <pgmlink/region_graph.h>


namespace ba = boost::adaptors;


namespace pgmlink {
  ////
  //// RegionGraph
  ////
  RegionGraph::RegionGraph() :
    maximum_label_(0u) {
    add(node_neighbors()).add(arc_dissimilarity()).add(node_label()).add(node_contains()).add(node_conflicts());
  }


  void RegionGraph::union_neighbors(const Region& r1,
                                    const Region& r2,
                                    const Region& new_region) {
    NeighborMap& neighbor_map = get(node_neighbors());
    const std::set<Region>& neighbors1 = neighbor_map[r1];
    const std::set<Region>& neighbors2 = neighbor_map[r2];

    std::set<Region>& neighbors_new = neighbor_map.get_value(new_region);

    boost::copy(neighbors1 | ba::filtered(std::bind1st(std::not_equal_to<Region>(), r2)),
                std::inserter(neighbors_new, neighbors_new.begin()));
    boost::copy(ba::filter(neighbors2,
                           std::bind1st(std::not_equal_to<Region>(), r1)),
                           std::inserter(neighbors_new, neighbors_new.begin()));
  }

  
  int RegionGraph::merge_regions(Region r1, Region r2) {
    if (!valid(r1)) {
      return 1;
    }
    if (!valid(r2)) {
      return 2;
    }
    Region new_region = addNode();
    ++maximum_label_;

    
    // DissimilarityMap& dissimilarity_map = get(arc_dissimilarity());
    LabelMap& labels_map = get(node_label());
    // ContainingMap& contains_map = get(node_contains());
    // ConflictMap& conflicts_map = get(node_conflicts());

    union_neighbors(r1, r2, new_region);
    labels_map.set(new_region, maximum_label_);
    

    
    return 0;
  }

  int RegionGraph::merge_regions(int, int) {
    

    return 0;
  }

}
