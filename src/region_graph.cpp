// stl
#include <set>
#include <algorithm>

// lemon
#include <lemon/list_graph.h> // lemon::ListGraph
#include <lemon/maps.h> // lemon::IterableValueMap

//pgmlink
#include <pgmlink/region_graph.h>



namespace pgmlink {
  ////
  //// RegionGraph
  ////
  RegionGraph::RegionGraph() :
    maximum_label_(0u) {
    add(node_neighbors()).add(arc_dissimilarity()).add(node_label()).add(node_contains()).add(node_conflicts());
  }

  
  int RegionGraph::merge_regions(Region r1, Region r2) {
    if (!valid(r1)) {
      return 1;
    }
    if (!valid(r2)) {
      return 2;
    }
    /*Region new_region = addNode();
    ++maximum_label_;

    NeighborMap& neighbors = get(node_neighbors());
    DissimilarityMap& dissimilarity = get(node_dissimilarity());
    LabelMap& labels = get(node_labels());
    ContainingMap& contains = get(node_contains());
    ConflictMap& conflicts = get(node_conflicts());

    neighbors.set(new_region, create_union(r1, r2));
    labels.set(new_region, max_label);*/
    

    
    return 0;
  }

  int RegionGraph::merge_regions(int, int) {
    

    return 0;
  }

}
