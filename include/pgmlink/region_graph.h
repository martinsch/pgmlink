#ifndef REGION_GRAPH_H
#define REGION_GRAPH_H

// stl
#include <set>
#include <map>

// boost
#include <boost/shared_ptr.hpp>

// lemon
#include <lemon/list_graph.h> // lemon::ListGraph
#include <lemon/maps.h> // lemon::IterableValueMap


// pgmlink
#include <pgmlink/graph.h>


namespace pgmlink {

  class RegionGraph;

  typedef boost::shared_ptr<RegionGraph> RegionGraphPtr;

  typedef unsigned label_type;
  

  ////
  //// node_neighbors
  ////  
  struct node_neighbors {};
  template <typename Graph>
  struct property_map<node_neighbors, Graph> {
    typedef lemon::IterableValueMap<Graph, typename Graph::Node, std::set<typename Graph::Node> > type;
    static const std::string name;
  };
  
  template <typename Graph>
  const std::string property_map<node_neighbors, Graph>::name = "node_neighbors";


  ////
  //// arc_dissimilarity
  ////
  struct arc_dissimilarity {};
  template <typename Graph>
  struct property_map<arc_dissimilarity, Graph> {
    typedef lemon::IterableValueMap<Graph, typename Graph::Arc, double> type;
    static const std::string name;
  };

  template <typename Graph>
  const std::string property_map<arc_dissimilarity, Graph>::name = "arc_dissimilarity";


  ////
  //// node_label
  ////
  struct node_label {};
  template <typename Graph>
  struct property_map<node_label, Graph> {
    typedef lemon::IterableValueMap<Graph, typename Graph::Node, label_type> type;
    static const std::string name;
  };

  template <typename Graph>
  const std::string property_map<node_label, Graph>::name = "node_label";


  ////
  //// node_contains
  ////
  struct node_contains {};
  template <typename Graph>
  struct property_map<node_contains, Graph> {
    typedef lemon::IterableValueMap<Graph, typename Graph::Node, std::set<typename Graph::Node> > type;
    static const std::string name;
  };

  template <typename Graph>
  const std::string property_map<node_contains, Graph>::name = "node_contains";


  ////
  //// node_confligts
  struct node_conflicts {};
  template <typename Graph>
  struct property_map<node_conflicts, Graph> {
    typedef lemon::IterableValueMap<Graph, typename Graph::Node, std::set<typename Graph::Node> > type;
    static const std::string name;
  };

  template <typename Graph>
  const std::string property_map<node_conflicts, Graph>::name = "node_conflicts";


  ////
  //// RegionGraph
  ////
  class RegionGraph : public PropertyGraph<lemon::ListGraph> {
  public:
    typedef base_graph::Node Region;
    typedef base_graph::Arc Dissimilarity;
    typedef property_map<node_neighbors, RegionGraph> NeighborMap;
    typedef property_map<arc_dissimilarity, RegionGraph> DissimilarityMap;
    typedef property_map<node_label, RegionGraph> LabelMap;
    typedef property_map<node_contains, RegionGraph> ContainingMap;
    typedef property_map<node_conflicts, RegionGraph> ConflictMap;
  private:
    label_type maximum_label_;
  public:
    RegionGraph();
    int merge_regions(Region r1, Region r2);
    int merge_regions(int label1, int label2);
  };

}


#endif /* REGION_GRAPH_H */
