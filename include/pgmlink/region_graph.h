#ifndef REGION_GRAPH_H
#define REGION_GRAPH_H

// boost
#include <boost/shared_ptr.hpp>

// pgmlink
#include <pgmlink/graph.h>


namespace pgmlink {
  
  class RegionGraph;

  typedef boost::shared_ptr<RegionGraph> RegionGraphPtr;

  ////
  //// RegionGraph
  ////
  class RegionGraph : public PropertyGraph<lemon::ListGraph> {
  public:
    typedef base_graph::Node Region;
    typedef base_graph::Arc Dissimilarity;
  private:
  public:
    int merge_regions(Region r1, Region r2);
  };

}


#endif /* REGION_GRAPH_H */
