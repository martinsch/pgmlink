#ifndef HANSLOVSKY_H
#define HANSLOVSKY_H

// stl headers
#include <vector>
#include <stdexcept>


// external headers
#include <lemon/maps.h>


// pgmlink headers
#include "pgmlink/hypotheses.h"
#include "pgmlink/event.h"
#include "pgmlink/traxels.h"

namespace pgmlink {
  ////
  //// MergerResolver
  ////
  class MergerResolver {
  private:
    HypothesesGraph* g_;
    
    MergerResolver() {}
    // template <typename ArcIterator>
    void collect_arcs(ArcIterator,
		      std::vector<HypothesesGraph::base_graph::Arc>&);
    
    void add_arcs_for_replacement_node(HypothesesGraph::Node,
				       Traxel,
				       std::vector<HypotheseGraph::base_graph::Arc>,
				       std::vector<HypotheseGraph::base_graph::Arc>);
    
    void deactivate_arcs(std::vector<HypotheseGraph::base_graph::Arc>);
    
    void refine_node(HypothesesGraph::Node,
		     std::size_t;)
      
      public:
      MergerResolver(HypothesesGraph* g) : g_(g) {}
    HypothesesGraph* resolve_mergers();
    
  };
}


#endif /* HANSLOVSKY_H */
