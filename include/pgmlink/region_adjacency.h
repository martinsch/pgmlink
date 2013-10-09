#ifndef REGION_ADJACENCY_H
#define REGION_ADJACENCY_H

#include <vector>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <set>
#include <boost/serialization/set.hpp>
#include <boost/shared_ptr.hpp>
#include <lemon/list_graph.h>
#include <lemon/maps.h>
#include <vigra/multi_array.hxx>
#include <vigra/multi_iterator.hxx>
#include <vigra/tinyvector.hxx>

#include "pgmlink/graph.h"
#include "pgmlink/log.h"
#include "pgmlink/pgmlink_export.h"

namespace pgmlink {
  ////
  //// RegionAdjacencyGraph
  ////

  // Properties of a RegionAdjacencyGraph

    // node_features
  struct node_features {};
  template <typename Graph>
    struct property_map<node_features, Graph> {
    typedef lemon::IterableValueMap< Graph, typename Graph::Node, std::map<std::string,double> > type;
    static const std::string name;
  };
  template <typename Graph>
    const std::string property_map<node_features,Graph>::name = "node_features";

  struct node_labels {};
    template <typename Graph>
      struct property_map<node_labels, Graph> {
      typedef lemon::IterableValueMap< Graph, typename Graph::Node, std::vector<int> > type;
      static const std::string name;
    };
    template <typename Graph>
      const std::string property_map<node_labels,Graph>::name = "node_labels";

	// edge_weight
	struct edge_weight {};
	template <typename Graph>
	  struct property_map<edge_weight, Graph> {
	  typedef lemon::IterableValueMap< Graph, typename Graph::Edge, double > type;
	  static const std::string name;
	};
	template <typename Graph>
	  const std::string property_map<edge_weight,Graph>::name = "edge_weight";


  class RegionAdjacencyGraph : public PropertyGraph<lemon::ListGraph> {
  public:
    RegionAdjacencyGraph() {
      // Properties attached to every RegionAdjacencyGraph
      add(node_features());
      add(edge_weight());
      add(node_labels());
    };

    typedef property_map<edge_weight, typename RegionAdjacencyGraph::base_graph>::type edge_weight_m;

    // use this instead of calling the parent graph directly
    RegionAdjacencyGraph::Node add_node( int label );
    RegionAdjacencyGraph::Node add_node( int label, std::map<std::string,double> features );
    RegionAdjacencyGraph::Edge add_edge( RegionAdjacencyGraph::Node n1, RegionAdjacencyGraph::Node n2 );

    double compute_edge_weight( RegionAdjacencyGraph::Node n1, RegionAdjacencyGraph::Node n2 );
    double get_edge_weight( RegionAdjacencyGraph::Edge e );
    edge_weight_m& get_edge_weight_map();
    std::map<int, RegionAdjacencyGraph::Node>& get_label_to_node_map();

    template<unsigned int N, typename T>
    void buildRegionAdjacencyGraph(const vigra::MultiArrayView<N, T> segmentImage, int background_value = 0);

	void merge_nodes_threshold(double threshold, bool greedy=false);

	std::vector<std::vector<int> > get_labels_vector();

  private:
    void incrementPerimeter(RegionAdjacencyGraph::Node n);
    void incrementIntersection(RegionAdjacencyGraph::Node n1, RegionAdjacencyGraph::Node n2);

    std::map<int, RegionAdjacencyGraph::Node> label_to_node_;
    std::map<RegionAdjacencyGraph::Node, int> node_perimeters_;
    std::map<RegionAdjacencyGraph::Node, std::map<RegionAdjacencyGraph::Node, int> > nodes_intersections_;

  };


/* IMPLEMENTATIONS */

template<unsigned int N, typename T>
void RegionAdjacencyGraph::buildRegionAdjacencyGraph(const vigra::MultiArrayView<N, T> segmentImage, int background_value) {
	LOG(logDEBUG) << "RegionAdjacencyGraph::getRegionAdjacencyGraph entered";
	std::map<int, std::set<int> > neighbors;

	LOG(logDEBUG2) << "looping over segment image";
	for (typename vigra::MultiArrayView<N, T>::const_iterator it = segmentImage.begin(); it != segmentImage.end(); ++it) {
		// look at neighboring pixels (right and bottom pixel) and check whether they are different
		int label = (int) *it;
		// if it is the background label, skip
		if (label == background_value) {
			continue;
		}
		// if the label is not in the graph yet, add a node
		RegionAdjacencyGraph::Node current_node = add_node(label);

		// check whether its neighboring pixels (bottom and right) are different
		vigra::TinyVector<long, N> coordinates = it.get<0>();
		for(unsigned int d = 0; d < N; ++d) {
			coordinates[d] += 1;
			if (segmentImage.isInside(coordinates)) {
				if ((int) segmentImage[coordinates] == background_value) {
					incrementPerimeter(current_node);
				} else if (segmentImage[coordinates] != *it) {
					incrementPerimeter(current_node);
					int neighbor_label = segmentImage[coordinates];
					LOG(logDEBUG4) << "found neighbor " << neighbor_label;
					RegionAdjacencyGraph::Node neighbor_node = add_node(neighbor_label);
					incrementPerimeter(neighbor_node);
					incrementIntersection(current_node, neighbor_node);
					neighbors[label].insert( neighbor_label );
					neighbors[neighbor_label].insert( label );
				}
			}
		}
	}

	LOG(logDEBUG2) << "looping over neighbors";
	for (std::map<int, std::set<int> >::const_iterator it_map = neighbors.begin(); it_map != neighbors.end(); ++it_map) {
		std::set<int> neighbors_set = neighbors[it_map->first];
		for (std::set<int>::const_iterator it = neighbors_set.begin(); it != neighbors_set.end(); ++it ) {
			RegionAdjacencyGraph::Node n1 = label_to_node_.find(*it)->second;
			RegionAdjacencyGraph::Node n2 = label_to_node_.find(it_map->first)->second;
			LOG(logDEBUG4) << "adding edge between " << *it << " and " << it_map->first ;
			assert(n1 != n2);
			add_edge(n1, n2);
		}
	}
}

} // namespace



#endif /* REGION_ADJACENCY_H */
