#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <algorithm>

#include "pgmlink/region_adjacency.h"
#include "pgmlink/log.h"

using namespace std;

namespace pgmlink {
  ////
  //// class RegionAdjacencyGraph
  ////
RegionAdjacencyGraph::Node RegionAdjacencyGraph::add_node( int label ) {
	std::map<int, RegionAdjacencyGraph::Node>::const_iterator it = label_to_node_.find(label);
	if (it != label_to_node_.end()) {
		LOG(logDEBUG4) << "add_node: node " << label << " already exists, not added ";
		return it->second;
	}
	LOG(logDEBUG4) << "add_node: add " << label << " to graph";
	RegionAdjacencyGraph::Node node = addNode();
	label_to_node_[label] = node;

	property_map<node_labels, RegionAdjacencyGraph::base_graph>::type& node_labels_map = get(node_labels());
	std::vector<int> label_list = node_labels_map[node];
	label_list.push_back(label);
	node_labels_map.set(node,label_list);

	return node;
}

RegionAdjacencyGraph::Node RegionAdjacencyGraph::add_node( int label, std::map<std::string,double> features ) {
//	typedef property_map<node_features, RegionAdjacencyGraph::base_graph>::type node_features_map;
//	lemon::IterableValueMap< Graph, typename Graph::Node, std::map<std::string,double> >
	property_map<node_features, RegionAdjacencyGraph::base_graph>::type& node_feature_map = get(node_features());
	RegionAdjacencyGraph::Node node = add_node(label);

	node_feature_map.set(node, features);
	return node;
}


RegionAdjacencyGraph::Edge RegionAdjacencyGraph::add_edge( RegionAdjacencyGraph::Node n1, RegionAdjacencyGraph::Node n2 ) {
	property_map<edge_weight, RegionAdjacencyGraph::base_graph>::type& edge_weight_map = get(edge_weight());
	for (IncEdgeIt e(*this, n1); e != lemon::INVALID; ++e) {
		if (this->v(e) == n2 || this->u(e) == n2) {
			LOG(logDEBUG4) << "edge already exists";
			return e;
		}
	}
	RegionAdjacencyGraph::Edge edge = addEdge(n1, n2);
	LOG(logDEBUG4) << "edge added";
	edge_weight_map.set(edge, compute_edge_weight(n1, n2));

	return edge;
}

double RegionAdjacencyGraph::compute_edge_weight( RegionAdjacencyGraph::Node n1, RegionAdjacencyGraph::Node n2 ) {
	if (this->id(n2) < this->id(n1)) {
		return compute_edge_weight(n2, n1);
	}
	assert(this->id(n1) < this->id(n2));

	std::map<RegionAdjacencyGraph::Node, std::map<RegionAdjacencyGraph::Node, int> >::iterator it = nodes_intersections_.find(n1);
	if (it == nodes_intersections_.end() ) {
		throw runtime_error("node1 not contained in intersection_map");
	}
	std::map<RegionAdjacencyGraph::Node, int>& inner_map = it->second;
	std::map<RegionAdjacencyGraph::Node, int>::iterator inner_it = inner_map.find(n2);
	if (inner_it == inner_map.end() ) {
		throw runtime_error("node2 not contained in intersection_map");
	}

	double intersect = inner_it->second;

	std::map<RegionAdjacencyGraph::Node, int>::iterator it_p = node_perimeters_.find(n1);
	if( it_p == node_perimeters_.end() ) {
		throw runtime_error("node1 not contained in perimeter_map");
	}
	int perimeter1 = it_p->second;

	it_p = node_perimeters_.find(n2);
	if( it_p == node_perimeters_.end() ) {
		throw runtime_error("node2 not contained in perimeter_map");
	}
	int perimeter2 = it_p->second;

	double weight = intersect/min(perimeter1,perimeter2);
	LOG(logDEBUG3) << "intersect = " << intersect << ", perimeter1 = " << perimeter1 << ", perimeter2 = " << perimeter2 << ", edge weight = " << weight;

	return weight;
}

double RegionAdjacencyGraph::get_edge_weight( RegionAdjacencyGraph::Edge e ) {
	property_map<edge_weight, RegionAdjacencyGraph::base_graph>::type& edge_weight_map = get(edge_weight());
	return edge_weight_map[e];
}

void RegionAdjacencyGraph::incrementPerimeter(RegionAdjacencyGraph::Node n) {
	std::map<RegionAdjacencyGraph::Node, int>::iterator it = node_perimeters_.find(n);
	if( it == node_perimeters_.end() ) {
		node_perimeters_[n] = 1;
		return;
	}
	it->second += 1;
}

void RegionAdjacencyGraph::incrementIntersection(RegionAdjacencyGraph::Node n1, RegionAdjacencyGraph::Node n2) {
	if (this->id(n2) < this->id(n1) ) { // only store node_intersections from smaller to bigger id to prevent redundancy
		incrementIntersection(n2, n1);
		return;
	}

	assert(this->id(n1) < this->id(n2));
	std::map<RegionAdjacencyGraph::Node, std::map<RegionAdjacencyGraph::Node, int> >::iterator it = nodes_intersections_.find(n1);
	if (it == nodes_intersections_.end() ) {
		std::map<RegionAdjacencyGraph::Node, int> inner_map;
		inner_map[n2] = 1;
		nodes_intersections_[n1] = inner_map;
		return;
	} else {
		std::map<RegionAdjacencyGraph::Node, int>& inner_map = it->second;
		std::map<RegionAdjacencyGraph::Node, int>::iterator inner_it = inner_map.find(n2);
		if (inner_it == inner_map.end() ) {
			inner_map[n2] = 1;
		} else {
			inner_it->second += 1;
		}
	}
}

property_map<edge_weight, RegionAdjacencyGraph::base_graph>::type& RegionAdjacencyGraph::get_edge_weight_map() {
	return get(edge_weight());
}

std::map<int, RegionAdjacencyGraph::Node>& RegionAdjacencyGraph::get_label_to_node_map() {
	return label_to_node_;
}


void RegionAdjacencyGraph::merge_nodes_threshold(double threshold, bool greedy) {
	property_map<edge_weight, RegionAdjacencyGraph::base_graph>::type& edge_weight_map = get(edge_weight());
	property_map<node_labels, RegionAdjacencyGraph::base_graph>::type& node_labels_map = get(node_labels());

	std::vector<RegionAdjacencyGraph::Edge> edges_to_contract;

	for (EdgeIt e(*this); e != lemon::INVALID; ++e) {
		if (edge_weight_map[e] < threshold) {
			continue;
		}
		edges_to_contract.push_back(e);
	}

	for (std::vector<RegionAdjacencyGraph::Edge>::const_iterator it = edges_to_contract.begin(); it != edges_to_contract.end(); ++it) {
		Edge e = *it;
		if (!this->valid(e)) { // might be invalid since cycles are removed by contract()
			LOG(logDEBUG) << "edge is invalid";
			continue;
		}
		Node u = this->u(e);
		Node v = this->v(e);
		if (!this->valid(u) || !this->valid(v)) {
			LOG(logDEBUG) << "one of the nodes is invalid!!!!";
			continue;
		}

		std::vector<int> label_list = node_labels_map[u];
		for (std::vector<int>::const_iterator lab_it = node_labels_map[v].begin(); lab_it != node_labels_map[v].end(); ++lab_it) {
			label_list.push_back(*lab_it);
		}
		node_labels_map.set(u, label_list);

		this->contract(u,v,greedy);

	}
}


std::vector<std::vector<int> > RegionAdjacencyGraph::get_labels_vector() {
	property_map<node_labels, RegionAdjacencyGraph::base_graph>::type& node_labels_map = get(node_labels());
	std::vector<std::vector<int> > result;
	for (NodeIt n(*this); n != lemon::INVALID; ++n) {
		result.push_back(node_labels_map[n]);
	}

	return result;
}




} /* namespace pgmlink */
