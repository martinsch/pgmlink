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

	create_region(node, false);

	return node;
}

int RegionAdjacencyGraph::find_region(std::vector<int> labels_sorted) {
	for(std::vector<Region>::iterator it = regions_.begin(); it != regions_.end(); ++it ) {
		if (it->contains_labels.size() != labels_sorted.size()) {
			continue;
		}
		std::vector<int>::const_iterator it1 = it->contains_labels.begin();
		std::vector<int>::const_iterator it2 = labels_sorted.begin();
		int count = 0;
		for( ; it2 != labels_sorted.end(); ++it1, ++it2 ) {
			assert(*it1 <= segment_max_id_ && "a region can only contain segments, not other regions (6)");
			assert(*it2 <= segment_max_id_ && "a region can only contain segments, not other regions (7)");
			if ( *it1 != *it2 ) {
				break;
			}
			++count;
		}
		// if the two maps are the same, then count == labels.size(), i.e.
		// the region is already contained in regions_ and its id will be returned
		if (count == labels_sorted.size()) {
			return it->id;
		}
	}
	return -1;
}


int RegionAdjacencyGraph::create_region(RegionAdjacencyGraph::Node n, bool with_check) {
	property_map<node_labels, RegionAdjacencyGraph::base_graph>::type& node_labels_map = get(node_labels());
	std::vector<int> labels = node_labels_map[n];
	std::sort(labels.begin(), labels.end());

	property_map<node_to_region_id, RegionAdjacencyGraph::base_graph>::type& node_to_region_id_map = get(node_to_region_id());
	if (with_check) { // check whether the region already exists
		int region_found = find_region(labels);
		if (region_found != -1 ) {
			LOG(logDEBUG3) << "region found, id = " << region_found;
			node_to_region_id_map.set(n, region_found);
			return region_found;
		}
	}

	Region reg = Region();
	reg.id = ++region_max_id_;
	LOG(logDEBUG3) << "create region with id = " << reg.id;
	for (std::vector<int>::const_iterator it = labels.begin(); it != labels.end(); ++it) {
		assert(*it <= segment_max_id_ && "a region can only contain segments, not other regions (1)");
		reg.contains_labels.push_back(*it);
	}
	// sort the vector
	std::sort(reg.contains_labels.begin(), reg.contains_labels.end());
	regions_.push_back(reg);
	assert(reg.id == (int) regions_.size());

	node_to_region_id_map.set(n, reg.id);

	return reg.id;
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

void RegionAdjacencyGraph::increment_perimeter(RegionAdjacencyGraph::Node n) {
	std::map<RegionAdjacencyGraph::Node, int>::iterator it = node_perimeters_.find(n);
	if( it == node_perimeters_.end() ) {
		node_perimeters_[n] = 1;
		return;
	}
	it->second += 1;
}

void RegionAdjacencyGraph::increment_intersection(RegionAdjacencyGraph::Node n1, RegionAdjacencyGraph::Node n2) {
	if (this->id(n2) < this->id(n1) ) { // only store node_intersections from smaller to bigger id to prevent redundancy
		increment_intersection(n2, n1);
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

bool RegionAdjacencyGraph::is_connected_component(RegionAdjacencyGraph::Node n) {
	IncEdgeIt e(*this, n);
	if (e == lemon::INVALID) {
		return true;
	}
	return false;
}

std::vector<int> RegionAdjacencyGraph::get_connected_component_ids() {
	assert(lemon::countEdges(*this) == 0 && "to get the connected components, all edges must have been contracted before:");
	return connected_component_region_ids_;
}

std::vector<Region> RegionAdjacencyGraph::get_regions() {
	return regions_;
}

void RegionAdjacencyGraph::contract_nodes(RegionAdjacencyGraph::Node a, RegionAdjacencyGraph::Node b, bool remove_loops ) {
	property_map<node_labels, RegionAdjacencyGraph::base_graph>::type& node_labels_map = get(node_labels());

	std::vector<int> label_list = node_labels_map[a];
	for (std::vector<int>::const_iterator lab_it = node_labels_map[b].begin(); lab_it != node_labels_map[b].end(); ++lab_it) {
		label_list.push_back(*lab_it);
	}
	node_labels_map.set(a, label_list);

	this->contract(a,b,remove_loops); // a and b are contracted to a
}

void RegionAdjacencyGraph::merge_nodes_threshold(double threshold, bool greedy) {
	LOG(logDEBUG) << "merge_nodes_threshold entered";
	// remove connected components
	std::vector<RegionAdjacencyGraph::Node> nodes_to_remove;
	for (NodeIt n(*this); n != lemon::INVALID; ++n) {
		if (is_connected_component(n)) {
			nodes_to_remove.push_back(n);
		}
	}
	int num_edges = lemon::countEdges(*this);
	LOG(logDEBUG1) << "removing " << nodes_to_remove.size() << " nodes (connected components)";
	for (std::vector<RegionAdjacencyGraph::Node>::const_iterator it = nodes_to_remove.begin(); it != nodes_to_remove.end(); ++it) {
		assert(is_connected_component(*it));
		this->erase(*it);
	}
	assert(num_edges == lemon::countEdges(*this));

	property_map<edge_weight, RegionAdjacencyGraph::base_graph>::type& edge_weight_map = get(edge_weight());

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

		contract_nodes(u,v,true); // true, if loop edges should be removed (i.e. edges where u=v)

	}

	// create regions (but check whether they exist already)
	for (NodeIt n(*this); n != lemon::INVALID; ++n) {
		create_region(n, true);
	}

	// mark new connected components
	property_map<node_to_region_id, RegionAdjacencyGraph::base_graph>::type& node_to_region_id_map = get(node_to_region_id());
	for(NodeIt n(*this); n != lemon::INVALID; ++n) {
		if (is_connected_component(n)) {
			int reg_id = node_to_region_id_map[n];
			connected_component_region_ids_.push_back(reg_id);
		}
	}
}

std::vector<Region> RegionAdjacencyGraph::get_affected_regions(int cc_id) {
	std::vector<Region> affected_regions;
	Region cc = regions_[cc_id - 1];
	assert(cc_id == cc.id);
	for (std::vector<int>::const_iterator label_it = cc.contains_labels.begin(); label_it != cc.contains_labels.end(); ++label_it ) {
		int label = *label_it;
		assert(label <= segment_max_id_ && "a region can only contain segments, not other regions (2)");
		for( std::vector<Region>::const_iterator reg_it = regions_.begin(); reg_it != regions_.end(); ++reg_it ) {
			if (std::find(reg_it->contains_labels.begin(), reg_it->contains_labels.end(), label) != reg_it->contains_labels.end()) {
				affected_regions.push_back(*reg_it);
			}
		}
	}
	return affected_regions;
}

std::vector<std::vector<int> > RegionAdjacencyGraph::get_conflicts_cc(std::vector<Region> affected_regions, int cc_id) {
	std::vector<std::vector<int> > conflict_sets;

	Region cc = regions_[cc_id - 1];
	assert(cc_id == cc.id);
	for(std::vector<int>::const_iterator label_it = cc.contains_labels.begin(); label_it != cc.contains_labels.end(); ++label_it) {
		assert(*label_it <= segment_max_id_ && "a region can only contain segments, not other regions (3)");
		std::vector<int> conflict_set;
		for (std::vector<Region>::const_iterator reg_it = affected_regions.begin(); reg_it != affected_regions.end(); ++reg_it ) {
			if (std::find(reg_it->contains_labels.begin(), reg_it->contains_labels.end(), *label_it) != reg_it->contains_labels.end()) {
				// add to the conflict set if it is not already in it (this could be the case since affected_regions has duplicates,
				// e.g. the connected component is added to the affected_regions n times, since each of the n segments adds it.
				// the check whether it is already in the set, is cheaper here than when adding affected regions
				if (std::find(conflict_set.begin(), conflict_set.end(), reg_it->id) == conflict_set.end()) {
					conflict_set.push_back(reg_it->id);
				}
			}
		}

		conflict_sets.push_back(conflict_set);
	}

	return conflict_sets;
}

std::map<int, std::vector<std::vector<int> > > RegionAdjacencyGraph::get_conflict_sets() {
	std::map<int, std::vector<std::vector<int> > > result;
	for (std::vector<int>::const_iterator id_it = connected_component_region_ids_.begin(); id_it != connected_component_region_ids_.end(); ++id_it) {
		int cc_id = *id_it;
		std::vector<Region> affected_regions = get_affected_regions(cc_id);
		std::vector<std::vector<int> > conflict_sets = get_conflicts_cc(affected_regions, cc_id);
		result[cc_id] = conflict_sets;
	}
	return result;
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
