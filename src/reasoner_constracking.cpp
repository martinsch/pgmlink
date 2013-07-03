#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <opengm/inference/lpcplex.hxx>
#include <opengm/datastructures/marray/marray.hxx>

#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/traxels.h"

using namespace std;

namespace pgmlink {
ConservationTracking::~ConservationTracking() {
//	if (pgm_ != NULL) {
//		delete pgm_;
//		pgm_ = NULL;
//	}
//   if (optimizer_ != NULL) {
//      delete optimizer_;
//      optimizer_ = NULL;
//   }
}

double ConservationTracking::forbidden_cost() const {
  return forbidden_cost_;
}


void ConservationTracking::formulate(const HypothesesGraph& hypotheses) {
LOG(logDEBUG) << "ConservationTracking::formulate: entered";
reset();
pgm_ = boost::shared_ptr<pgm::OpengmModelDeprecated>(new pgm::OpengmModelDeprecated());

HypothesesGraph const *graph;
if (with_tracklets_) {
	LOG(logINFO) << "ConservationTracking::formulate: generating tracklet graph";
	tracklet2traxel_node_map_ = generateTrackletGraph2(hypotheses, tracklet_graph_);
	graph = &tracklet_graph_;
} else {
	graph = &hypotheses;
}

LOG(logDEBUG) << "ConservationTracking::formulate: add_transition_nodes";
add_transition_nodes(*graph);
LOG(logDEBUG) << "ConservationTracking::formulate: add_appearance_nodes";
add_appearance_nodes(*graph);
LOG(logDEBUG) << "ConservationTracking::formulate: add_disappearance_nodes";
add_disappearance_nodes(*graph);

LOG(logDEBUG) << "ConservationTracking::formulate: add_division_nodes";
if (with_divisions_) {
	add_division_nodes(*graph);
}
pgm::OpengmModelDeprecated::ogmGraphicalModel* model = pgm_->Model();

LOG(logDEBUG) << "ConservationTracking::formulate: add_finite_factors";
add_finite_factors(*graph);
LOG(logDEBUG) << "ConservationTracking::formulate: finished add_finite_factors";
typedef opengm::LPCplex<pgm::OpengmModelDeprecated::ogmGraphicalModel,pgm::OpengmModelDeprecated::ogmAccumulator> cplex_optimizer;
cplex_optimizer::Parameter param;
param.verbose_ = true;
param.integerConstraint_ = true;
param.epGap_ = ep_gap_;
LOG(logDEBUG) << "ConservationTracking::formulate ep_gap = " << param.epGap_;

optimizer_ = new cplex_optimizer(*model, param);

LOG(logDEBUG) << "ConservationTracking::formulate: add_constraints";
add_constraints(*graph);

LOG(logINFO) << "number_of_transition_nodes_ = " << number_of_transition_nodes_;
LOG(logINFO) << "number_of_appearance_nodes_ = " << number_of_appearance_nodes_;
LOG(logINFO) << "number_of_disappearance_nodes_ = " << number_of_disappearance_nodes_;
LOG(logINFO) << "number_of_division_nodes_ = " << number_of_division_nodes_;

}

void ConservationTracking::infer() {
opengm::InferenceTermination status = optimizer_->infer();
if (status != opengm::NORMAL) {
	throw std::runtime_error(
			"GraphicalModel::infer(): optimizer terminated abnormally");
}
}

void ConservationTracking::conclude(HypothesesGraph& g) {
// extract solution from optimizer
vector<pgm::OpengmModelDeprecated::ogmInference::LabelType> solution;
opengm::InferenceTermination status = optimizer_->arg(solution);
if (status != opengm::NORMAL) {
	throw runtime_error(
			"GraphicalModel::infer(): solution extraction terminated abnormally");
}

// add 'active' properties to graph
g.add(node_active2()).add(arc_active()).add(division_active());

property_map<node_active2, HypothesesGraph::base_graph>::type& active_nodes =
		g.get(node_active2());
property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs =
		g.get(arc_active());
property_map<division_active, HypothesesGraph::base_graph>::type& division_nodes =
			g.get(division_active());
if (!with_tracklets_) {
	tracklet_graph_.add(tracklet_intern_arc_ids()).add(traxel_arc_id());
}
property_map<tracklet_intern_arc_ids, HypothesesGraph::base_graph>::type& tracklet_arc_id_map =
		tracklet_graph_.get(tracklet_intern_arc_ids());
property_map<traxel_arc_id, HypothesesGraph::base_graph>::type& traxel_arc_id_map =
			tracklet_graph_.get(traxel_arc_id());

//	vector<size_t> count_objects(max_number_objects_+1,0);

for (HypothesesGraph::ArcIt a(g); a!=lemon::INVALID; ++a) {
	active_arcs.set(a, false);
}

// write state after inference into 'active'-property maps
// the node is also active if its appearance node is active
for (std::map<HypothesesGraph::Node, size_t>::const_iterator it =
		app_node_map_.begin(); it != app_node_map_.end(); ++it) {
	if (with_tracklets_) {
		// set state of tracklet nodes
		std::vector<HypothesesGraph::Node> traxel_nodes = tracklet2traxel_node_map_[it->first];
		for(std::vector<HypothesesGraph::Node>::const_iterator tr_n_it = traxel_nodes.begin(); tr_n_it != traxel_nodes.end(); ++tr_n_it) {
//				++count_objects[solution[it->second]];
			HypothesesGraph::Node n = *tr_n_it;
			active_nodes.set(n, solution[it->second]);
		}
		// set state of tracklet internal arcs
		std::vector<int> arc_ids = tracklet_arc_id_map[it->first];
		for(std::vector<int>::const_iterator arc_id_it = arc_ids.begin(); arc_id_it != arc_ids.end(); ++arc_id_it) {
			HypothesesGraph::Arc a = g.arcFromId(*arc_id_it);
			assert(active_arcs[a] == false);
			if(solution[it->second] > 0) {
				active_arcs.set(a, true);
				assert(active_nodes[g.source(a)] == solution[it->second] && "tracklet internal arcs must have the same flow as their connected nodes");
				assert(active_nodes[g.target(a)] == solution[it->second] && "tracklet internal arcs must have the same flow as their connected nodes");
			}
		}
	} else {
//			++count_objects[solution[it->second]];
		active_nodes.set(it->first, solution[it->second]);
	}
}

// the node is also active if its disappearance node is active
for (std::map<HypothesesGraph::Node, size_t>::const_iterator it =
		dis_node_map_.begin(); it != dis_node_map_.end(); ++it) {

	if (solution[it->second] > 0) {
		if (with_tracklets_) {
			// set state of tracklet nodes
			std::vector<HypothesesGraph::Node> traxel_nodes = tracklet2traxel_node_map_[it->first];
			for(std::vector<HypothesesGraph::Node>::const_iterator tr_n_it = traxel_nodes.begin(); tr_n_it != traxel_nodes.end(); ++tr_n_it) {
				HypothesesGraph::Node n = *tr_n_it;
				if (active_nodes[n] == 0) {
//						++count_objects[solution[it->second]];
//						--count_objects[0];
					active_nodes.set(n, solution[it->second]);
				} else {
					assert(active_nodes[n] == solution[it->second]);
				}
			}
			// set state of tracklet internal arcs
			std::vector<int> arc_ids = tracklet_arc_id_map[it->first];
			for(std::vector<int>::const_iterator arc_id_it = arc_ids.begin(); arc_id_it != arc_ids.end(); ++arc_id_it) {
				HypothesesGraph::Arc a = g.arcFromId(*arc_id_it);
				if (solution[it->second] > 0) {
					active_arcs.set(a, true);
					assert(active_nodes[g.source(a)] == solution[it->second] && "tracklet internal arcs must have the same flow as their connected nodes");
					assert(active_nodes[g.target(a)] == solution[it->second] && "tracklet internal arcs must have the same flow as their connected nodes");
				}
			}
		} else {
			if (active_nodes[it->first] == 0) {
//					++count_objects[solution[it->second]];
//					--count_objects[0];
				active_nodes.set(it->first, solution[it->second]);
			} else {
				assert(active_nodes[it->first] == solution[it->second]);
			}
		}
	}
}

for (std::map<HypothesesGraph::Arc, size_t>::const_iterator it =
		arc_map_.begin(); it != arc_map_.end(); ++it) {
	if (solution[it->second] >= 1) {
		if (with_tracklets_) {
			active_arcs.set(g.arcFromId((traxel_arc_id_map[it->first])), true);
		} else {
			active_arcs.set(it->first, true);
		}
	}
}
// initialize division node map
if (with_divisions_) {
	for (std::map<HypothesesGraph::Node, size_t>::const_iterator it =
			div_node_map_.begin(); it != div_node_map_.end(); ++it) {
		division_nodes.set(it->first, false);
	}
	for (std::map<HypothesesGraph::Node, size_t>::const_iterator it =
			div_node_map_.begin(); it != div_node_map_.end(); ++it) {
		if (solution[it->second] >=1) {
			if (with_tracklets_) {
				// set division property for the last node in the tracklet
				division_nodes.set(tracklet2traxel_node_map_[it->first].back(), true);
			} else {
				division_nodes.set(it->first, true);
			}
		}
	}
}

//	LOG(logINFO) << "ConservationTracking::conclude: number of objects in node:";
//	for (size_t i = 0; i<=max_number_objects_; ++i) {
//		LOG(logINFO) << "   " << i << ": " << count_objects[i];
//	}
}


const std::map<HypothesesGraph::Arc, size_t>& ConservationTracking::get_arc_map() const {
return arc_map_;
}

void ConservationTracking::reset() {
if (optimizer_ != NULL) {
	delete optimizer_;
	optimizer_ = NULL;
}
arc_map_.clear();
div_node_map_.clear();
app_node_map_.clear();
dis_node_map_.clear();
}


void ConservationTracking::add_appearance_nodes(const HypothesesGraph& g) {
size_t count = 0;
for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
	bool hasOutarcs = false;
	for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
		hasOutarcs = true;
		break;
	}
	if (!hasOutarcs) {
		continue;
	}
	pgm_->Model()->addVariable(max_number_objects_+1);
	app_node_map_[n] = pgm_->Model()->numberOfVariables() - 1;
	assert(pgm_->Model()->numberOfLabels(app_node_map_[n]) == max_number_objects_ + 1);
	LOG(logDEBUG4) << "appearance node added with id " << app_node_map_[n];
	++count;
}
number_of_appearance_nodes_ = count;
}

void ConservationTracking::add_disappearance_nodes(const HypothesesGraph& g) {
size_t count = 0;
for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
	bool hasInarcs = false;
	for (HypothesesGraph::InArcIt a(g, n); a != lemon::INVALID; ++a) {
		hasInarcs = true;
		break;
	}
	if (!hasInarcs) {
		bool hasOutarcs = false;
		for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
			hasOutarcs = true;
			break;
		}
		if (hasOutarcs) {
			continue;
		}
		// else: add the disappearance variable, otherwise single nodes/tracklets without
		// incoming and outgoing arcs won't be treated as variables
	}
	pgm_->Model()->addVariable(max_number_objects_+1);
	dis_node_map_[n] = pgm_->Model()->numberOfVariables() - 1;
	assert(pgm_->Model()->numberOfLabels(dis_node_map_[n]) == max_number_objects_ + 1);
	LOG(logDEBUG4) << "disappearance node added with id " << dis_node_map_[n];
	++count;
}
number_of_disappearance_nodes_ = count;
}

void ConservationTracking::add_transition_nodes(const HypothesesGraph& g) {
size_t count = 0;
for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
	pgm_->Model()->addVariable(max_number_objects_+1);
	arc_map_[a] = pgm_->Model()->numberOfVariables() - 1;
	assert(pgm_->Model()->numberOfLabels(arc_map_[a]) == max_number_objects_ + 1);
	LOG(logDEBUG4) << "transition node added with id " << arc_map_[a];
	++count;
}
number_of_transition_nodes_ = count;
}

void ConservationTracking::add_division_nodes(const HypothesesGraph& g) {
size_t count = 0;
for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
	size_t number_of_outarcs = 0;
	for(HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
		++number_of_outarcs;
	}
	if (number_of_outarcs > 1) {
		pgm_->Model()->addVariable(2);
		div_node_map_[n] = pgm_->Model()->numberOfVariables() - 1;
		assert(pgm_->Model()->numberOfLabels(div_node_map_[n]) == 2);
		LOG(logDEBUG4) << "division node added with id " << div_node_map_[n];
		++count;
	}
}
number_of_division_nodes_ = count;
}

namespace{
double get_transition_prob(double distance, size_t state, double alpha) {
double prob = exp(- distance / alpha);
if (state == 0) {
	return 1 - prob;
}
return prob;
}
}
void ConservationTracking::add_finite_factors(const HypothesesGraph& g) {
LOG(logDEBUG) << "ConservationTracking::add_finite_factors: entered";
property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = g.get(node_tracklet());
property_map<tracklet_intern_dist, HypothesesGraph::base_graph>::type& tracklet_intern_dist_map = g.get(tracklet_intern_dist());

////
//// add detection factors
////
LOG(logDEBUG) << "ConservationTracking::add_finite_factors: add detection factors";
for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
	size_t num_vars = 0;
	vector<size_t> vi;
	vector<double> cost;

	bool hasOutarcs = false;
	for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
		hasOutarcs = true;
		break;
	}
	bool hasInarcs = false;
	for (HypothesesGraph::InArcIt a(g, n); a != lemon::INVALID; ++a) {
		hasInarcs = true;
		break;
	}
	int node_begin_time = -1;
	int node_end_time = -1;
	if (with_tracklets_) {
		node_begin_time = tracklet_map[n].front().Timestep;
		node_end_time = tracklet_map[n].back().Timestep;
	} else {
		node_begin_time = traxel_map[n].Timestep;
		node_end_time = traxel_map[n].Timestep;
	}
	bool singleNodeTrack = !hasInarcs && !hasOutarcs;

	if (app_node_map_.count(n) > 0) {
		vi.push_back(app_node_map_[n]);
		if (node_begin_time == g.earliest_timestep()) {
			// pay no appearance costs in the first timestep
			cost.push_back(0.);
		} else {
			if(with_tracklets_){
				cost.push_back(appearance_cost_(tracklet_map[n].front()));
			} else {
				cost.push_back(appearance_cost_(traxel_map[n]));
			}
		}
		++num_vars;
	}
	if (dis_node_map_.count(n) > 0) {
		vi.push_back(dis_node_map_[n]);
		double c = 0;
		if (singleNodeTrack) {
			// single node tracks only have a disappearance node but not an appearance node
			// also consider tracklets here which can reach over multiple time steps
			if (node_end_time != g.latest_timestep()) {
				if(with_tracklets_){
					c += disappearance_cost_(tracklet_map[n].back());
				}else {
					c += disappearance_cost_(traxel_map[n]);
				}
			}
			if (node_begin_time != g.earliest_timestep()) {
				// single node tracks do not have an appearance node, hence handle appearance here
				if(with_tracklets_){
					c += appearance_cost_(tracklet_map[n].front());
				}else {
					c += appearance_cost_(traxel_map[n]);
				}
			}
		} else if (node_end_time != g.latest_timestep()) {
			if(with_tracklets_){
				c += disappearance_cost_(tracklet_map[n].back());
			}else {
				c += disappearance_cost_(traxel_map[n]);
			}
		}
		cost.push_back(c);
		++num_vars;
	}

	// convert vector to array
	vector<size_t> coords(num_vars, 0); // number of variables
	// ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
	pgm::OpengmExplicitFactor<double> table(vi.begin(), vi.end(), 0, (max_number_objects_ + 1));
	for (size_t state = 0; state <= max_number_objects_; ++state) {
		double energy = 0;
		if (with_tracklets_) {
			// add all detection factors of the internal nodes
			for (std::vector<Traxel>::const_iterator trax_it = tracklet_map[n].begin(); trax_it != tracklet_map[n].end(); ++trax_it){
				energy += detection_(*trax_it, state);
			}
			// add all transition factors of the internal arcs
			for (std::vector<double>::const_iterator intern_dist_it = tracklet_intern_dist_map[n].begin();
					intern_dist_it != tracklet_intern_dist_map[n].end(); ++intern_dist_it) {
				energy += transition_(get_transition_prob(*intern_dist_it, state, transition_parameter_));
			}
		} else {
			energy = detection_(traxel_map[n], state);
		}
		LOG(logDEBUG2) << "ConservationTracking::add_finite_factors: detection[" << state <<
						"] = " << energy;
		for (size_t var_idx = 0; var_idx < num_vars; ++var_idx) {
			coords[var_idx] = state;
			// if only one of the variables is > 0, then it is an appearance in this time frame
			// or a disappearance in the next timeframe. Hence, add the cost of appearance/disappearance
			// to the detection cost
			table.set_value(coords, energy+state*cost[var_idx]);
			coords[var_idx] = 0;
			LOG(logDEBUG4) << "ConservationTracking::add_finite_factors: var_idx " << var_idx <<
								" = " << energy;
		}
		// also this energy if both variables have the same state
		if (num_vars == 2) {
			coords[0] = state;
			coords[1] = state;
			// only pay detection energy if both variables are on
			table.set_value(coords, energy);
			coords[0] = 0;
			coords[1] = 0;

			LOG(logDEBUG4) << "ConservationTracking::add_finite_factors: var_idxs 0 and var_idx 1 = " << energy;
		}
	}

	LOG(logDEBUG3) << "ConservationTracking::add_finite_factors: adding table to pgm";
	table.add_to(*(pgm_->Model()));
}

////
//// add transition factors
////
LOG(logDEBUG) << "ConservationTracking::add_finite_factors: add transition factors";
property_map<arc_distance, HypothesesGraph::base_graph>::type& arc_distances = g.get(arc_distance());
for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
	size_t vi[] = { arc_map_[a] };
	vector<size_t> coords(1, 0); // number of variables
	// ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
	pgm::OpengmExplicitFactor<double> table(vi, vi + 1, 0, (max_number_objects_ + 1));
	for (size_t state = 0; state <= max_number_objects_; ++state) {
		double energy = transition_(get_transition_prob(arc_distances[a], state, transition_parameter_));
		LOG(logDEBUG2) << "ConservationTracking::add_finite_factors: transition[" << state <<
				"] = " << energy;
		coords[0] = state;
		table.set_value(coords, energy);
		coords[0] = 0;
	}
	table.add_to(*(pgm_->Model()));
}

////
//// add division factors
////
if (with_divisions_) {
	LOG(logDEBUG) << "ConservationTracking::add_finite_factors: add division factors";
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		if (div_node_map_.count(n) == 0) {
			continue;
		}
		size_t vi[] = { div_node_map_[n] };
		vector<size_t> coords(1, 0); // number of variables
		// ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
		pgm::OpengmExplicitFactor<double> table(vi, vi + 1, 0, 2);
		for (size_t state = 0; state <= 1; ++state) {
			double energy = 0;
			if (with_tracklets_) {
				energy = division_(tracklet_map[n].back(), state);
			} else {
				energy = division_(traxel_map[n], state);
			}
			LOG(logDEBUG2) << "ConservationTracking::add_finite_factors: division[" << state <<
					"] = " << energy;
			coords[0] = state;
			table.set_value(coords, energy);
			coords[0] = 0;
		}
		table.add_to(*(pgm_->Model()));
	}
}
}

size_t ConservationTracking::cplex_id(size_t opengm_id, size_t state) {
size_t number_of_nodes_with_multiple_states = number_of_transition_nodes_+
		number_of_appearance_nodes_ + number_of_disappearance_nodes_;
if (opengm_id <= (number_of_nodes_with_multiple_states)) {
	return (max_number_objects_ + 1) * opengm_id + state;
} else if (opengm_id <= (number_of_nodes_with_multiple_states + number_of_division_nodes_)) {
	return (max_number_objects_ + 1) * (number_of_nodes_with_multiple_states) +
			2 * (opengm_id - number_of_nodes_with_multiple_states) + state;
}
throw std::runtime_error("cplex_id(): open_gm id does not exist");
}


void ConservationTracking::add_constraints(const HypothesesGraph& g) {
LOG(logDEBUG) << "ConservationTracking::add_constraints: entered";
typedef opengm::LPCplex<pgm::OpengmModelDeprecated::ogmGraphicalModel,pgm::OpengmModelDeprecated::ogmAccumulator> cplex;


property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = g.get(node_tracklet());

LOG(logDEBUG) << "ConservationTracking::add_constraints: transitions";
for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
	std::stringstream traxel_names_ss;
	for (std::vector<Traxel>::const_iterator trax_it = tracklet_map[n].begin(); trax_it != tracklet_map[n].end(); ++trax_it){
		traxel_names_ss << trax_it->Id << " ";
	}
	std::string traxel_names = traxel_names_ss.str();



	vector<size_t> cplex_idxs, cplex_idxs2;
	vector<int> coeffs, coeffs2;

	////
	//// outgoing transitions
	////
	size_t num_outarcs = 0;
	// couple detection and transitions: Y_ij <= App_i
	for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
		assert(app_node_map_.count(n) > 0 && "this node should be contained in app_node_map_ since it has outgoing arcs");
		for (size_t nu = 0; nu < max_number_objects_; ++nu) {
			for (size_t mu = nu+1; mu <= max_number_objects_; ++mu) {
				cplex_idxs.clear();
				coeffs.clear();
				coeffs.push_back(1);
				cplex_idxs2.push_back(cplex_id(app_node_map_[n],nu));
				coeffs.push_back(1);
				cplex_idxs.push_back(cplex_id(arc_map_[a],mu));
				// 0 <= App_i[nu] + Y_ij[mu] <= 1  forall mu>nu
				optimizer_->addConstraint(cplex_idxs.begin(),
							cplex_idxs.end(), coeffs.begin(), 0, 1);
				LOG(logDEBUG3) << "ConservationTracking::add_constraints:" <<
						" Y_ij <= App_i added for Traxel " << traxel_names << ", "
						<< "n = " << app_node_map_[n] << ", a = " << arc_map_[a] << ", nu = " << nu << ", mu = " << mu;
			}
		}
		++num_outarcs;
	}


	int div_cplex_id = -1;
	if (with_divisions_ && div_node_map_.count(n) > 0) {
		LOG(logDEBUG3) << "div_node_map_[n] = " << div_node_map_[n];
		LOG(logDEBUG3) << "number_of_transition_nodes_ = " << number_of_transition_nodes_;
		LOG(logDEBUG3) << "number_of_division_nodes_ = " << number_of_division_nodes_;
		div_cplex_id = cplex_id(div_node_map_[n], 1);
	}
	if (num_outarcs > 0) {
		// couple transitions: sum(Y_ij) = D_i + App_i
		cplex_idxs.clear();
		coeffs.clear();
		for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
			for(size_t nu = 1; nu <= max_number_objects_; ++nu) {
				coeffs.push_back(nu);
				cplex_idxs.push_back(cplex_id(arc_map_[a],nu));
			}
		}
		if (div_cplex_id != -1) {
			cplex_idxs.push_back(div_cplex_id);
			coeffs.push_back(-1);
		}
		for (size_t nu = 1; nu <= max_number_objects_; ++nu) {
			coeffs.push_back(-nu);
			cplex_idxs.push_back(cplex_id(app_node_map_[n],nu));
		}

		// 0 <= sum_nu [ sum_j( nu * Y_ij[nu] ) ] - [ sum_nu nu * X_i[nu] + D_i[1] + sum_nu nu * App_i[nu] ]<= 0
		optimizer_->addConstraint(cplex_idxs.begin(),
				cplex_idxs.end(), coeffs.begin(), 0, 0);
		LOG(logDEBUG3) << "ConservationTracking::add_constraints:" <<
				" sum(Y_ij) = D_i + App_i added for Traxel " << traxel_names << ", "	<< "n = " << app_node_map_[n];

	}

	if (div_cplex_id != -1) {
		// couple detection and division: D_i = 1 => App_i = 1
		assert(app_node_map_.count(n) > 0 && "this node should be contained in app_node_map_ since it may divide");
		cplex_idxs.clear();
		coeffs.clear();

		cplex_idxs.push_back(div_cplex_id);
		coeffs.push_back(1);

		cplex_idxs.push_back(cplex_id(app_node_map_[n],1));
		coeffs.push_back(-1);

		// -1 <= D_i[1] - App_i[1] <= 0
		optimizer_->addConstraint(cplex_idxs.begin(),
				cplex_idxs.end(), coeffs.begin(), -1, 0);
		LOG(logDEBUG3) << "ConservationTracking::add_constraints:" <<
				" D_i=1 => App_i =1 added for  Traxel " << traxel_names << ", " << "n = " << app_node_map_[n] << ", d = " << div_node_map_[n];


		// couple divsion and transition: D_1 = 1 => sum_k(Y_ik) = 2
		cplex_idxs2.clear();
		coeffs2.clear(); // -m <= 2 * D_i[1] - sum_j ( Y_ij[1] ) <= 0
		cplex_idxs2.push_back(div_cplex_id);
		coeffs2.push_back(2);

		for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
			for (size_t nu = 2; nu <= max_number_objects_; ++nu) {
				// D_i[1] = 1 => Y_ij[nu] = 0 forall nu > 1
				cplex_idxs.clear();
				coeffs.clear();
				cplex_idxs.push_back(div_cplex_id);
				coeffs.push_back(1);
				cplex_idxs.push_back(cplex_id(arc_map_[a], nu));
				coeffs.push_back(1);

				// 0 <= D_i[1] + Y_ij[nu] <= 1 forall nu>1
			   optimizer_->addConstraint(cplex_idxs.begin(),
						cplex_idxs.end(), coeffs.begin(), 0, 1);
				LOG(logDEBUG3) << "ConservationTracking::add_constraints:" <<
						" D_i=1 => Y_i[nu]=0 added for Traxel " << traxel_names << ", " << "d = " << div_node_map_[n] << ", y = " << arc_map_[a] << ", nu = " << nu;

			}

			cplex_idxs2.push_back(cplex_id(arc_map_[a], 1));
			coeffs2.push_back(-1);
		}

		// -m <= 2 * D_i[1] - sum_j (Y_ij[1]) <= 0
		optimizer_->addConstraint(cplex_idxs2.begin(),
				cplex_idxs2.end(), coeffs2.begin(), -int(max_number_objects_), 0);
		LOG(logDEBUG3) << "ConservationTracking::add_constraints:" <<
				" D_i = 1 => sum_k(Y_ik) = 2 added for Traxel " << traxel_names << ", " << "d = " << div_node_map_[n];
	}


	////
	//// incoming transitions
	////
	// couple transitions: sum_k(Y_kj) = Dis_j
	cplex_idxs.clear();
	coeffs.clear();

	size_t num_inarcs = 0;
	for (HypothesesGraph::InArcIt a(g, n); a != lemon::INVALID; ++a) {
		for (size_t nu = 1; nu <= max_number_objects_; ++nu) {
			cplex_idxs.push_back(cplex_id(arc_map_[a],nu));
			coeffs.push_back(nu);
		}
		++num_inarcs;
	}

	if (num_inarcs > 0) {
		assert(dis_node_map_.count(n) > 0 && "this node should be contained in dis_node_map_ since it has incoming arcs");
		for (size_t nu = 1; nu <= max_number_objects_; ++nu) {
			cplex_idxs.push_back(cplex_id(dis_node_map_[n],nu));
			coeffs.push_back(-nu);
		}

		// 0 <= sum_nu [ nu * sum_i (Y_ij[nu] ) ] - sum_nu ( nu * X_j[nu] ) - sum_nu ( nu * Dis_j[nu] ) <= 0
		optimizer_->addConstraint(cplex_idxs.begin(),
				cplex_idxs.end(), coeffs.begin(), 0, 0);
		LOG(logDEBUG3) << "ConservationTracking::add_constraints:" <<
				" sum_k(Y_kj) = Dis_j added for Traxel " << traxel_names << ", " << "n = " << dis_node_map_[n];
	}


	////
	//// disappearance/appearance coupling
	////
	if (app_node_map_.count(n) > 0 && dis_node_map_.count(n) > 0) {
		for (size_t nu = 1; nu <= max_number_objects_; ++nu) {
			cplex_idxs.clear();
			coeffs.clear();

			cplex_idxs.push_back(cplex_id(app_node_map_[n],nu));
			coeffs.push_back(1);

			cplex_idxs.push_back(cplex_id(dis_node_map_[n],nu));
			coeffs.push_back(-1);

			cplex_idxs.push_back(cplex_id(dis_node_map_[n],0));
			coeffs.push_back(-1);

			// A_i[nu] = 1 => V_i[nu] = 1 v V_i[0] = 1
			// -1 <= App_i[nu] - ( Dis_i[nu] + Dis_i[0] ) <= 0 forall nu > 0
			optimizer_->addConstraint(cplex_idxs.begin(),
					cplex_idxs.end(), coeffs.begin(), -1, 0);
			LOG(logDEBUG3) << "ConservationTracking::add_constraints:" <<
					" A_i[nu] = 1 => V_i[nu] = 1 v V_i[0] = 1 added for Traxel " << traxel_names << ", " << "n = " << app_node_map_[n];
		}

		for (size_t nu = 1; nu <= max_number_objects_; ++nu) {
			cplex_idxs.clear();
			coeffs.clear();

			cplex_idxs.push_back(cplex_id(dis_node_map_[n],nu));
			coeffs.push_back(1);

			cplex_idxs.push_back(cplex_id(app_node_map_[n],nu));
			coeffs.push_back(-1);

			cplex_idxs.push_back(cplex_id(app_node_map_[n],0));
			coeffs.push_back(-1);

			// V_i[nu] = 1 => A_i[nu] = 1 v A_i[0] = 1
			// -1 <= Dis_i[nu] - ( App_i[nu] + App_i[0] ) <= 0 forall nu > 0
			optimizer_->addConstraint(cplex_idxs.begin(),
					cplex_idxs.end(), coeffs.begin(), -1, 0);
			LOG(logDEBUG3) << "ConservationTracking::add_constraints:" <<
					" V_i[nu] = 1 => A_i[nu] = 1 v A_i[0] = 1 added for Traxel " << traxel_names << ", " << "n = " << app_node_map_[n];
		}
	}


	if (!with_misdetections_allowed_) {
	  cplex_idxs.clear();
	  coeffs.clear();
	  if (dis_node_map_.count(n) > 0) {
		  cplex_idxs.push_back(cplex_id(dis_node_map_[n],0));
		  coeffs.push_back(1);
	  }
	  if (app_node_map_.count(n) > 0) {
		  cplex_idxs.push_back(cplex_id(app_node_map_[n],0));
		  coeffs.push_back(1);
	  }

	  // V_i[0] = 0 => 1 <= A_i[0]
	  // A_i[0] = 0 => 1 <= V_i[0]
	  // V_i <= m, A_i <= m
	  // 0 <= Dis_i[0] + App_i[0] <= 0
	  optimizer_->addConstraint(cplex_idxs.begin(),
								cplex_idxs.end(), coeffs.begin(), 0, 0);
	  LOG(logDEBUG3) << "ConservationTracking::add_constraints:" <<
						" A_i[0] + V_i[0] = 0 added for Traxel " << traxel_names;
	}

	if (!with_appearance_ && (dis_node_map_.count(n) > 0)) {
		cplex_idxs.clear();
		coeffs.clear();
		cplex_idxs.push_back(cplex_id(dis_node_map_[n],0));
		coeffs.push_back(1);
		// V_i[0] != 0
		// 1 <= V_i[0] <= m
		optimizer_->addConstraint(cplex_idxs.begin(),
									cplex_idxs.end(), coeffs.begin(), 1, max_number_objects_ );
		  LOG(logDEBUG3) << "ConservationTracking::add_constraints:" <<
							" V_i[0] => 1 added for Traxel " << traxel_names << ", " << "n = " << dis_node_map_[n];
	}

	if (!with_appearance_ && (app_node_map_.count(n) > 0)) {
		cplex_idxs.clear();
		coeffs.clear();
		cplex_idxs.push_back(cplex_id(app_node_map_[n],0));
		coeffs.push_back(1);
		// A_i[0] != 0
		// 1 <= A_i[0] <= m
		optimizer_->addConstraint(cplex_idxs.begin(),
									cplex_idxs.end(), coeffs.begin(), 1, max_number_objects_ );
		  LOG(logDEBUG3) << "ConservationTracking::add_constraints:" <<
							" A_i[0] => 1 added for Traxel " << traxel_names << ", " << "n = " << app_node_map_[n];
	}
  }


}






//void ConservationTracking::fix_detections(const HypothesesGraph& g) {
////   typedef opengm::LPCplex<pgm::OpengmModelDeprecated, opengm::Minimizer> cplex;
////	typedef opengm::LPCplex<OpengmModelDeprecated::ogmGraphicalModel,OpengmModelDeprecated::ogmAccumulator> cplex;
//	typedef opengm::LPCplex<pgm::OpengmModelDeprecated::ogmGraphicalModel,pgm::OpengmModelDeprecated::ogmAccumulator> cplex;
//	property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
//	property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = g.get(node_tracklet());
//	property_map<tracklet_intern_dist, HypothesesGraph::base_graph>::type& tracklet_intern_dist_map = g.get(tracklet_intern_dist());
//
//	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
//		vector<size_t> cplex_idxs;
//		vector<int> coeffs;
//		vector<double> energies;
//		for (size_t state = 0; state <= max_number_objects_; ++state) {
//			if (with_tracklets_) {
//				double e = 0;
//				// add all detection factors of the internal nodes
//				for (std::vector<Traxel>::const_iterator trax_it = tracklet_map[n].begin(); trax_it != tracklet_map[n].end(); ++trax_it){
//					e += detection_(*trax_it, state);
//				}
//				// add all transition factors of the internal arcs
//				for (std::vector<double>::const_iterator intern_dist_it = tracklet_intern_dist_map[n].begin();
//						intern_dist_it != tracklet_intern_dist_map[n].end(); ++intern_dist_it) {
//					e += transition_(get_transition_prob(*intern_dist_it, state));
//				}
//				energies.push_back(e);
//			} else {
//				energies.push_back((double) detection_(traxel_map[n], state));
//			}
//		}
//		size_t max_state = std::min_element(energies.begin(), energies.end()) - energies.begin();
//
//		// if it is most probable that it is a misdetection (max_state == 0), then fix all X_i, App_i, Dis_i to 0
//		// otherwise, add constraint that one of those must be greater than 0
//		// (Detections cannot be fixed to their max_state, example: 1 --- 2 --- 1 would be invalid)
//
//		if (max_state == 0) {
////			cplex_idxs.push_back(cplex_id(node_map_[n], 0));
////			coeffs.push_back(1);
////			// X_i[0] = 1
////			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
////					cplex_idxs.end(), coeffs.begin(), 1, 1);
////			LOG(logDEBUG3) << "ConservationTracking::fix_detections:" <<
////								" X_i set to 0 for n = " << node_map_[n];
////			if (with_appearance_) {
////				cplex_idxs.clear();
////				coeffs.clear();
////				cplex_idxs.push_back(cplex_id(app_node_map_[n],0));
////				coeffs.push_back(1);
////				// App_i[0] = 1
////				dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
////						cplex_idxs.end(), coeffs.begin(), 1, 1);
////				LOG(logDEBUG3) << "ConservationTracking::fix_detections:" <<
////									" App_i set to 0 for n = " << node_map_[n];
////			}
////			if (with_disappearance_) {
////				cplex_idxs.clear();
////				coeffs.clear();
////				cplex_idxs.push_back(cplex_id(dis_node_map_[n],0));
////				coeffs.push_back(1);
////				// Dis_i[0] = 1
////				dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
////						cplex_idxs.end(), coeffs.begin(), 1, 1);
////				LOG(logDEBUG3) << "ConservationTracking::fix_detections:" <<
////									" Dis_i set to 0 for n = " << node_map_[n];
////			}
//
////		}
//		} else {  // max_state > 0
//			// one of X_i, App_i or Dis_i must be greater than 0
//			cplex_idxs.push_back(cplex_id(node_map_[n], 0));
//			coeffs.push_back(1);
//			size_t num_vars = 1;
//
//			if (with_appearance_) {
//				cplex_idxs.push_back(cplex_id(app_node_map_[n],0));
//				coeffs.push_back(1);
//				++num_vars;
//			}
//			if (with_disappearance_) {
//				cplex_idxs.push_back(cplex_id(dis_node_map_[n],0));
//				coeffs.push_back(1);
//				++num_vars;
//			}
//			// 2 <= 1*X_i[0] + 1*App_i[0] + 1*Dis_i[0] <= 2
////			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
//			optimizer_->addConstraint(cplex_idxs.begin(),
//					cplex_idxs.end(), coeffs.begin(), num_vars-1, num_vars-1);
//			LOG(logDEBUG3) << "ConservationTracking::fix_detections:" <<
//								" X_i != 0 added for " << "n = " << node_map_[n];
//		}
//	}
//}


} /* namespace pgmlink */
