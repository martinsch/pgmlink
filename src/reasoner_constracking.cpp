#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <opengm/inference/lpcplex.hxx>
#include <opengm/datastructures/marray/marray.hxx>

#include "pgmlink/graphical_model.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/traxels.h"

using namespace std;

namespace Tracking {
SingleTimestepTraxelConservation::~SingleTimestepTraxelConservation() {
	if (pgm_ != NULL) {
		delete pgm_;
		pgm_ = NULL;
	}
	if (optimizer_ != NULL) {
		delete optimizer_;
		optimizer_ = NULL;
	}
}

double SingleTimestepTraxelConservation::forbidden_cost() const {
	return forbidden_cost_;
}

bool SingleTimestepTraxelConservation::with_constraints() const {
	return with_constraints_;
}

void SingleTimestepTraxelConservation::formulate(const HypothesesGraph& hypotheses) {
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: entered";
	reset();
	pgm_ = new OpengmModel();

	HypothesesGraph const *graph;
	if (with_tracklets_) {
		LOG(logINFO) << "SingleTimestepTraxelConservation::formulate: generating tracklet graph";
		tracklet2traxel_node_map_ = generateTrackletGraph2(hypotheses,tracklet_graph_);
		graph = &tracklet_graph_;
	} else {
		graph = &hypotheses;
	}
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_detection_nodes";
	add_detection_nodes(*graph);
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_transition_nodes";
	add_transition_nodes(*graph);
	if (with_appearance_) {
		LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_appearance_nodes";
		add_appearance_nodes(*graph);
	}
	if (with_disappearance_) {
		LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_disappearance_nodes";
		add_disappearance_nodes(*graph);
	}
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_division_nodes";
	add_division_nodes(*graph);

	OpengmModel::ogmGraphicalModel* model = pgm_->Model();

	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_finite_factors";
	add_finite_factors(*graph);

	typedef opengm::LPCplex<OpengmModel::ogmGraphicalModel,OpengmModel::ogmAccumulator> cplex_optimizer;
	cplex_optimizer::Parameter param;
	param.verbose_ = true;
	param.integerConstraint_ = true;
	param.epGap_ = ep_gap_;
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate ep_gap = " << param.epGap_;

	optimizer_ = new cplex_optimizer(*model, param);

	if (with_constraints_) {
		LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_constraints";
		add_constraints(*graph);
	}

	if (fixed_detections_) {
		LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: fix_detections";
		assert(with_appearance_ || with_disappearance_);
		fix_detections(*graph);
	}
}

void SingleTimestepTraxelConservation::infer() {
	opengm::InferenceTermination status = optimizer_->infer();
	if (status != opengm::NORMAL) {
		throw std::runtime_error(
				"GraphicalModel::infer(): optimizer terminated abnormally");
	}
}

void SingleTimestepTraxelConservation::conclude(HypothesesGraph& g) {
	// extract solution from optimizer
	vector<OpengmModel::ogmInference::LabelType> solution;
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

	vector<size_t> count_objects(max_number_objects_+1,0);

	for (HypothesesGraph::ArcIt a(g); a!=lemon::INVALID; ++a) {
		active_arcs.set(a, false);
	}

	// write state after inference into 'active'-property maps
	for (std::map<HypothesesGraph::Node, size_t>::const_iterator it =
			node_map_.begin(); it != node_map_.end(); ++it) {
		if (with_tracklets_) {
			// set state of tracklet nodes
			std::vector<HypothesesGraph::Node> traxel_nodes = tracklet2traxel_node_map_[it->first];
			for(std::vector<HypothesesGraph::Node>::const_iterator tr_n_it = traxel_nodes.begin(); tr_n_it != traxel_nodes.end(); ++tr_n_it) {
				HypothesesGraph::Node n = *tr_n_it;
				active_nodes.set(n, solution[it->second]);
				++count_objects[solution[it->second]];
			}
			// set state of tracklet internal arcs
			std::vector<int> arc_ids = tracklet_arc_id_map[it->first];
			LOG(logDEBUG) << "arc_ids.size() = " << arc_ids.size();
			LOG(logDEBUG) << "traxel_nodes.size() = " << traxel_nodes.size();
			assert(arc_ids.size() == traxel_nodes.size() -1);
			for(std::vector<int>::const_iterator arc_id_it = arc_ids.begin(); arc_id_it != arc_ids.end(); ++arc_id_it) {
				HypothesesGraph::Arc a = g.arcFromId(*arc_id_it);
				if (solution[it->second] > 0) {
					active_arcs.set(a, true);
				}
			}
		} else {
			++count_objects[solution[it->second]];
			active_nodes.set(it->first, solution[it->second]);
		}
	}
	// the node is also active if its appearance node is active
	for (std::map<HypothesesGraph::Node, size_t>::const_iterator it =
			app_node_map_.begin(); it != app_node_map_.end(); ++it) {
		if (solution[it->second] > 0) {
			if (with_tracklets_) {
				// set state of tracklet nodes
				std::vector<HypothesesGraph::Node> traxel_nodes = tracklet2traxel_node_map_[it->first];
				for(std::vector<HypothesesGraph::Node>::const_iterator tr_n_it = traxel_nodes.begin(); tr_n_it != traxel_nodes.end(); ++tr_n_it) {
					++count_objects[solution[it->second]];
					--count_objects[0];
					HypothesesGraph::Node n = *tr_n_it;
					assert(active_nodes[n] == 0);
					active_nodes.set(n, solution[it->second]);
				}
				// set state of tracklet internal arcs
				std::vector<int> arc_ids = tracklet_arc_id_map[it->first];
				for(std::vector<int>::const_iterator arc_id_it = arc_ids.begin(); arc_id_it != arc_ids.end(); ++arc_id_it) {
					HypothesesGraph::Arc a = g.arcFromId(*arc_id_it);
					assert(active_arcs[a] == false);
					if(solution[it->second] > 0) {
						active_arcs.set(a, true);
					}
				}
			} else {
				++count_objects[solution[it->second]];
				--count_objects[0];
				assert(active_nodes[it->first]==0);
				active_nodes.set(it->first, solution[it->second]);
			}
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
					++count_objects[solution[it->second]];
					--count_objects[0];
					HypothesesGraph::Node n = *tr_n_it;
					assert(active_nodes[n] == 0);
					active_nodes.set(n, solution[it->second]);
				}
				// set state of tracklet internal arcs
				std::vector<int> arc_ids = tracklet_arc_id_map[it->first];
				for(std::vector<int>::const_iterator arc_id_it = arc_ids.begin(); arc_id_it != arc_ids.end(); ++arc_id_it) {
					HypothesesGraph::Arc a = g.arcFromId(*arc_id_it);
					assert(active_arcs[a] == false);
					if (solution[it->second] > 0) {
						active_arcs.set(a, true);
					}
				}
			} else {
				++count_objects[solution[it->second]];
				--count_objects[0];
				assert(active_nodes[it->first]==0);
				active_nodes.set(it->first, solution[it->second]);
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
	for (std::map<HypothesesGraph::Node, size_t>::const_iterator it =
			node_map_.begin(); it != node_map_.end(); ++it) {
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

	LOG(logINFO) << "SingleTimestepTraxelConservation::conclude: number of objects in node:";
	for (size_t i = 0; i<=max_number_objects_; ++i) {
		LOG(logINFO) << "   " << i << ": " << count_objects[i];
	}
}

const OpengmModel* SingleTimestepTraxelConservation::get_graphical_model() const {
	return pgm_;
}

const std::map<HypothesesGraph::Node, size_t>& SingleTimestepTraxelConservation::get_node_map() const {
	return node_map_;
}

const std::map<HypothesesGraph::Arc, size_t>& SingleTimestepTraxelConservation::get_arc_map() const {
	return arc_map_;
}

void SingleTimestepTraxelConservation::reset() {
	if (pgm_ != NULL) {
		delete pgm_;
		pgm_ = NULL;
	}
	if (optimizer_ != NULL) {
		delete optimizer_;
		optimizer_ = NULL;
	}
	node_map_.clear();
	arc_map_.clear();
	div_node_map_.clear();
	app_node_map_.clear();
	dis_node_map_.clear();
}

void SingleTimestepTraxelConservation::add_detection_nodes(const HypothesesGraph& g) {
	size_t count = 0;
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		pgm_->Model()->addVariable(max_number_objects_+1);
		node_map_[n] = pgm_->Model()->numberOfVariables() - 1;
		assert(pgm_->Model()->numberOfLabels(node_map_[n]) == max_number_objects_ + 1);
		LOG(logDEBUG4) << "detection node added with id " << node_map_[n];
		++count;
	}
	number_of_detection_nodes_ = count;
}

void SingleTimestepTraxelConservation::add_appearance_nodes(const HypothesesGraph& g) {
	size_t count = 0;
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		pgm_->Model()->addVariable(max_number_objects_+1);
		app_node_map_[n] = pgm_->Model()->numberOfVariables() - 1;
		assert(pgm_->Model()->numberOfLabels(app_node_map_[n]) == max_number_objects_ + 1);
		LOG(logDEBUG4) << "appearance node added with id " << app_node_map_[n];
		++count;
	}
	number_of_appearance_nodes_ = count;
}

void SingleTimestepTraxelConservation::add_disappearance_nodes(const HypothesesGraph& g) {
	size_t count = 0;
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		pgm_->Model()->addVariable(max_number_objects_+1);
		dis_node_map_[n] = pgm_->Model()->numberOfVariables() - 1;
		assert(pgm_->Model()->numberOfLabels(dis_node_map_[n]) == max_number_objects_ + 1);
		LOG(logDEBUG4) << "detection node added with id " << dis_node_map_[n];
		++count;
	}
	number_of_disappearance_nodes_ = count;
}

void SingleTimestepTraxelConservation::add_transition_nodes(const HypothesesGraph& g) {
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

void SingleTimestepTraxelConservation::add_division_nodes(const HypothesesGraph& g) {
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
double get_transition_prob(double distance, size_t state) {
	double prob = exp(- distance / 5);
	if (state == 0) {
		return 1 - prob;
	}
	return prob;
}
}
void SingleTimestepTraxelConservation::add_finite_factors(const HypothesesGraph& g) {
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::add_finite_factors: entered";
	property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
	property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = g.get(node_tracklet());
	property_map<tracklet_intern_dist, HypothesesGraph::base_graph>::type& tracklet_intern_dist_map = g.get(tracklet_intern_dist());

	////
	//// add detection factors
	////
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::add_finite_factors: add detection factors";
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		size_t num_vars = 1;
		if (with_appearance_) ++num_vars;
		if (with_disappearance_) ++num_vars;
		vector<size_t> vi;
		vi.push_back(node_map_[n]);
		if (with_appearance_) {
			vi.push_back(app_node_map_[n]);
		}
		if (with_disappearance_) {
			vi.push_back(dis_node_map_[n]);
		}
		// convert vector to array
//		size_t* viarray = &vi[0];
		vector<size_t> coords(num_vars, 0); // number of variables
		// ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
		OpengmExplicitFactor<double> table(vi.begin(), vi.end(), 0, (max_number_objects_ + 1));
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
					energy += transition_(get_transition_prob(*intern_dist_it, state));
				}
			} else {
				energy = detection_(traxel_map[n], state);
			}
			LOG(logDEBUG2) << "SingleTimestepTraxelConservation::add_finite_factors: detection[" << state <<
							"] = " << energy;
			for (size_t var_idx = 0; var_idx < num_vars; ++var_idx) {
				coords[var_idx] = state;
				table.set_value(coords, energy);
				coords[var_idx] = 0;
				LOG(logDEBUG4) << "SingleTimestepTraxelConservation::add_finite_factors: var_idx " << var_idx <<
									" = " << energy;
			}
		}

		LOG(logDEBUG3) << "SingleTimestepTraxelConservation::add_finite_factors: adding table to pgm";
		table.add_to(*pgm_);
	}

	////
	//// add transition factors
	////
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::add_finite_factors: add transition factors";
	property_map<arc_distance, HypothesesGraph::base_graph>::type& arc_distances = g.get(arc_distance());
	for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
		size_t vi[] = { arc_map_[a] };
		vector<size_t> coords(1, 0); // number of variables
		// ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
		OpengmExplicitFactor<double> table(vi, vi + 1, 0, (max_number_objects_ + 1));
		for (size_t state = 0; state <= max_number_objects_; ++state) {
			double energy = transition_(get_transition_prob(arc_distances[a], state));
			LOG(logDEBUG2) << "SingleTimestepTraxelConservation::add_finite_factors: transition[" << state <<
					"] = " << energy;
			coords[0] = state;
			table.set_value(coords, energy);
			coords[0] = 0;
		}
		table.add_to(*pgm_);
	}

	////
	//// add division factors
	////
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::add_finite_factors: add division factors";
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		if (div_node_map_.count(n) == 0) {
			continue;
		}
		size_t vi[] = { div_node_map_[n] };
		vector<size_t> coords(1, 0); // number of variables
		// ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
		OpengmExplicitFactor<double> table(vi, vi + 1, 0, 2);
		for (size_t state = 0; state <= 1; ++state) {
			double energy = 0;
			if (with_tracklets_) {
				energy = division_(tracklet_map[n].back(), state);
			} else {
				energy = division_(traxel_map[n], state);
			}
			LOG(logDEBUG2) << "SingleTimestepTraxelConservation::add_finite_factors: division[" << state <<
					"] = " << energy;
			coords[0] = state;
			table.set_value(coords, energy);
			coords[0] = 0;
		}
		table.add_to(*pgm_);
	}
}

size_t SingleTimestepTraxelConservation::cplex_id(size_t opengm_id, size_t state) {
	size_t number_of_nodes_with_multiple_states = number_of_detection_nodes_ + number_of_transition_nodes_+
			number_of_appearance_nodes_ + number_of_disappearance_nodes_;
	if (opengm_id <= (number_of_nodes_with_multiple_states)) {
		return (max_number_objects_ + 1) * opengm_id + state;
	} else if (opengm_id <= (number_of_nodes_with_multiple_states + number_of_division_nodes_)) {
		return (max_number_objects_ + 1) * (number_of_nodes_with_multiple_states) +
				2 * (opengm_id - number_of_nodes_with_multiple_states) + state;
	}
	throw std::runtime_error("cplex_id(): open_gm id does not exist");
}


void SingleTimestepTraxelConservation::add_constraints(const HypothesesGraph& g) {
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::add_constraints: entered";
	typedef opengm::LPCplex<OpengmModel::ogmGraphicalModel,OpengmModel::ogmAccumulator> cplex;

	LOG(logDEBUG) << "SingleTimestepTraxelConservation::add_constraints: transitions";
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		vector<size_t> cplex_idxs, cplex_idxs2;
		vector<int> coeffs, coeffs2;

		////
		//// outgoing transitions
		////
		size_t num_outarcs = 0;
		// couple detection and transitions: Y_ij <= X_i + App_i
		for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
			for (size_t nu = 0; nu < max_number_objects_; ++nu) {
				for (size_t mu = nu+1; mu <= max_number_objects_; ++mu) {
					cplex_idxs.clear();
					coeffs.clear();
					coeffs.push_back(1);
					cplex_idxs.push_back(cplex_id(node_map_[n],nu));
					if (with_appearance_) {
						coeffs.push_back(1);
						cplex_idxs2.push_back(cplex_id(app_node_map_[n],nu));
					}
					coeffs.push_back(1);
					cplex_idxs.push_back(cplex_id(arc_map_[a],mu));
					// 0 <= X_i[nu] + App_i[nu] + Y_ij[mu] <= 1  forall mu>nu
					dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
								cplex_idxs.end(), coeffs.begin(), 0, 1);
					LOG(logDEBUG3) << "SingleTimestepTraxelConservation::add_constraints:" <<
							" Y_ij <= X_i + App_i added for "
							<< "n = " << node_map_[n] << ", a = " << arc_map_[a] << ", nu = " << nu << ", mu = " << mu;
				}
			}
			++num_outarcs;
		}


		int div_cplex_id = -1;
		if (div_node_map_.count(n) > 0) {
			LOG(logDEBUG3) << "div_node_map_[n] = " << div_node_map_[n];
			LOG(logDEBUG3) << "number_of_detection_nodes_ = " << number_of_detection_nodes_;
			LOG(logDEBUG3) << "number_of_transition_nodes_ = " << number_of_transition_nodes_;
			LOG(logDEBUG3) << "number_of_division_nodes_ = " << number_of_division_nodes_;
			div_cplex_id = cplex_id(div_node_map_[n], 1);
		}
		if (num_outarcs > 0) {
			// couple transitions: sum(Y_ij) = D_i + X_i + App_i
			cplex_idxs.clear();
			coeffs.clear();
			for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
				for(size_t nu = 1; nu <= max_number_objects_; ++nu) {
					coeffs.push_back(nu);
					cplex_idxs.push_back(cplex_id(arc_map_[a],nu));
				}
			}
			for (size_t nu = 1; nu <= max_number_objects_; ++nu) {
				coeffs.push_back(-nu);
				cplex_idxs.push_back(cplex_id(node_map_[n],nu));
			}
			if (div_cplex_id != -1) {
				cplex_idxs.push_back(div_cplex_id);
				coeffs.push_back(-1);
			}
			if (with_appearance_) {
				for (size_t nu = 1; nu <= max_number_objects_; ++nu) {
					coeffs.push_back(-nu);
					cplex_idxs.push_back(cplex_id(app_node_map_[n],nu));
				}
			}

			// 0 <= sum_nu [ sum_j( nu * Y_ij[nu] ) ] - [ sum_nu nu * X_i[nu] + D_i[1] + sum_nu nu * App_i[nu] ]<= 0
			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
					cplex_idxs.end(), coeffs.begin(), 0, 0);
			LOG(logDEBUG3) << "SingleTimestepTraxelConservation::add_constraints:" <<
					" sum(Y_ij) = D_i + X_i + App_i added for "	<< "n = " << node_map_[n];

		}

		if (div_cplex_id != -1) {
			// couple detection and division: D_i = 1 => X_i + App_i = 1
			cplex_idxs.clear();
			coeffs.clear();

			cplex_idxs.push_back(div_cplex_id);
			coeffs.push_back(1);

			cplex_idxs.push_back(cplex_id(node_map_[n],1));
			coeffs.push_back(-1);

			if (with_appearance_) {
				cplex_idxs.push_back(cplex_id(app_node_map_[n],1));
				coeffs.push_back(-1);
			}
			// -1 <= D_i[1] - X_i[1] - App_i[1] <= 0
			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
					cplex_idxs.end(), coeffs.begin(), -1, 0);
			LOG(logDEBUG3) << "SingleTimestepTraxelConservation::add_constraints:" <<
					" D_i=1 => X_i + App_i =1 added for " << "n = " << node_map_[n] << ", d = " << div_node_map_[n];


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
					dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
							cplex_idxs.end(), coeffs.begin(), 0, 1);
					LOG(logDEBUG3) << "SingleTimestepTraxelConservation::add_constraints:" <<
							" D_i=1 => Y_i[nu]=0 added for " << "d = " << div_node_map_[n] << ", y = " << arc_map_[a] << ", nu = " << nu;

				}

				cplex_idxs2.push_back(cplex_id(arc_map_[a], 1));
				coeffs2.push_back(-1);
			}

			// -m <= 2 * D_i[1] - sum_j (Y_ij[1]) <= 0
			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs2.begin(),
					cplex_idxs2.end(), coeffs2.begin(), -int(max_number_objects_), 0);
			LOG(logDEBUG3) << "SingleTimestepTraxelConservation::add_constraints:" <<
					" D_i = 1 => sum_k(Y_ik) = 2 added for " << "d = " << div_node_map_[n];
		}


		////
		//// incoming transitions
		////
		// couple transitions: sum_k(Y_kj) = X_j + Dis_j
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
			for (size_t nu = 1; nu <= max_number_objects_; ++nu) {
				cplex_idxs.push_back(cplex_id(node_map_[n],nu));
				coeffs.push_back(-nu);
			}

			if (with_disappearance_) {
				for (size_t nu = 1; nu <= max_number_objects_; ++nu) {
					cplex_idxs.push_back(cplex_id(dis_node_map_[n],nu));
					coeffs.push_back(-nu);
				}
			}

			// 0 <= sum_nu [ nu * sum_i (Y_ij[nu] ) ] - sum_nu ( nu * X_j[nu] ) - sum_nu ( nu * Dis_j[nu] ) <= 0
			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
					cplex_idxs.end(), coeffs.begin(), 0, 0);
			LOG(logDEBUG3) << "SingleTimestepTraxelConservation::add_constraints:" <<
					" sum_k(Y_kj) = X_j + Dis_j added for " << "n = " << node_map_[n];
		}


		////
		//// disappearance/appearance coupling
		////
		if (with_appearance_ || with_disappearance_) {
			cplex_idxs.clear();
			coeffs.clear();

			for (size_t nu = 1; nu <= max_number_objects_; ++nu) {
				cplex_idxs.push_back(cplex_id(node_map_[n],nu));
				coeffs.push_back(1);
			}

			if (with_appearance_) {
				for (size_t nu = 1; nu <= max_number_objects_; ++nu) {
					cplex_idxs.push_back(cplex_id(app_node_map_[n],nu));
					coeffs.push_back(1);
				}
			}

			if (with_disappearance_) {
				for (size_t nu = 1; nu <= max_number_objects_; ++nu) {
					cplex_idxs.push_back(cplex_id(dis_node_map_[n],nu));
					coeffs.push_back(1);
				}
			}

			// 0 <= sum_nu ( X_i[nu] ) + sum_nu ( App_i[nu] ) + sum_nu ( Dis_i[nu] ) <= 1
			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
					cplex_idxs.end(), coeffs.begin(), 0, 1);
			LOG(logDEBUG3) << "SingleTimestepTraxelConservation::add_constraints:" <<
					" X_i>0 v` App_i>0 v` Dis_i>0 added for " << "n = " << node_map_[n];
		}
	}


//	if (with_appearance_) {
//		// set appearances in last time step zero
//		typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
//		node_timestep_map_t& node_timestep_map = g.get(node_timestep());
//		property_map<node_traxel, HypothesesGraph::base_graph>::type& node_traxel_map = g.get(node_traxel());
//		int t = g.latest_timestep();
//		vector<int> cplex_idxs;
//		vector<size_t> coeffs;
//		for(node_timestep_map_t::ItemIt n(node_timestep_map, t); n!=lemon::INVALID; ++n) {
//			assert(node_traxel_map[n].Timestep == g.latest_timestep());
//			coeffs.push_back(1);
//			cplex_idxs.push_back(cplex_id(app_node_map_[n],0));
//			// 1 <= App_i[0]<= 1
//			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
//						cplex_idxs.end(), coeffs.begin(), 1, 1);
//			LOG(logDEBUG3) << "SingleTimestepTraxelConservation::add_constraints:" <<
//					" 0 <= App_i <= 0 added for " << "n = " << node_map_[n]  << ", node.Id = " << node_traxel_map[n].Id;
//		}
//	}

//	if (with_disappearance_) {
//		// set disappearances in last time step zero
//		typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
//		node_timestep_map_t& node_timestep_map = g.get(node_timestep());
//		property_map<node_traxel, HypothesesGraph::base_graph>::type& node_traxel_map = g.get(node_traxel());
//		int t = g.latest_timestep();
//		vector<int> cplex_idxs;
//		vector<size_t> coeffs;
//		for(node_timestep_map_t::ItemIt n(node_timestep_map, t); n!=lemon::INVALID; ++n) {
//			assert(node_traxel_map[n].Timestep == g.latest_timestep());
//			coeffs.push_back(1);
//			cplex_idxs.push_back(cplex_id(dis_node_map_[n],1));
//			// 1 <= Dis_i[0]<= 1
//			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
//						cplex_idxs.end(), coeffs.begin(), 1, 1);
//			LOG(logDEBUG3) << "SingleTimestepTraxelConservation::add_constraints:" <<
//					" 0 <= Dis_i <= 0 added for " << "n = " << node_map_[n] << ", node.Id = " << node_traxel_map[n].Id;
//		}
//	}
}



void SingleTimestepTraxelConservation::fix_detections(const HypothesesGraph& g) {
	typedef opengm::LPCplex<OpengmModel::ogmGraphicalModel,OpengmModel::ogmAccumulator> cplex;
	property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
	property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = g.get(node_tracklet());
	property_map<tracklet_intern_dist, HypothesesGraph::base_graph>::type& tracklet_intern_dist_map = g.get(tracklet_intern_dist());

	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		vector<size_t> cplex_idxs;
		vector<int> coeffs;
		vector<double> energies;
		for (size_t state = 0; state <= max_number_objects_; ++state) {
			if (with_tracklets_) {
				double e = 0;
				// add all detection factors of the internal nodes
				for (std::vector<Traxel>::const_iterator trax_it = tracklet_map[n].begin(); trax_it != tracklet_map[n].end(); ++trax_it){
					e += detection_(*trax_it, state);
				}
				// add all transition factors of the internal arcs
				for (std::vector<double>::const_iterator intern_dist_it = tracklet_intern_dist_map[n].begin();
						intern_dist_it != tracklet_intern_dist_map[n].end(); ++intern_dist_it) {
					e += transition_(get_transition_prob(*intern_dist_it, state));
				}
				energies.push_back(e);
			} else {
				energies.push_back((double) detection_(traxel_map[n], state));
			}
		}
		size_t max_state = std::min_element(energies.begin(), energies.end()) - energies.begin();

		// if it is most probable that it is a misdetection (max_state == 0), then fix all X_i, App_i, Dis_i to 0
		// otherwise, add constraint that one of those must be greater than 0
		// (Detections cannot be fixed to their max_state, example: 1 --- 2 --- 1 would be invalid)

		if (max_state == 0) {
//			cplex_idxs.push_back(cplex_id(node_map_[n], 0));
//			coeffs.push_back(1);
//			// X_i[0] = 1
//			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
//					cplex_idxs.end(), coeffs.begin(), 1, 1);
//			LOG(logDEBUG3) << "SingleTimestepTraxelConservation::fix_detections:" <<
//								" X_i set to 0 for n = " << node_map_[n];
//			if (with_appearance_) {
//				cplex_idxs.clear();
//				coeffs.clear();
//				cplex_idxs.push_back(cplex_id(app_node_map_[n],0));
//				coeffs.push_back(1);
//				// App_i[0] = 1
//				dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
//						cplex_idxs.end(), coeffs.begin(), 1, 1);
//				LOG(logDEBUG3) << "SingleTimestepTraxelConservation::fix_detections:" <<
//									" App_i set to 0 for n = " << node_map_[n];
//			}
//			if (with_disappearance_) {
//				cplex_idxs.clear();
//				coeffs.clear();
//				cplex_idxs.push_back(cplex_id(dis_node_map_[n],0));
//				coeffs.push_back(1);
//				// Dis_i[0] = 1
//				dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
//						cplex_idxs.end(), coeffs.begin(), 1, 1);
//				LOG(logDEBUG3) << "SingleTimestepTraxelConservation::fix_detections:" <<
//									" Dis_i set to 0 for n = " << node_map_[n];
//			}

//		}
		} else {  // max_state > 0
			// one of X_i, App_i or Dis_i must be greater than 0
			cplex_idxs.push_back(cplex_id(node_map_[n], 0));
			coeffs.push_back(1);
			size_t num_vars = 1;

			if (with_appearance_) {
				cplex_idxs.push_back(cplex_id(app_node_map_[n],0));
				coeffs.push_back(1);
				++num_vars;
			}
			if (with_disappearance_) {
				cplex_idxs.push_back(cplex_id(dis_node_map_[n],0));
				coeffs.push_back(1);
				++num_vars;
			}
			// 2 <= 1*X_i[0] + 1*App_i[0] + 1*Dis_i[0] <= 2
			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
					cplex_idxs.end(), coeffs.begin(), num_vars-1, num_vars-1);
			LOG(logDEBUG3) << "SingleTimestepTraxelConservation::fix_detections:" <<
								" X_i != 0 added for " << "n = " << node_map_[n];
		}
	}
}


} /* namespace Tracking */
