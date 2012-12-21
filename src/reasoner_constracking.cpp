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

	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_detection_nodes";
	add_detection_nodes(hypotheses);
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_transition_nodes";
	add_transition_nodes(hypotheses);
	if (with_appearance_) {
		LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_appearance_nodes";
		add_appearance_nodes(hypotheses);
	}
	if (with_disappearance_) {
		LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_disappearance_nodes";
		add_disappearance_nodes(hypotheses);
	}
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_division_nodes";
	add_division_nodes(hypotheses);

	OpengmModel::ogmGraphicalModel* model = pgm_->Model();

	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_finite_factors";
	add_finite_factors(hypotheses);

	typedef opengm::LPCplex<OpengmModel::ogmGraphicalModel,OpengmModel::ogmAccumulator> cplex_optimizer;
	cplex_optimizer::Parameter param;
	param.verbose_ = true;
	param.integerConstraint_ = true;
	param.epGap_ = ep_gap_;
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate ep_gap = " << param.epGap_;

	optimizer_ = new cplex_optimizer(*model, param);

	if (with_constraints_) {
		LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_constraints";
		add_constraints(hypotheses);
	}

	if (fixed_detections_) {
		assert(with_appearance_ || with_disappearance_);
		LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: fix_detections";
		fix_detections(hypotheses);
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

	vector<size_t> count_objects(max_number_objects_+1,0);

	// write state after inference into 'active'-property maps
	for (std::map<HypothesesGraph::Node, size_t>::const_iterator it =
			node_map_.begin(); it != node_map_.end(); ++it) {
		++count_objects[solution[it->second]];
		active_nodes.set(it->first, solution[it->second]);
	}
	// the node is also active if its appearance node is active
	for (std::map<HypothesesGraph::Node, size_t>::const_iterator it =
			app_node_map_.begin(); it != app_node_map_.end(); ++it) {
		++count_objects[solution[it->second]];
		if (solution[it->second] > 0) {
			active_nodes.set(it->first, solution[it->second]);
		}
	}
	// the node is also active if its disappearance node is active
	for (std::map<HypothesesGraph::Node, size_t>::const_iterator it =
			dis_node_map_.begin(); it != dis_node_map_.end(); ++it) {
		++count_objects[solution[it->second]];
		if (solution[it->second] > 0) {
			active_nodes.set(it->first, solution[it->second]);
		}
	}

	for (std::map<HypothesesGraph::Arc, size_t>::const_iterator it =
			arc_map_.begin(); it != arc_map_.end(); ++it) {
		bool state = false;
		if (solution[it->second] >= 1)
			state = true;
		active_arcs.set(it->first, state);
	}
	// initialize division node map
	for (std::map<HypothesesGraph::Node, size_t>::const_iterator it =
			node_map_.begin(); it != node_map_.end(); ++it) {
		division_nodes.set(it->first, false);
	}
	for (std::map<HypothesesGraph::Node, size_t>::const_iterator it =
			div_node_map_.begin(); it != div_node_map_.end(); ++it) {
		bool state = false;
		if (solution[it->second] >=1)
			state = true;
		division_nodes.set(it->first, state);
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
			double energy = detection_(traxel_map[n], state);
			LOG(logDEBUG2) << "SingleTimestepTraxelConservation::add_finite_factors: detection[" << state <<
							"] = " << energy;
			for (size_t var_idx = 0; var_idx < num_vars; ++var_idx) {
				coords[var_idx] = state;
				table.set_value(coords, energy);
				coords[var_idx] = 0;
				LOG(logDEBUG3) << "SingleTimestepTraxelConservation::add_finite_factors: var_idx " << var_idx <<
									"] = " << energy;
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
			double energy = division_(traxel_map[n], state);
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
							" Y_ij <= X_ij + App_i added for "
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

			// 0 <= sum_nu [ nu * sum_i (Y_ij[nu] ) ] - sum_nu ( nu * X_j[nu] ) <= 0
			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
					cplex_idxs.end(), coeffs.begin(), 0, 0);
			LOG(logDEBUG3) << "SingleTimestepTraxelConservation::add_constraints:" <<
					" sum_k(Y_kj) = X_j added for "	<< "n = " << node_map_[n];
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
}



void SingleTimestepTraxelConservation::fix_detections(const HypothesesGraph& g) {
	typedef opengm::LPCplex<OpengmModel::ogmGraphicalModel,OpengmModel::ogmAccumulator> cplex;
	property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());

	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		vector<size_t> cplex_idxs;
		vector<int> coeffs;
		vector<double> energies;
		for (size_t state = 0; state <= max_number_objects_; ++state) {
			energies.push_back((double) detection_(traxel_map[n], state));
		}
		size_t max_state = std::min_element(energies.begin(), energies.end()) - energies.begin();
		cplex_idxs.push_back(cplex_id(node_map_[n], max_state));
		coeffs.push_back(1);
		if (with_appearance_) {
			cplex_idxs.push_back(cplex_id(app_node_map_[n],max_state));
			coeffs.push_back(1);
		}
		if (with_disappearance_) {
			cplex_idxs.push_back(cplex_id(dis_node_map_[n],max_state));
			coeffs.push_back(1);
		}
		// 1 <= 1*X_i[max_state] + 1*App_i[max_state] + 1*Dis_i[max_state] <= val
		dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
				cplex_idxs.end(), coeffs.begin(), 1, 1);
		LOG(logDEBUG3) << "SingleTimestepTraxelConservation::fix_detections:" <<
							" fix_Detection to " << max_state << " added for " << "n = " << node_map_[n];
	}
}


} /* namespace Tracking */
