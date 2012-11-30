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
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_division_nodes";
	add_division_nodes(hypotheses);

//	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_finite_factors";
//	add_finite_factors(hypotheses);

	OpengmModel::ogmGraphicalModel* model = pgm_->Model();

	typedef opengm::LPCplex<OpengmModel::ogmGraphicalModel,OpengmModel::ogmAccumulator> cplex_optimizer;
	cplex_optimizer::Parameter param;
	param.verbose_ = true;
	param.integerConstraint_ = true;
	param.epGap_ = ep_gap_;
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate ep_gap = " << param.epGap_;

	optimizer_ = new cplex_optimizer(*model, param);

	add_finite_factors(hypotheses);

	if (with_constraints_) {
		LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: add_constraints";
		add_constraints(hypotheses);
	}

	if (fixed_detections_) {
		LOG(logDEBUG) << "SingleTimestepTraxelConservation::formulate: fix_detections";
		fix_detections(hypotheses, 1);
	}
}

void SingleTimestepTraxelConservation::infer() {
	opengm::InferenceTermination status = optimizer_->infer();
	if (status != opengm::NORMAL) {
		throw std::runtime_error(
				"GraphicalModel::infer(): optimizer terminated unnormally");
	}
}

void SingleTimestepTraxelConservation::conclude(HypothesesGraph& g) {
	// extract solution from optimizer
	vector<OpengmModel::ogmInference::LabelType> solution;
	opengm::InferenceTermination status = optimizer_->arg(solution);
	if (status != opengm::NORMAL) {
		throw runtime_error(
				"GraphicalModel::infer(): solution extraction terminated unnormally");
	}

	// add 'active' properties to graph
	g.add(node_active()).add(arc_active());
	property_map<node_active, HypothesesGraph::base_graph>::type& active_nodes =
			g.get(node_active());
	property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs =
			g.get(arc_active());

	// write state after inference into 'active'-property maps
	for (std::map<HypothesesGraph::Node, size_t>::const_iterator it =
			node_map_.begin(); it != node_map_.end(); ++it) {
		bool state = false;
		if (solution[it->second] == 1)
			state = true;
		active_nodes.set(it->first, state);
	}
	for (std::map<HypothesesGraph::Arc, size_t>::const_iterator it =
			arc_map_.begin(); it != arc_map_.end(); ++it) {
		bool state = false;
		if (solution[it->second] == 1)
			state = true;
		active_arcs.set(it->first, state);
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
}

void SingleTimestepTraxelConservation::add_detection_nodes(const HypothesesGraph& g) {
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		pgm_->Model()->addVariable(max_number_objects_+1);
		node_map_[n] = pgm_->Model()->numberOfVariables() - 1;
	}
}

void SingleTimestepTraxelConservation::add_transition_nodes(const HypothesesGraph& g) {
	for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
		pgm_->Model()->addVariable(max_number_objects_+1);
		arc_map_[a] = pgm_->Model()->numberOfVariables() - 1;
	}
}

void SingleTimestepTraxelConservation::add_division_nodes(const HypothesesGraph& g) {
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		size_t number_of_outarcs = 0;
		for(HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
			++number_of_outarcs;
		}
		if (number_of_outarcs > 1) {
			pgm_->Model()->addVariable(2);
			div_node_map_[n] = pgm_->Model()->numberOfVariables() - 1;
		}
	}
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
		size_t vi[] = { node_map_[n] };
		vector<size_t> coords(1, 0); // number of variables
		OpengmExplicitFactor<double> table(vi, vi + 1);
		for (size_t state = 0; state <= max_number_objects_; ++state) {
			LOG(logDEBUG3) << "1";
			double energy = detection_(traxel_map[n], state);
			LOG(logDEBUG2) << "SingleTimestepTraxelConservation::add_finite_factors: detection[" << state <<
					"] = " << energy;
			coords[0] = state;
			LOG(logDEBUG3) << "2";
			table.set_value(coords, energy);
			LOG(logDEBUG3) << "3";
			coords[0] = 0;
			LOG(logDEBUG3) << "4";
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
		OpengmExplicitFactor<double> table(vi, vi + 1);
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
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::add_finite_factors: add detection factors";
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		if (div_node_map_.count(n) == 0) {
			continue;
		}
		size_t vi[] = { div_node_map_[n] };
		vector<size_t> coords(1, 0); // number of variables
		OpengmExplicitFactor<double> table(vi, vi + 1);
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

namespace {
inline size_t cplex_id(size_t opengm_id) {
	return 2 * opengm_id + 1;
}
}

void SingleTimestepTraxelConservation::add_constraints(const HypothesesGraph& g) {
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::add_constraints: entered";
	typedef opengm::LPCplex<OpengmModel::ogmGraphicalModel,OpengmModel::ogmAccumulator> cplex;

	////
	//// outgoing transitions
	////
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::add_constraints: outgoing transitions";
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		vector<size_t> cplex_idxs;
		vector<int> coeffs;

		// couple detection and transitions: Y_ij <= X_i
		for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
			cplex_idxs.clear();
			coeffs.clear();

			cplex_idxs.push_back(cplex_id(node_map_[n]));
			cplex_idxs.push_back(cplex_id(arc_map_[a]));
			coeffs.push_back(1);
			coeffs.push_back(-1);
			// 0 <= X_i - Y_i <= m
			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
				cplex_idxs.end(), coeffs.begin(), 0, max_number_objects_);
		}

		// couple transitions: sum(Y_ij) = D_i + X_i
		cplex_idxs.clear();
		coeffs.clear();
		for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
			cplex_idxs.push_back(cplex_id(arc_map_[a]));
			coeffs.push_back(1);
		}
		int div_cplex_id = -1;
		if (div_node_map_.count(n) > 0) {
			div_cplex_id = cplex_id(div_node_map_[n]);
			cplex_idxs.push_back(div_cplex_id);
			coeffs.push_back(-1);
		}
		cplex_idxs.push_back(cplex_id(node_map_[n]));
		coeffs.push_back(-1);

		// 0 <= sum_j(Y_ij) - X_i - D_i <= 0
		dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
				cplex_idxs.end(), coeffs.begin(), 0, 0);


		if (div_cplex_id != -1) {
			// couple detection and division: D_i = 1 => X_i = 1
			cplex_idxs.clear();
			coeffs.clear(); // coeffs: D_i = 1 => X_i >= 1
			vector<int> coeffs2; // coeffs2: D_1 = 1 => X_i <=1

			cplex_idxs.push_back(div_cplex_id);
			coeffs.push_back(1);
			coeffs2.push_back(max_number_objects_);

			cplex_idxs.push_back(cplex_id(node_map_[n]));
			coeffs.push_back(-1);
			coeffs2.push_back(1);

			// -m <= D_i - X_i <= 0
			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
					cplex_idxs.end(), coeffs.begin(), -max_number_objects_, 0);

			// 0 <= m * D_i + X_i <= m+1
			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
					cplex_idxs.end(), coeffs2.begin(), 0, max_number_objects_+1);


			// couple divsion and transition: D_1 = 1 => sum_k(Y_ik) = 2
			cplex_idxs.clear();
			coeffs.clear();
			cplex_idxs.push_back(div_cplex_id);
			coeffs.push_back(2);

			for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
				cplex_idxs.push_back(arc_map_[a]);
				coeffs.push_back(-1);
			}

			// -m <= 2 * D_i - sum_j(Y_ij) <= 0
			dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
					cplex_idxs.end(), coeffs.begin(), - max_number_objects_, 0);

		}


	}

	////
	//// incoming transitions
	////
	LOG(logDEBUG) << "SingleTimestepTraxelConservation::add_constraints: incoming transitions";
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		// couple transitions: sum_k(Y_kj) = X_j
		vector<size_t> cplex_idxs;
		vector<int> coeffs;
		for (HypothesesGraph::InArcIt a(g, n); a != lemon::INVALID; ++a) {
			cplex_idxs.push_back(cplex_id(arc_map_[a]));
			coeffs.push_back(1);
		}
		cplex_idxs.push_back(cplex_id(node_map_[n]));
		coeffs.push_back(-1);

		// 0 <= sum_k(Y_kj) - X_j <= 0
		dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
				cplex_idxs.end(), coeffs.begin(), 0, 0);
	}
}



void SingleTimestepTraxelConservation::fix_detections(const HypothesesGraph& g,
		size_t val) {
	typedef opengm::LPCplex<OpengmModel::ogmGraphicalModel,OpengmModel::ogmAccumulator> cplex;
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
		vector<size_t> cplex_idxs;
		cplex_idxs.push_back(cplex_id(node_map_[n]));
		vector<int> coeffs;
		coeffs.push_back(1);
		// val <= 1*detection <= val
		dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(),
				cplex_idxs.end(), coeffs.begin(), val, val);
	}
}


} /* namespace Tracking */
