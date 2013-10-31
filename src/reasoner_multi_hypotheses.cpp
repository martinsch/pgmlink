// stl
#include <vector>
#include <stdexcept>
#include <cassert>

// boost
#include <boost/shared_ptr.hpp>

// opengm
#include <opengm/inference/lpcplex.hxx>
#include <opengm/datastructures/marray/marray.hxx>

// pgmlink
#include "pgmlink/log.h"
#include "pgmlink/pgm.h"
#include "pgmlink/traxels.h"
#include "pgmlink/multi_hypotheses_graph.h"
#include "pgmlink/reasoner_multi_hypotheses.h"
#include "pgmlink/pgm_multi_hypotheses.h"



namespace pgmlink {

////
//// class MultiHypotheses
////

MultiHypotheses::~MultiHypotheses() {
  reset();
}


// double MultiHypotheses::forbidden() const {
//   return builder_->forbidden_cost();
// }


bool MultiHypotheses::with_constraints() const {
  return with_constraints_;
}


void MultiHypotheses::formulate( const MultiHypothesesGraph& hypotheses ) {
  LOG(logDEBUG) << "MultiHypotheses::formlate: entered";
  reset();

  // build the model
  linking_model_ = boost::shared_ptr<pgm::multihypotheses::Model>(builder_->build(hypotheses));

  // refine the model with hard constraints
  pgm::OpengmLPCplex::Parameter param;
  param.verbose_ = true;
  param.integerConstraint_ = true;
  param.epGap_ = ep_gap_;
  param.timeLimit_ = cplex_timeout_;
  LOG(logDEBUG) << "MultiHypotheses::formulate ep_gap = " << param.epGap_;
  pgm::OpengmLPCplex* cplex = new pgm::OpengmLPCplex(*(linking_model_->opengm_model), param);
  LOG(logDEBUG) << "MultiHypotheses:formulate created model";
  optimizer_ = cplex;

  if (with_constraints_) {
    LOG(logDEBUG) << "MultiHypotheses::formulate: add hard constraints";
    builder_->add_hard_constraints( *linking_model_, hypotheses, *cplex );
  }
  if (builder_->has_hierarchical_counting_factor() || builder_->has_counting_incoming_factor() ) {
	  builder_->add_count_hard_constraints( *linking_model_, hypotheses, *cplex );
  }
}


double MultiHypotheses::infer() {
  opengm::InferenceTermination status = optimizer_->infer();
  if (status != opengm::NORMAL) {
    throw std::runtime_error("GraphicalModel::infer(): optimizer terminated unnormaly");
  }
  return optimizer_->value();
}


void MultiHypotheses::conclude( MultiHypothesesGraph& g ) {
  // extract solution from optimizer
  std::vector<pgm::OpengmLPCplex::LabelType> solution;
  opengm::InferenceTermination status = optimizer_->arg(solution);
  if (status != opengm::NORMAL) {
    throw std::runtime_error("GraphicalModel::infer(): solution extraction terminated unnormally");
  }

  MultiHypothesesGraph::ContainedRegionsMap& regions = g.get(node_regions_in_component());
  MultiHypothesesGraph::node_timestep_map& timesteps = g.get(node_timestep());
  for (MultiHypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
    std::vector<Traxel>& traxels = regions.get_value(n);
    if (builder_->timestep_range_specified() &&
        (timesteps[n] < builder_->first_timestep() || timesteps[n] > builder_->last_timestep())) {
      continue;
    }
    for (std::vector<Traxel>::iterator t = traxels.begin(); t != traxels.end(); ++t) {

      // set traxels active
      if (builder_->has_detection_vars()) {
        float state = 0.;
        trax_var_map::const_iterator it = linking_model_->var_of_trax().find(*t);
        assert(it != linking_model_->var_of_trax().end());
        if (solution[it->second] == 1) {
          state = 1.;
        }
        t->features["active"] = feature_array(1, state);
        
      } else {
        t->features["active"] = feature_array(1, 1.);
      }

      // add edges lists
      for (MultiHypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
        if (builder_->timestep_range_specified() &&
            (timesteps[g.target(a)] > builder_->last_timestep())) {
          continue;
        }
        const std::vector<Traxel>& neighbors = regions[g.target(a)];
        for (std::vector<Traxel>::const_iterator neighbor = neighbors.begin();
             neighbor != neighbors.end();
             ++neighbor) {
          arc_var_map::const_iterator it = linking_model_->var_of_arc().find(pgm::multihypotheses::Model::TraxelArc(*t, *neighbor));
          // assert(it != linking_model_->var_of_arc().end());
          if (it != linking_model_->var_of_arc().end() && solution[it->second] == 1) {
            t->features["outgoing"].push_back(neighbor->Id);
          }
        }
      }
      assert(t->features["outgoing"].size() <= 2);

      // add parent if any
      for (MultiHypothesesGraph::InArcIt a(g, n); a != lemon::INVALID; ++a) {
        if (builder_->timestep_range_specified() &&
            (timesteps[g.source(a)] < builder_->first_timestep())) {
          continue;
        }
        const std::vector<Traxel>& neighbors = regions[g.source(a)];
        for (std::vector<Traxel>::const_iterator neighbor = neighbors.begin();
             neighbor != neighbors.end();
             ++neighbor) {
          arc_var_map::const_iterator it = linking_model_->var_of_arc().find(pgm::multihypotheses::Model::TraxelArc(*neighbor, *t));
          // assert(it != linking_model_->var_of_arc().end());
          if (it != linking_model_->var_of_arc().end() && solution[it->second] == 1) {
            t->features["parent"].push_back(neighbor->Id);
          }
        }
      }
      assert(t->features["parent"].size() <= 1);
    }
  }
}


const pgm::OpengmModel* MultiHypotheses::get_graphical_model() const {
  return linking_model_->opengm_model.get();
}


const MultiHypotheses::trax_var_map& MultiHypotheses::get_trax_map() const {
  return linking_model_->var_of_trax();
}


const MultiHypotheses::arc_var_map& MultiHypotheses::get_arc_map() const {
  return linking_model_->var_of_arc();
}

void MultiHypotheses::reset() {
  if (optimizer_ != NULL) {
    delete optimizer_;
    optimizer_ = NULL;
  }
}


}


