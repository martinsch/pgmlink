// stl
#include <vector>
#include <stdexcept>
#include <cassert>

// boost
#include <boost/shared_ptr.hpp>

// opengm
#include <memory.h>
#include <string.h>
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

  // add 'active' properties to graph
  g.add(node_active()).add(arc_active());
  property_map<node_active, HypothesesGraph::base_graph>::type& active_nodes = g.get(node_active());
  property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = g.get(arc_active());
  MultiHypothesesGraph::TraxelMap& traxel_map = g.get(node_traxel());

  // write state after inference into 'active'-property maps
  if(builder_->has_detection_vars()) {
    for(pgm::multihypotheses::Model::node_var_map::const_iterator it = linking_model_->var_of_node().begin(); it != linking_model_->var_of_node().end(); ++it) {
      bool state = false;
      if(solution[it->second] == 1) state = true;
      active_nodes.set(it->first, state);
      LOG(logDEBUG4) << "GraphicalModel::conclude detection " << traxel_map[it->first] << " has state " << state;
    }
  } else {
    // if we have no detection vars we set all nodes to active by default
    active_nodes.setAll(true); 
  }
  for(pgm::multihypotheses::Model::arc_var_map::const_iterator it = linking_model_->var_of_arc().begin(); it != linking_model_->var_of_arc().end(); ++it) {
    bool state = false;
    if(solution[it->second] == 1) state = true;
    active_arcs.set(it->first, state);
  }
}


const pgm::OpengmModel* MultiHypotheses::get_graphical_model() const {
  return linking_model_->opengm_model.get();
}


const MultiHypotheses::node_var_map& MultiHypotheses::get_trax_map() const {
  return linking_model_->var_of_node();
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


