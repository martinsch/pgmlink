#include <algorithm>
#include <cassert>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>
#include <boost/scoped_ptr.hpp>
#include <string.h>
#include <memory.h>

#ifdef WITH_GUROBI
#include <opengm/inference/lpgurobi.hxx>
#else
#include <opengm/inference/lpcplex.hxx>
#endif

#include <opengm/datastructures/marray/marray.hxx>

#include "pgmlink/pgm.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_pgm.h"
#include "pgmlink/pgm_chaingraph.h"
#include "pgmlink/traxels.h"

//#include <ostream>

using namespace std;

namespace pgmlink {
  ////
  //// class Chaingraph
  ////

  Chaingraph::~Chaingraph() {
    if(optimizer_ != NULL) {
	delete optimizer_;
	optimizer_ = NULL;
    }

    delete builder_;
    builder_ = NULL;
  }

double Chaingraph::forbidden_cost() const {
  return builder_->forbidden_cost();
}

bool Chaingraph::with_constraints() const {
  return with_constraints_;
}

void Chaingraph::formulate( const HypothesesGraph& hypotheses ) {
    LOG(logDEBUG) << "Chaingraph::formulate: entered";
    reset();

    // build the model
    linking_model_ = boost::shared_ptr<pgm::chaingraph::Model>(builder_->build(hypotheses));

    // refine the model with hard constraints
    pgm::OpengmLPCplex::Parameter param;
    param.verbose_ = true;
    param.integerConstraint_ = true;
    param.epGap_ = ep_gap_;
    param.timeLimit_ = cplex_timeout_;
    LOG(logDEBUG) << "Chaingraph::formulate ep_gap = " << param.epGap_;
    pgm::OpengmLPCplex* cplex = new pgm::OpengmLPCplex(*(linking_model_->opengm_model), param);
    optimizer_ = cplex; // opengm::Inference optimizer_

    if (with_constraints_) {
      LOG(logDEBUG) << "Chaingraph::formulate: add_constraints";
      builder_->add_hard_constraints( *linking_model_ , hypotheses, *cplex );
    }
    
    if (fixed_detections_) {
      LOG(logDEBUG) << "Chaingraph::formulate: fix_detections";
      builder_->fix_detections( *linking_model_, hypotheses, *cplex );
    }
}



void Chaingraph::infer() {
    opengm::InferenceTermination status = optimizer_->infer();
    if(status != opengm::NORMAL) {
        throw std::runtime_error("GraphicalModel::infer(): optimizer terminated unnormally");
    }
}


void Chaingraph::conclude( HypothesesGraph& g ) {
    // extract solution from optimizer
  vector<pgm::OpengmLPCplex::LabelType> solution;
    opengm::InferenceTermination status = optimizer_->arg(solution);
    if(status != opengm::NORMAL) {
	throw runtime_error("GraphicalModel::infer(): solution extraction terminated unnormally");
    }

    // add 'active' properties to graph
    g.add(node_active()).add(arc_active());
    property_map<node_active, HypothesesGraph::base_graph>::type& active_nodes = g.get(node_active());
    property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = g.get(arc_active());

    // write state after inference into 'active'-property maps
    if(builder_->has_detection_vars()) {
      for(pgm::chaingraph::Model::node_var_map::const_iterator it = linking_model_->var_of_node().begin(); it != linking_model_->var_of_node().end(); ++it) {
	bool state = false;
	if(solution[it->second] == 1) state = true;
	active_nodes.set(it->first, state);
      }
    } else {
      // if we have no detection vars we set all nodes to active by default
      active_nodes.setAll(true); 
    }
    for(pgm::chaingraph::Model::arc_var_map::const_iterator it = linking_model_->var_of_arc().begin(); it != linking_model_->var_of_arc().end(); ++it) {
	bool state = false;
	if(solution[it->second] == 1) state = true;
	active_arcs.set(it->first, state);
    }
}

  const pgm::OpengmModel* Chaingraph::get_graphical_model() const {
    return linking_model_->opengm_model.get();
  }

  const Chaingraph::node_var_map& Chaingraph::get_node_map() const {
    return linking_model_->var_of_node();
  }

  const Chaingraph::arc_var_map& Chaingraph::get_arc_map() const {
    return linking_model_->var_of_arc();
  }

void Chaingraph::reset() {
    if(optimizer_ != NULL) {
	delete optimizer_;
	optimizer_ = NULL;
    }
}

} /* namespace pgmlink */ 
