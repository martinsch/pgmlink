#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <utility>
#include <opengm/inference/lpcplex.hxx>
#include <opengm/datastructures/marray/marray.hxx>

#include "pgmlink/graphical_model.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_opengm.h"
#include "pgmlink/traxels.h"

//#include <ostream>

using namespace std;

namespace pgmlink {
  namespace pgm {
    ////
    //// class FactorEntry
    ////
    inline void FactorEntry::set( OpengmModel::ValueType v ) {
      if( entry_ == NULL ) {
	throw runtime_error("FactorEntry: is invalid");
      }
      *entry_ = v;
    }

    inline OpengmModel::ValueType FactorEntry::get() const {
      if( entry_ == NULL ) {
	throw runtime_error("FactorEntryPtr: is invalid");
      }
      return *entry_;
    } 
    
    ////
    //// class ChaingraphModelBuilder
    ////
    ChaingraphModelBuilder& ChaingraphModelBuilder::hypotheses(shared_ptr<const HypothesesGraph> g ) {
      if(g == NULL) {
	throw invalid_argument("ChaingraphModelBuilder::hypotheses(): null pointer");
      }
      hypotheses_ = g;
      return *this;
    }
    
    ChaingraphModelBuilder& ChaingraphModelBuilder::appearance( function<double (const Traxel&)> f ) {
      if(!f) {
	throw invalid_argument("ChaingraphModelBuilder::appearance(): empty function");
      }
      appearance_ = f;
      return *this;
    }

    ChaingraphModelBuilder& ChaingraphModelBuilder::disappearance( function<double (const Traxel&)> f ) {
      if(!f) {
	throw invalid_argument("ChaingraphModelBuilder::disappearance(): empty function");
      }
      disappearance_ = f;
      return *this;
    }

    ChaingraphModelBuilder& ChaingraphModelBuilder::move( function<double (const Traxel&,const Traxel&)> f ) {
      if(!f) {
	throw invalid_argument("ChaingraphModelBuilder::move(): empty function");
      }
      move_ = f;
      return *this;
    }

    ChaingraphModelBuilder& ChaingraphModelBuilder::with_detection_vars( function<double (const Traxel&)> detection,
									 function<double (const Traxel&)> non_detection) {
      if(!(detection && non_detection)) {
	throw invalid_argument("ChaingraphModelBuilder::with_detection_vars(): empty function");
      }

      with_detection_vars_ = true;
      detection_ = detection;
      non_detection_ = non_detection;
      return *this;
    }

    ChaingraphModelBuilder& ChaingraphModelBuilder::without_detection_vars() {
      with_detection_vars_ = false;
      detection_ = NULL;
      non_detection_ = NULL;
      return *this;
    }

    ChaingraphModelBuilder& ChaingraphModelBuilder::with_divisions( function<double (const Traxel&,const Traxel&,const Traxel&)> division) {
      if(!division) {
	throw invalid_argument("ChaingraphModelBuilder::division(): empty function");
      }
      with_divisions_ = true;
      division_ = division;
      return *this;
    }

    ChaingraphModelBuilder& ChaingraphModelBuilder::without_divisions() {
      with_divisions_ = false;
      division_ = NULL;
      return *this;
    }

    boost::shared_ptr<TrackingModel> ChaingraphModelBuilder::build() const {
      using boost::shared_ptr;
      using std::map;

      LOG(logDEBUG) << "ChaingraphModelBuilder::build: entered";
      if( !with_detection_vars_ ) {
	throw std::runtime_error("ChaingraphModelBuilder::build(): option without detection vars not yet implemented");
      }
      if( !with_divisions_ ) {
	throw std::runtime_error("ChaingraphModelBuilder::build(): option without divisions not yet implemented");
      }

      shared_ptr<TrackingModel> model( new TrackingModel() );
      
      if( with_detection_vars_ ) {
	add_detection_vars( *model );
      }
      add_assignment_vars( *model );

      if( with_detection_vars_ ) {
	for(HypothesesGraph::NodeIt n(*hypotheses_); n!=lemon::INVALID; ++n) {
	  add_detection_factor( *model, n );
	}
      }

      for(HypothesesGraph::NodeIt n(*hypotheses_); n!=lemon::INVALID; ++n) {
	add_outgoing_factor( *model, n );
 	add_incoming_factor( *model, n );
      }

      return model;
    }
    namespace {
      inline size_t cplex_id(size_t opengm_id) {
	return 2*opengm_id + 1;
      }
    }

    void ChaingraphModelBuilder::add_hard_constraints( const TrackingModel& m, const HypothesesGraph& hypotheses, OpengmLPCplex& cplex ) {
      LOG(logDEBUG) << "Chaingraph::add_constraints: entered";
      ////
      //// outgoing transitions
      ////
      LOG(logDEBUG) << "Chaingraph::add_constraints: outgoing transitions";
      for(HypothesesGraph::NodeIt n(hypotheses); n!=lemon::INVALID; ++n) {
	// couple detection and transitions
	for(HypothesesGraph::OutArcIt a(hypotheses, n); a!=lemon::INVALID; ++a) {
	  couple(m, n, a, cplex);
	}
	// couple transitions
	vector<size_t> cplex_idxs;
	for(HypothesesGraph::OutArcIt a(hypotheses, n); a!=lemon::INVALID; ++a) {
	  cplex_idxs.push_back(cplex_id(m.arc_var.find(a)->second));
	}
	vector<int> coeffs(cplex_idxs.size(), 1);
	// 0 <= 1*transition + ... + 1*transition <= 2
	cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 2);
      }

      ////
      //// incoming transitions
      ////
      LOG(logDEBUG) << "Chaingraph::add_constraints: incoming transitions";
      for(HypothesesGraph::NodeIt n(hypotheses); n!=lemon::INVALID; ++n) {
	// couple detection and transitions
	for(HypothesesGraph::InArcIt a(hypotheses, n); a!=lemon::INVALID; ++a) {
	  couple(m, n, a, cplex);
	}
	    
	// couple transitions
	vector<size_t> cplex_idxs;
	for(HypothesesGraph::InArcIt a(hypotheses, n); a!=lemon::INVALID; ++a) {
	  cplex_idxs.push_back(cplex_id(m.arc_var.find(a)->second));
	}
	vector<int> coeffs(cplex_idxs.size(), 1);
	// 0 <= 1*transition + ... + 1*transition <= 1
	cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 1);
      }
    }

    inline void ChaingraphModelBuilder::add_detection_vars( TrackingModel& m ) const {
      for(HypothesesGraph::NodeIt n(*hypotheses_); n!=lemon::INVALID; ++n) {
	m.opengm_model->addVariable(2);
	m.node_var[n] = m.opengm_model->numberOfVariables() - 1; 
      }
    }

    inline void ChaingraphModelBuilder::add_assignment_vars( TrackingModel& m ) const {
      for(HypothesesGraph::ArcIt a(*hypotheses_); a!=lemon::INVALID; ++a) {
	m.opengm_model->addVariable(2);
	m.arc_var[a] = m.opengm_model->numberOfVariables() - 1; 
      }
    }

    void ChaingraphModelBuilder::add_detection_factor( TrackingModel& m, const HypothesesGraph::Node& n) const {
      property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = hypotheses_->get(node_traxel());
      size_t vi[] = {m.node_var[n]};
      vector<size_t> coords(1,0);
      OpengmExplicitFactor<double> table( vi, vi+1 );

      coords[0] = 0;
      table.set_value( coords, non_detection_(traxel_map[n]) ); 

      coords[0] = 1;
      table.set_value( coords, detection_(traxel_map[n]) );

      table.add_to( *(m.opengm_model) );
    }

    inline void ChaingraphModelBuilder::add_outgoing_factor( TrackingModel& m, 
			      const HypothesesGraph::Node& n ) const {
      using namespace std;

      LOG(logDEBUG) << "ChaingraphModelBuilder::add_outgoing_factor(): entered";
      property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = hypotheses_->get(node_traxel());
      // collect and count outgoing arcs
      vector<HypothesesGraph::Arc> arcs; 
      vector<size_t> vi; 		// opengm variable indeces
      vi.push_back(m.node_var[n]); // first detection node, remaining will be transition nodes
      int count = 0;
      for(HypothesesGraph::OutArcIt a(*hypotheses_, n); a != lemon::INVALID; ++a) {
	arcs.push_back(a);
	vi.push_back(m.arc_var[a]);
	++count;
      }

      // construct factor
      if(count == 0) {
	// build value table
	size_t table_dim = 1; 		// only one detection var
	std::vector<size_t> coords;
	OpengmExplicitFactor<double> table( vi );

	// opportunity
	coords = std::vector<size_t>(table_dim, 0); 		// (0)
	table.set_value( coords, opportunity_cost_ );

	// disappearance 
	coords = std::vector<size_t>(table_dim, 0);
	coords[0] = 1; 						// (1)
	table.set_value( coords, disappearance_(traxel_map[n]) );

	table.add_to( *m.opengm_model );

      } else if(count == 1) {
	// no division possible
	size_t table_dim = 2; 		// detection var + 1 * transition var
	std::vector<size_t> coords;
	OpengmExplicitFactor<double> table( vi, forbidden_cost_ );

	// opportunity configuration
	coords = std::vector<size_t>(table_dim, 0); // (0,0)
	table.set_value( coords, opportunity_cost_ );

	// disappearance configuration
	coords = std::vector<size_t>(table_dim, 0);
	coords[0] = 1; // (1,0)
	table.set_value( coords, disappearance_(traxel_map[n]) );

	// move configurations
	coords = std::vector<size_t>(table_dim, 1);
	// (1,1)
	table.set_value( coords, move_(traxel_map[n], traxel_map[hypotheses_->target(arcs[0])]) );

	table.add_to( *m.opengm_model );

      } else {
	// build value table
	size_t table_dim = count + 1; 		// detection var + n * transition var
	std::vector<size_t> coords;
	OpengmExplicitFactor<double> table( vi, forbidden_cost_ );

	// opportunity configuration
	coords = std::vector<size_t>(table_dim, 0); // (0,0,...,0)
	table.set_value( coords, opportunity_cost_ );

	// disappearance configuration
	coords = std::vector<size_t>(table_dim, 0);
	coords[0] = 1; // (1,0,...,0)
	table.set_value( coords, disappearance_(traxel_map[n]) );

	// move configurations
	coords = std::vector<size_t>(table_dim, 0);
	coords[0] = 1;
	// (1,0,0,0,1,0,0)
	for(size_t i = 1; i < table_dim; ++i) {
	  coords[i] = 1; 
	  table.set_value( coords, move_(traxel_map[n], traxel_map[hypotheses_->target(arcs[i-1])]) );
	  coords[i] = 0; // reset coords
	}
      
	// division configurations
	coords = std::vector<size_t>(table_dim, 0);
	coords[0] = 1;
	// (1,0,0,1,0,1,0,0) 
	for(unsigned int i = 1; i < table_dim - 1; ++i) {
	  for(unsigned int j = i+1; j < table_dim; ++j) {
	    coords[i] = 1;
	    coords[j] = 1;
	    table.set_value(coords, division_(traxel_map[n],
					      traxel_map[hypotheses_->target(arcs[i-1])],
					      traxel_map[hypotheses_->target(arcs[j-1])]
					      ));
	  
	    // reset
	    coords[i] = 0;
	    coords[j] = 0;
	  }
	}

	table.add_to( *m.opengm_model );      

      }   
      LOG(logDEBUG) << "ChaingraphModelBuilder::add_outgoing_factor(): leaving";
    }

    inline void ChaingraphModelBuilder::add_incoming_factor( TrackingModel& m,
			      const HypothesesGraph::Node& n) const {
      using namespace std;

      LOG(logDEBUG) << "ChaingraphModelBuilder::add_incoming_factor(): entered";
      property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = hypotheses_->get(node_traxel());
      // collect and count incoming arcs
      vector<size_t> vi; // opengm variable indeces
      int count = 0;
      for(HypothesesGraph::InArcIt a(*hypotheses_, n); a != lemon::INVALID; ++a) {
	vi.push_back(m.arc_var[a]);
	++count;
      }
      vi.push_back(m.node_var[n]); 
      std::reverse(vi.begin(), vi.end());
    
      //// construct factor
      // build value table
      size_t table_dim = count + 1; // detection var + n * transition var
      OpengmExplicitFactor<double> table( vi, forbidden_cost_ );
      std::vector<size_t> coords;

      // allow opportunity configuration
      // (0,0,...,0)
      coords = std::vector<size_t>(table_dim, 0);
      table.set_value( coords, 0 );

      // appearance configuration
      coords = std::vector<size_t>(table_dim, 0);
      coords[0] = 1; // (1,0,...,0)
      table.set_value( coords, appearance_(traxel_map[n]) );
      assert(table.get_value( coords ) == appearance_(traxel_map[n]));

      // allow move configurations
      coords = std::vector<size_t>(table_dim, 0);
      coords[0] = 1;
      // (1,0,0,0,1,0,0)
      for(size_t i = 1; i < table_dim; ++i) {
	coords[i] = 1; 
	table.set_value( coords, 0 );
	coords[i] = 0; // reset coords
      }

      table.add_to( *m.opengm_model );
      LOG(logDEBUG) << "ChaingraphModelBuilder::add_incoming_factor(): leaving";
    }

    void ChaingraphModelBuilder::couple(const TrackingModel& m, const HypothesesGraph::Node& n, const HypothesesGraph::Arc& a, OpengmLPCplex& cplex ) {
      vector<size_t> cplex_idxs; 
      cplex_idxs.push_back(cplex_id(m.node_var.find(n)->second));
      cplex_idxs.push_back(cplex_id(m.arc_var.find(a)->second));
      vector<int> coeffs;
      coeffs.push_back(1);
      coeffs.push_back(-1);
      // 0 <= 1*detection - 1*transition <= 1
      cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin() , 0, 1);
    }

    void ChaingraphModelBuilder::fix_detections( const TrackingModel& m, const HypothesesGraph& g, OpengmLPCplex& cplex ) {
      for(HypothesesGraph::NodeIt n(g); n!=lemon::INVALID; ++n) {
	vector<size_t> cplex_idxs; 
	cplex_idxs.push_back(cplex_id(m.node_var.find(n)->second));
	vector<int> coeffs;
	coeffs.push_back(1);
	// 1 <= 1*detection <= 1
	cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin() , 1, 1);
      }
    }
    
  } /* namespace pgm */



  ////
  //// class Chaingraph
  ////

  Chaingraph::~Chaingraph() {
    if(optimizer_ != NULL) {
	delete optimizer_;
	optimizer_ = NULL;
    }
  }

double Chaingraph::forbidden_cost() const {
    return forbidden_cost_;
}

bool Chaingraph::with_constraints() const {
    return with_constraints_;
}

namespace {
   struct NullDeleter {template<typename T> void operator()(T*) {} };
}

void Chaingraph::formulate( const HypothesesGraph& hypotheses ) {
    LOG(logDEBUG) << "Chaingraph::formulate: entered";
    reset();

    // configure the model builder
    shared_ptr<const HypothesesGraph> g(&hypotheses, NullDeleter() );
    pgm::ChaingraphModelBuilder builder =
      pgm::ChaingraphModelBuilder(g, appearance_,disappearance_,move_,opportunity_cost_, forbidden_cost_)
      .with_detection_vars( detection_, non_detection_ )
      .with_divisions( division_ );

    // build the model
    tracking_model_ = builder.build();

    // refine the model with hard constraints
    pgm::OpengmLPCplex::Parameter param;
    param.verbose_ = true;
    param.integerConstraint_ = true;
    param.epGap_ = ep_gap_;
    LOG(logDEBUG) << "Chaingraph::formulate ep_gap = " << param.epGap_;
    pgm::OpengmLPCplex* cplex = new pgm::OpengmLPCplex(*(tracking_model_->opengm_model), param);
    optimizer_ = cplex; // opengm::Inference optimizer_

    if (with_constraints_) {
      LOG(logDEBUG) << "Chaingraph::formulate: add_constraints";
      pgm::ChaingraphModelBuilder::add_hard_constraints( *tracking_model_ , hypotheses, *cplex );
    }
    
    if (fixed_detections_) {
      LOG(logDEBUG) << "Chaingraph::formulate: fix_detections";
      pgm::ChaingraphModelBuilder::fix_detections( *tracking_model_, hypotheses, *cplex );
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
    vector<OpengmModel<>::ogmInference::LabelType> solution;
    opengm::InferenceTermination status = optimizer_->arg(solution);
    if(status != opengm::NORMAL) {
	throw runtime_error("GraphicalModel::infer(): solution extraction terminated unnormally");
    }

    // add 'active' properties to graph
    g.add(node_active()).add(arc_active());
    property_map<node_active, HypothesesGraph::base_graph>::type& active_nodes = g.get(node_active());
    property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = g.get(arc_active());

    // write state after inference into 'active'-property maps
    for(std::map<HypothesesGraph::Node, size_t>::const_iterator it = tracking_model_->node_var.begin(); it != tracking_model_->node_var.end(); ++it) {
	bool state = false;
	if(solution[it->second] == 1) state = true;
	active_nodes.set(it->first, state);
    }
    for(std::map<HypothesesGraph::Arc, size_t>::const_iterator it = tracking_model_->arc_var.begin(); it != tracking_model_->arc_var.end(); ++it) {
	bool state = false;
	if(solution[it->second] == 1) state = true;
	active_arcs.set(it->first, state);
    }
}

  const OpengmModel<>::ogmGraphicalModel* Chaingraph::get_graphical_model() const {
    return tracking_model_->opengm_model.get();
  }

  const std::map<HypothesesGraph::Node, size_t>& Chaingraph::get_node_map() const {
    return tracking_model_->node_var;
  }

  const std::map<HypothesesGraph::Arc, size_t>& Chaingraph::get_arc_map() const {
    return tracking_model_->arc_var;
  }

void Chaingraph::reset() {
    if(optimizer_ != NULL) {
	delete optimizer_;
	optimizer_ = NULL;
    }
}

} /* namespace pgmlink */ 
