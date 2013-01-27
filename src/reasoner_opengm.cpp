#include <algorithm>
#include <cassert>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>
#include <boost/scoped_ptr.hpp>
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
    void FactorEntry::set( OpengmModel::ValueType v ) {
      if( entry_ == NULL ) {
	throw runtime_error("FactorEntry: is invalid");
      }
      *entry_ = v;
    }

    OpengmModel::ValueType FactorEntry::get() const {
      if( entry_ == NULL ) {
	throw runtime_error("FactorEntryPtr: is invalid");
      }
      return *entry_;
    } 

    ////
    //// class ChaingraphModel
    ////
    ChaingraphModel::ChaingraphModel() : opengm_model( new OpengmModel() ) {
      init();
    }

    ChaingraphModel::ChaingraphModel(shared_ptr<OpengmModel> m,
				     const node_var_map& node_var,
				     const arc_var_map& arc_var
				     )
      : opengm_model(m) {
      node_var_.left = node_var;
      arc_var_.left = arc_var;
      init();
      }

    const ChaingraphModel::node_var_map& ChaingraphModel::var_of_node() const {
      return node_var_.left;
    }

    const ChaingraphModel::var_node_map& ChaingraphModel::node_of_var() const {
      return node_var_.right;
    }

    const ChaingraphModel::arc_var_map& ChaingraphModel::var_of_arc() const {
      return arc_var_.left;
    }

    const ChaingraphModel::var_arc_map& ChaingraphModel::arc_of_var() const
    {
      return arc_var_.right;      
    }

    ChaingraphModel::var_t ChaingraphModel::var_of_node(node_t e) const {
      node_var_map::const_iterator it = var_of_node().find(e);
      if(it!=var_of_node().end()) {
	return it->second;
      } else {
	throw std::out_of_range("ChaingraphModel::var_of_node(): key does not exist");
      }
    }

    ChaingraphModel::var_t ChaingraphModel::var_of_arc(arc_t e) const {
      arc_var_map::const_iterator it = var_of_arc().find(e);
      if(it!=var_of_arc().end()) {
	return it->second;
      } else {
	throw std::out_of_range("ChaingraphModel::var_of_arc(): key does not exist");
      }
    }

    ChaingraphModel::node_t ChaingraphModel::node_of_var(var_t e) const {
      var_node_map::const_iterator it = node_of_var().find(e);
      if(it!=node_of_var().end()) {
	return it->second;
      } else {
	throw std::out_of_range("ChaingraphModel::node_of_var(): key does not exist");
      }
    }

    ChaingraphModel::arc_t ChaingraphModel::arc_of_var(var_t e) const {
      var_arc_map::const_iterator it = arc_of_var().find(e);
      if(it!=arc_of_var().end()) {
	return it->second;
      } else {
	throw std::out_of_range("ChaingraphModel::arc_of_var(): key does not exist");
      }
    }

    ChaingraphModel::VarCategory ChaingraphModel::var_category(var_t e) const {
      if(arc_of_var().count(e)) {
	return ChaingraphModel::arc_var;
      } else if(node_of_var().count(e)) {
	return ChaingraphModel::node_var;
      } else {
	throw std::out_of_range("ChaingraphModel::var_category(): key does not exist");
      }
    }


    void ChaingraphModel::init() {
      weight_map[det_weight] = vector<OpengmModel::IndexType>();
      weight_map[mov_weight] = vector<OpengmModel::IndexType>();
      weight_map[div_weight] = vector<OpengmModel::IndexType>();
      weight_map[app_weight] = vector<OpengmModel::IndexType>();
      weight_map[dis_weight] = vector<OpengmModel::IndexType>();
      weight_map[opp_weight] = vector<OpengmModel::IndexType>();
    }

    
    ////
    //// class ChaingraphModelBuilder
    ////
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

    namespace {
      inline size_t cplex_id(size_t opengm_id) {
	return 2*opengm_id + 1;
      }
    }
    void ChaingraphModelBuilder::add_hard_constraints( const ChaingraphModel& m, const HypothesesGraph& hypotheses, OpengmLPCplex& cplex ) {
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
	  cplex_idxs.push_back(cplex_id(m.var_of_arc(a)));
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
	  cplex_idxs.push_back(cplex_id(m.var_of_arc(a)));
	}
	vector<int> coeffs(cplex_idxs.size(), 1);
	// 0 <= 1*transition + ... + 1*transition <= 1
	cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 1);
      }
    }

    void ChaingraphModelBuilder::fix_detections( const ChaingraphModel& m, const HypothesesGraph& g, OpengmLPCplex& cplex ) {
      for(HypothesesGraph::NodeIt n(g); n!=lemon::INVALID; ++n) {
	vector<size_t> cplex_idxs; 
	cplex_idxs.push_back(cplex_id(m.var_of_node(n)));
	vector<int> coeffs;
	coeffs.push_back(1);
	// 1 <= 1*detection <= 1
	cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin() , 1, 1);
      }
    }

    inline void ChaingraphModelBuilder::add_detection_vars( const HypothesesGraph& hypotheses, ChaingraphModel& m ) const {
      for(HypothesesGraph::NodeIt n(hypotheses); n!=lemon::INVALID; ++n) {
	m.opengm_model->addVariable(2);
	m.node_var_.left.insert(ChaingraphModel::node_var_map::value_type(n, m.opengm_model->numberOfVariables() - 1));
      }
    }

    inline void ChaingraphModelBuilder::add_assignment_vars( const HypothesesGraph& hypotheses, ChaingraphModel& m ) const {
      for(HypothesesGraph::ArcIt a(hypotheses); a!=lemon::INVALID; ++a) {
	m.opengm_model->addVariable(2);
	m.arc_var_.left.insert(ChaingraphModel::arc_var_map::value_type(a, m.opengm_model->numberOfVariables() - 1));
      }
    }

    void ChaingraphModelBuilder::couple(const ChaingraphModel& m, const HypothesesGraph::Node& n, const HypothesesGraph::Arc& a, OpengmLPCplex& cplex ) {
      vector<size_t> cplex_idxs; 
      cplex_idxs.push_back(cplex_id(m.var_of_node(n)));
      cplex_idxs.push_back(cplex_id(m.var_of_arc(a)));
      vector<int> coeffs;
      coeffs.push_back(1);
      coeffs.push_back(-1);
      // 0 <= 1*detection - 1*transition <= 1
      cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin() , 0, 1);
    }



    ////
    //// class TrainableChaingraphModelBuilder
    ////
    TrainableChaingraphModelBuilder* TrainableChaingraphModelBuilder::clone() const {
      return new TrainableChaingraphModelBuilder(*this);
    }

    ChaingraphModel* TrainableChaingraphModelBuilder::build(const HypothesesGraph& hypotheses) const {

      if( !has_detection_vars() ) {
	throw std::runtime_error("TrainableChaingraphModelBuilder::build: option without detection vars not yet implemented");
      }
      if( !has_divisions() ) {
	throw std::runtime_error("TrainableChaingraphModelBuilder::build: option without divisions not yet implemented");
      }

      //// setup the model
      ChaingraphModel* model( new ChaingraphModel() );

      // we need six weights; one per event
      assert(model->opengm_model->numberOfWeights() == 0);
      model->opengm_model->increaseNumberOfWeights(6);

      // assign the weight ids to event types
      model->weight_map[ChaingraphModel::det_weight].push_back(0);
      model->weight_map[ChaingraphModel::mov_weight].push_back(1);
      model->weight_map[ChaingraphModel::div_weight].push_back(2);
      model->weight_map[ChaingraphModel::app_weight].push_back(3);
      model->weight_map[ChaingraphModel::dis_weight].push_back(4);
      model->weight_map[ChaingraphModel::opp_weight].push_back(5);

      
      if( has_detection_vars() ) {
	add_detection_vars( hypotheses, *model );
      }
      add_assignment_vars( hypotheses, *model );

      if( has_detection_vars() ) {
      	for(HypothesesGraph::NodeIt n(hypotheses); n!=lemon::INVALID; ++n) {
      	  add_detection_factor( hypotheses, *model, n );
      	}
      }

      for(HypothesesGraph::NodeIt n(hypotheses); n!=lemon::INVALID; ++n) {
      	add_outgoing_factor( hypotheses, *model, n );
      	add_incoming_factor( hypotheses, *model, n );
      }

      return model;
    }

    void TrainableChaingraphModelBuilder::add_detection_factor( const HypothesesGraph& hypotheses, ChaingraphModel& m, const HypothesesGraph::Node& n ) const {
      property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = hypotheses.get(node_traxel());
      std::vector<size_t> var_indices;
      var_indices.push_back(m.var_of_node(n));
      size_t shape[] = {2};

      size_t indicate[] = {0};
      OpengmWeightedFeature<OpengmModel::ValueType>(var_indices, shape, shape+1, indicate, non_detection()(traxel_map[n]) )
      	.add_as_feature_to( *(m.opengm_model), m.weight_map[ChaingraphModel::det_weight].front() );

      indicate[0] = 1;
      OpengmWeightedFeature<OpengmModel::ValueType>(var_indices, shape, shape+1, indicate, detection()(traxel_map[n]) )
      	.add_as_feature_to( *(m.opengm_model), m.weight_map[ChaingraphModel::det_weight].front() );
    }

    namespace {
      std::vector<size_t> DecToBin(size_t number)
      {
	if ( number == 0 ) return std::vector<size_t>(1,0);
	if ( number == 1 ) return std::vector<size_t>(1,1);
	
	if ( number % 2 == 0 ) {
	  std::vector<size_t> ret = DecToBin(number / 2);
	  ret.push_back(0);
	  return ret;
	}
	else {
	  std::vector<size_t> ret = DecToBin(number / 2);
	  ret.push_back(1);
	  return ret;
	}
      }

      int BinToDec(std::vector<size_t> number)
      {
	int result = 0, pow = 1;
	for ( int i = number.size() - 1; i >= 0; --i, pow <<= 1 )
	  result += number[i] * pow;
	return result;
      }
    }

    inline void TrainableChaingraphModelBuilder::add_outgoing_factor( const HypothesesGraph& hypotheses,
								      ChaingraphModel& m, 
								      const HypothesesGraph::Node& n) const {
      using namespace std;
      property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = hypotheses.get(node_traxel());

      LOG(logDEBUG) << "TrainableChaingraphModelBuilder::add_outgoing_factor(): entered";
      // collect and count outgoing arcs
      vector<HypothesesGraph::Arc> arcs; 
      vector<size_t> vi; 		// opengm variable indeces
      vi.push_back(m.var_of_node(n)); // first detection node, remaining will be transition nodes
      int count = 0;
      for(HypothesesGraph::OutArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
	arcs.push_back(a);
	vi.push_back(m.var_of_arc(a));
	++count;
      }

      // construct factor
      if(count == 0) {
	size_t shape[] = {2}; // only one detection var
	size_t indicate[] = {0};

	OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape, shape+1, indicate, opportunity_cost() )
	  .add_as_feature_to( *(m.opengm_model), m.weight_map[ChaingraphModel::opp_weight].front() );

	indicate[0] = 1;
	OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape, shape+1, indicate, disappearance()(traxel_map[n]) )
	  .add_as_feature_to( *(m.opengm_model), m.weight_map[ChaingraphModel::dis_weight].front() );

      } else if(count == 1) {
	// no division possible
	assert( vi.size() == 2 );
	size_t shape[] = {2,2};	// detection var + 1 * transition var
	std::vector<size_t> coords;

	// opportunity configuration
	coords = std::vector<size_t>(2, 0); // (0,0)
	OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape, shape+2, coords.begin(), opportunity_cost() )
	  .add_as_feature_to( *(m.opengm_model), m.weight_map[ChaingraphModel::opp_weight].front() );

	// disappearance configuration
	coords = std::vector<size_t>(2, 0);
	coords[0] = 1; // (1,0)
	OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape, shape+2, coords.begin(), disappearance()(traxel_map[n]) )
	  .add_as_feature_to( *(m.opengm_model), m.weight_map[ChaingraphModel::dis_weight].front() );

	// move configurations
	coords = std::vector<size_t>(2, 1); // (1,1)
	OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape, shape+2, coords.begin(), move()(traxel_map[n], traxel_map[hypotheses.target(arcs[0])]) )
	  .add_as_feature_to( *(m.opengm_model), m.weight_map[ChaingraphModel::mov_weight].front() );

	// forbidden configuration
	coords = std::vector<size_t>(2, 0); // (0,1)
	coords[1] = 1;
	OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape, shape+2, coords.begin(), forbidden_cost() )
	  .add_to( *(m.opengm_model) );

      } else {
	size_t table_dim = count + 1; 		// detection var + n * transition var
	std::vector<size_t> coords;
	std::vector<size_t> shape(table_dim, 2);

	std::set<size_t > entries;
	for(size_t i = 0; i < static_cast<size_t>(std::pow(2, table_dim)); ++i) {
	  entries.insert( entries.end(), i );
	}
	
	// opportunity configuration
	coords = std::vector<size_t>(table_dim, 0); // (0,0,...,0)
	size_t check = entries.erase(BinToDec(coords));
	assert(check == 1);
	OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape.begin(), shape.end(), coords.begin(), opportunity_cost() )
	  .add_as_feature_to( *(m.opengm_model), m.weight_map[ChaingraphModel::opp_weight].front() );

	// disappearance configuration
	coords = std::vector<size_t>(table_dim, 0);
	coords[0] = 1; // (1,0,...,0)
	check = entries.erase(BinToDec(coords));
	assert(check == 1);
	OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape.begin(), shape.end(), coords.begin(), disappearance()(traxel_map[n]) )
	  .add_as_feature_to( *(m.opengm_model), m.weight_map[ChaingraphModel::dis_weight].front() );

	// move configurations
	coords = std::vector<size_t>(table_dim, 0);
	coords[0] = 1;
	// (1,0,0,0,1,0,0)
	for(size_t i = 1; i < table_dim; ++i) {
	  coords[i] = 1; 
	  check = entries.erase(BinToDec(coords));
	  assert(check == 1);
	  OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape.begin(), shape.end(), coords.begin(), move()(traxel_map[n], traxel_map[hypotheses.target(arcs[i-1])]) )
	    .add_as_feature_to( *(m.opengm_model), m.weight_map[ChaingraphModel::mov_weight].front() );

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
	    check = entries.erase(BinToDec(coords));
	    assert(check == 1);
	    OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape.begin(), shape.end(), coords.begin(), division()(traxel_map[n],
					      traxel_map[hypotheses.target(arcs[i-1])],
					      traxel_map[hypotheses.target(arcs[j-1])]
					      ) )
	      .add_as_feature_to( *(m.opengm_model), m.weight_map[ChaingraphModel::div_weight].front() );
	  
	    // reset
	    coords[i] = 0;
	    coords[j] = 0;
	  }
	}

	// forbidden configurations
	for(std::set<size_t>::iterator it = entries.begin(); it != entries.end(); ++it) {
	  coords = DecToBin(*it);
	  if(coords.size() < table_dim ) {
	    coords.insert(coords.begin(), table_dim-coords.size(), 0);
	  }
	  assert( coords.size() == table_dim );
	  OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape.begin(), shape.end(), coords.begin(), forbidden_cost())
	    .add_to( *(m.opengm_model) );
	}

      }   
      LOG(logDEBUG) << "TrainableChaingraphModelBuilder::add_outgoing_factor(): leaving";
    }

    inline void TrainableChaingraphModelBuilder::add_incoming_factor( const HypothesesGraph& hypotheses,
								      ChaingraphModel& m,
								      const HypothesesGraph::Node& n ) const {
      using namespace std;
      property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = hypotheses.get(node_traxel());

      LOG(logDEBUG) << "TrainableChaingraphModelBuilder::add_incoming_factor(): entered";
      // collect and count incoming arcs
      vector<size_t> vi; // opengm variable indeces
      int count = 0;
      for(HypothesesGraph::InArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
	vi.push_back(m.var_of_arc(a));
	++count;
      }
      vi.push_back(m.var_of_node(n)); 
      std::reverse(vi.begin(), vi.end());
    
      //// construct factor
      size_t table_dim = count + 1; 		// detection var + n * transition var
      std::vector<size_t> coords;
      std::vector<size_t> shape(table_dim, 2);
      
      std::set<size_t > entries;
      for(size_t i = 0; i < static_cast<size_t>(std::pow(2, table_dim)); ++i) {
	entries.insert( entries.end(), i );
      }


      // allow opportunity configuration
      // (0,0,...,0)
      coords = std::vector<size_t>(table_dim, 0);
      size_t check = entries.erase(BinToDec(coords));
      assert(check == 1);

      // appearance configuration
      coords = std::vector<size_t>(table_dim, 0);
      coords[0] = 1; // (1,0,...,0)
	check = entries.erase(BinToDec(coords));
	assert(check == 1);
	OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape.begin(), shape.end(), coords.begin(), appearance()(traxel_map[n]) )
	  .add_as_feature_to( *(m.opengm_model), m.weight_map[ChaingraphModel::app_weight].front() );

      // allow move configurations
      coords = std::vector<size_t>(table_dim, 0);
      coords[0] = 1;
      // (1,0,0,0,1,0,0)
      for(size_t i = 1; i < table_dim; ++i) {
	coords[i] = 1; 
	size_t check = entries.erase(BinToDec(coords));
	assert(check == 1);
	coords[i] = 0; // reset coords
      }

      // forbidden configurations
      for(std::set<size_t>::iterator it = entries.begin(); it != entries.end(); ++it) {
	coords = DecToBin(*it);
	if(coords.size() < table_dim ) {
	  coords.insert(coords.begin(), table_dim-coords.size(), 0);
	}
	assert( coords.size() == table_dim );
	OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape.begin(), shape.end(), coords.begin(), forbidden_cost())
	  .add_to( *(m.opengm_model) );
      }
      
      LOG(logDEBUG) << "TrainableChaingraphModelBuilder::add_incoming_factor(): leaving";
    }



    ////
    //// class ChaingraphModelBuilderECCV12
    ////
    ChaingraphModelBuilderECCV12* ChaingraphModelBuilderECCV12::clone() const {
      return new ChaingraphModelBuilderECCV12(*this);
    }

    ChaingraphModel* ChaingraphModelBuilderECCV12::build(const HypothesesGraph& hypotheses) const {
      using boost::shared_ptr;
      using std::map;

      LOG(logDEBUG) << "ChaingraphModelBuilderECCV12::build: entered";
      if( !has_detection_vars() ) {
	throw std::runtime_error("ChaingraphModelBuilderECCV12::build(): option without detection vars not yet implemented");
      }
      if( !has_divisions() ) {
	throw std::runtime_error("ChaingraphModelBuilderECCV12::build(): option without divisions not yet implemented");
      }

      ChaingraphModel* model( new ChaingraphModel() );
      
      if( has_detection_vars() ) {
	add_detection_vars( hypotheses, *model );
      }
      add_assignment_vars( hypotheses, *model );


      if( has_detection_vars() ) {
	for(HypothesesGraph::NodeIt n(hypotheses); n!=lemon::INVALID; ++n) {
	  add_detection_factor( hypotheses, *model, n );
	}
      }

      for(HypothesesGraph::NodeIt n(hypotheses); n!=lemon::INVALID; ++n) {
	add_outgoing_factor( hypotheses, *model, n );
 	add_incoming_factor( hypotheses, *model, n );
      }

      return model;
    }

    void ChaingraphModelBuilderECCV12::add_detection_factor( const HypothesesGraph& hypotheses, ChaingraphModel& m, const HypothesesGraph::Node& n) const {
      property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = hypotheses.get(node_traxel());

      size_t vi[] = {m.var_of_node(n)};
      vector<size_t> coords(1,0);
      OpengmExplicitFactor<double> table( vi, vi+1 );

      coords[0] = 0;
      table.set_value( coords, non_detection()(traxel_map[n]) ); 

      coords[0] = 1;
      table.set_value( coords, detection()(traxel_map[n]) );

      table.add_to( *(m.opengm_model) );
    }

    inline void ChaingraphModelBuilderECCV12::add_outgoing_factor( const HypothesesGraph& hypotheses,
								   ChaingraphModel& m, 
								   const HypothesesGraph::Node& n
								   ) const {
      using namespace std;
      property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = hypotheses.get(node_traxel());

      LOG(logDEBUG) << "ChaingraphModelBuilderECCV12::add_outgoing_factor(): entered";
      // collect and count outgoing arcs
      vector<HypothesesGraph::Arc> arcs; 
      vector<size_t> vi; 		// opengm variable indeces
      vi.push_back(m.var_of_node(n)); // first detection node, remaining will be transition nodes
      int count = 0;
      for(HypothesesGraph::OutArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
	arcs.push_back(a);
	vi.push_back(m.var_of_arc(a));
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
	table.set_value( coords, opportunity_cost() );

	// disappearance 
	coords = std::vector<size_t>(table_dim, 0);
	coords[0] = 1; 						// (1)
	table.set_value( coords, disappearance()(traxel_map[n]) );

	table.add_to( *m.opengm_model );

      } else if(count == 1) {
	// no division possible
	size_t table_dim = 2; 		// detection var + 1 * transition var
	std::vector<size_t> coords;
	OpengmExplicitFactor<double> table( vi, forbidden_cost() );

	// opportunity configuration
	coords = std::vector<size_t>(table_dim, 0); // (0,0)
	table.set_value( coords, opportunity_cost() );

	// disappearance configuration
	coords = std::vector<size_t>(table_dim, 0);
	coords[0] = 1; // (1,0)
	table.set_value( coords, disappearance()(traxel_map[n]) );

	// move configurations
	coords = std::vector<size_t>(table_dim, 1);
	// (1,1)
	table.set_value( coords, move()(traxel_map[n], traxel_map[hypotheses.target(arcs[0])]) );

	table.add_to( *m.opengm_model );

      } else {
	// build value table
	size_t table_dim = count + 1; 		// detection var + n * transition var
	std::vector<size_t> coords;
	OpengmExplicitFactor<double> table( vi, forbidden_cost() );

	// opportunity configuration
	coords = std::vector<size_t>(table_dim, 0); // (0,0,...,0)
	table.set_value( coords, opportunity_cost() );

	// disappearance configuration
	coords = std::vector<size_t>(table_dim, 0);
	coords[0] = 1; // (1,0,...,0)
	table.set_value( coords, disappearance()(traxel_map[n]) );

	// move configurations
	coords = std::vector<size_t>(table_dim, 0);
	coords[0] = 1;
	// (1,0,0,0,1,0,0)
	for(size_t i = 1; i < table_dim; ++i) {
	  coords[i] = 1; 
	  table.set_value( coords, move()(traxel_map[n], traxel_map[hypotheses.target(arcs[i-1])]) );
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
	    table.set_value(coords, division()(traxel_map[n],
					      traxel_map[hypotheses.target(arcs[i-1])],
					      traxel_map[hypotheses.target(arcs[j-1])]
					      ));
	  
	    // reset
	    coords[i] = 0;
	    coords[j] = 0;
	  }
	}

	table.add_to( *m.opengm_model );      

      }   
      LOG(logDEBUG) << "ChaingraphModelBuilderECCV12::add_outgoing_factor(): leaving";
    }

    inline void ChaingraphModelBuilderECCV12::add_incoming_factor( const HypothesesGraph& hypotheses,
								   ChaingraphModel& m,
								   const HypothesesGraph::Node& n) const {
      using namespace std;
      property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = hypotheses.get(node_traxel());

      LOG(logDEBUG) << "ChaingraphModelBuilderECCV12::add_incoming_factor(): entered";
      // collect and count incoming arcs
      vector<size_t> vi; // opengm variable indeces
      int count = 0;
      for(HypothesesGraph::InArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
	vi.push_back(m.var_of_arc(a));
	++count;
      }
      vi.push_back(m.var_of_node(n)); 
      std::reverse(vi.begin(), vi.end());
    
      //// construct factor
      // build value table
      size_t table_dim = count + 1; // detection var + n * transition var
      OpengmExplicitFactor<double> table( vi, forbidden_cost() );
      std::vector<size_t> coords;

      // allow opportunity configuration
      // (0,0,...,0)
      coords = std::vector<size_t>(table_dim, 0);
      table.set_value( coords, 0 );

      // appearance configuration
      coords = std::vector<size_t>(table_dim, 0);
      coords[0] = 1; // (1,0,...,0)
      table.set_value( coords, appearance()(traxel_map[n]) );
      assert(table.get_value( coords ) == appearance()(traxel_map[n]));

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
      LOG(logDEBUG) << "ChaingraphModelBuilderECCV12::add_incoming_factor(): leaving";
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
    linking_model_ = boost::shared_ptr<pgm::ChaingraphModel>(builder_->build(hypotheses));

    // refine the model with hard constraints
    pgm::OpengmLPCplex::Parameter param;
    param.verbose_ = true;
    param.integerConstraint_ = true;
    param.epGap_ = ep_gap_;
    LOG(logDEBUG) << "Chaingraph::formulate ep_gap = " << param.epGap_;
    pgm::OpengmLPCplex* cplex = new pgm::OpengmLPCplex(*(linking_model_->opengm_model), param);
    optimizer_ = cplex; // opengm::Inference optimizer_

    if (with_constraints_) {
      LOG(logDEBUG) << "Chaingraph::formulate: add_constraints";
      pgm::ChaingraphModelBuilder::add_hard_constraints( *linking_model_ , hypotheses, *cplex );
    }
    
    if (fixed_detections_) {
      LOG(logDEBUG) << "Chaingraph::formulate: fix_detections";
      pgm::ChaingraphModelBuilder::fix_detections( *linking_model_, hypotheses, *cplex );
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
    for(pgm::ChaingraphModel::node_var_map::const_iterator it = linking_model_->var_of_node().begin(); it != linking_model_->var_of_node().end(); ++it) {
	bool state = false;
	if(solution[it->second] == 1) state = true;
	active_nodes.set(it->first, state);
    }
    for(pgm::ChaingraphModel::arc_var_map::const_iterator it = linking_model_->var_of_arc().begin(); it != linking_model_->var_of_arc().end(); ++it) {
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
