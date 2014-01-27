// stl
#include <algorithm>
#include <cassert>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>
#include <iterator>
#include <functional>
#include <sstream>

// boost
#include <boost/shared_ptr.hpp>

// opengm
#include <opengm/inference/lpcplex.hxx>
#include <opengm/datastructures/marray/marray.hxx>

// pgmlink
#include "pgmlink/pgm.h"
#include "pgmlink/log.h"
#include "pgmlink/traxels.h"
#include "pgmlink/pgm_multi_hypotheses.h"
#include "pgmlink/multi_hypotheses_graph.h"
#include "pgmlink/nearest_neighbors.h"
#include "pgmlink/feature.h"

//#include <ostream>

// using namespace std;

namespace pgmlink {
namespace pgm {
namespace multihypotheses {
////
//// class multihypotheses::Model
////
Model::Model() : opengm_model( new OpengmModel() ) {
  init();
}

Model::Model(shared_ptr<OpengmModel> m,
             const node_var_map& node_var,
             const arc_var_map& arc_var
             )
    : opengm_model(m) {
  node_var_.left = node_var;
  arc_var_.left = arc_var;
  init();
}

const Model::node_var_map& Model::var_of_node() const {
  return node_var_.left;
}

const Model::var_node_map& Model::node_of_var() const {
  return node_var_.right;
}

const Model::arc_var_map& Model::var_of_arc() const {
  return arc_var_.left;
}

const Model::var_arc_map& Model::arc_of_var() const
{
  return arc_var_.right;      
}

Model::var_t Model::var_of_node(const node_t& e) const {
  node_var_map::const_iterator it = var_of_node().find(e);
  if(it!=var_of_node().end()) {
    return it->second;
  } else {
    throw std::out_of_range("MultiHypothesesModel::var_of_node(): key does not exist");
  }
}

Model::var_t Model::var_of_arc(const arc_t& e) const {
  arc_var_map::const_iterator it = var_of_arc().find(e);
  if(it!=var_of_arc().end()) {
    return it->second;
  } else {
    throw std::out_of_range("MultiHypothesesModel::var_of_arc(): key does not exist");
  }
}

Model::node_t Model::node_of_var(var_t e) const {
  var_node_map::const_iterator it = node_of_var().find(e);
  if(it!=node_of_var().end()) {
    return it->second;
  } else {
    throw std::out_of_range("MultiHypothesesModel::node_of_var(): key does not exist");
  }
}

Model::arc_t Model::arc_of_var(var_t e) const {
  var_arc_map::const_iterator it = arc_of_var().find(e);
  if(it!=arc_of_var().end()) {
    return it->second;
  } else {
    throw std::out_of_range("MultiHypothesesModel::arc_of_var(): key does not exist");
  }
}

Model::VarCategory Model::var_category(var_t e) const {
  if(arc_of_var().count(e)) {
    return Model::arc_var;
  } else if(node_of_var().count(e)) {
    return Model::node_var;
  } else {
    throw std::out_of_range("MultiHypothesesModel::var_category(): key does not exist");
  }
}


void Model::init() {
  weight_map[det_weight] = vector<OpengmModel::IndexType>();
  weight_map[mov_weight] = vector<OpengmModel::IndexType>();
  weight_map[div_weight] = vector<OpengmModel::IndexType>();
  weight_map[app_weight] = vector<OpengmModel::IndexType>();
  weight_map[dis_weight] = vector<OpengmModel::IndexType>();
}


////
//// class multihypotheses::ModelBuilder
////
ModelBuilder& ModelBuilder::appearance( function<double (const Traxel&)> app) {
  if (!app) {
    throw std::invalid_argument("MultiHypothesesModelBuilder::appearance(): empty function");
  }
  appearance_ = app;
  return *this;
}


ModelBuilder& ModelBuilder::disappearance( function<double (const Traxel&)> dis) {
  if (!dis) {
    throw std::invalid_argument("MultiHypothesesModelBuilder::disappearance(): empty function");
  }
  disappearance_ = dis;
  return *this;
}


ModelBuilder& ModelBuilder::move(function<double (const Traxel&, const Traxel&, feature_type)> mov) {
  if (!mov) {
    throw std::invalid_argument("MultiHypothesesModelBuilder::move(): empty function");
  }
  move_ = mov;
  return *this;
}


ModelBuilder& ModelBuilder::count(function<double (feature_type)> count) {
  if (!count) {
    throw std::invalid_argument("MultiHypothesesModelBuilder::count(): empty function");
  }
  count_ = count;
  return *this;
}


ModelBuilder& ModelBuilder::with_detection_vars(  function<double (const Traxel&, size_t)> detection,
                                                  function<double (const Traxel&, size_t)> non_detection) {
  if (!(detection && non_detection)) {
    throw std::invalid_argument("MultiHypothesesModelBuilder::with_detection_vars(): empty function");
  }

  with_detection_vars_ = true;
  detection_ = detection;
  non_detection_ = non_detection;
  return *this;
}


ModelBuilder& ModelBuilder::without_detection_vars() {
  with_detection_vars_ = false;
  detection_ = NULL;
  non_detection_ = NULL;
  return *this;
}

ModelBuilder& ModelBuilder::with_one_active_per_component_constraint(bool check) {
   with_one_active_per_component_constraint_ = check;
   return *this;
}

ModelBuilder& ModelBuilder::with_transition_parameter(int alpha) {
   transition_parameter_ = alpha;
   return *this;
}

ModelBuilder& ModelBuilder::with_divisions( function<double (const Traxel&, const Traxel&, const Traxel&, feature_type)> division ) {
  if (!division) {
    throw std::invalid_argument("MultiHypothesesModelBuilder::with_divisions(): empty function");
  }
  with_divisions_ = true;
  division_ = division;
  return *this;
}


ModelBuilder& ModelBuilder::without_divisions() {
  with_divisions_ = false;
  division_ = NULL;
  return *this;
}


ModelBuilder& ModelBuilder::with_classifier_priors( function<double (const Traxel&, const Traxel&, feature_type)> move,
                                                    function<double (feature_type)> count ) {
  with_classifier_priors_ = true;
  move_ = move;
  count_ = count;
  return *this;
}


ModelBuilder& ModelBuilder::without_classifier_priors() {
  with_classifier_priors_ = false;
  move_ = NULL;
  count_ = NULL;
  return *this;
}


ModelBuilder& ModelBuilder::with_maximal_conflict_cliques(bool check) {
  with_maximal_conflict_cliques_ = check;
  return *this;
}


ModelBuilder& ModelBuilder::with_hierarchical_counting_factor(bool check) {
	with_hierarchical_counting_factor_ = check;
	return *this;
}

ModelBuilder& ModelBuilder::with_counting_incoming_factor(bool check) {
	with_counting_incoming_factor_ = check;
	return *this;
}

ModelBuilder& ModelBuilder::with_maximum_arcs(unsigned maximum_outgoing_arcs) {
  maximum_outgoing_arcs_ = maximum_outgoing_arcs;
  with_maximum_arcs_ = true;
  return *this;
}


ModelBuilder& ModelBuilder::without_maximum_arcs() {
  with_maximum_arcs_ = false;
  return *this;
}


ModelBuilder& ModelBuilder::with_conflict_factors( function<double (const Traxel&, size_t)> detection) {
  if (!detection) {
    throw std::invalid_argument("MultiHypothesesModelBuilder::with_conflict_factors(): empty function");
  }
  with_detection_vars();
  detection_ = detection;
  with_conflict_factors_ = true;
  return *this;
}


ModelBuilder& ModelBuilder::with_timestep_range(int first, int last) {
  if (last < first) {
    throw std::runtime_error("Last timestep must not be smaller than first!");
  }
  first_timestep_ = first;
  last_timestep_ = last;
  with_timestep_range_ = true;
  return *this;
}


ModelBuilder& ModelBuilder::without_timestep_range() {
  with_timestep_range_ = false;
  return *this;
}


size_t ModelBuilder::cplex_id(OpengmLPCplex& cplex, const size_t opengm_id, const size_t state) const {
    return cplex.lpNodeVi(opengm_id, state);
}

size_t ModelBuilder::cplex_id(OpengmLPCplex& cplex, const size_t opengm_id) const {
    return cplex_id(cplex, opengm_id, 1);
}


void ModelBuilder::add_hard_constraints(const Model& m, const MultiHypothesesGraph& hypotheses, OpengmLPCplex& cplex) {
  LOG(logDEBUG) << "MultiHypotheses::add_hard_constraints: entered";

  ////
  //// outgoing transitions
  ////
  LOG(logDEBUG) << "MultiHypotheses::add_hard_constraints: outgoing transitions";
  for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
    if (has_detection_vars()) {
      for (MultiHypothesesGraph::OutArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
        couple(m, n, a, cplex);
      }
    }

    // couple assignments
    std::vector<size_t> cplex_idxs;
    for (HypothesesGraph::OutArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
      cplex_idxs.push_back(cplex_id(cplex, m.var_of_arc(a)));
    }
    couple_outgoing(cplex_idxs, cplex);
  }

  ////
  //// incoming transitions
  ////
  for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
    if (has_detection_vars()) {
      for (MultiHypothesesGraph::InArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
        couple(m, n, a, cplex);
      }
    }

    // couple transitions
    std::vector<size_t> cplex_idxs;
    for (MultiHypothesesGraph::InArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
      cplex_idxs.push_back(cplex_id(cplex, m.var_of_arc(a)));
    }
    couple_incoming(cplex_idxs, cplex);
  }

  ////
  //// conflicts
  ////
  if (has_detection_vars()) {
    couple_conflicts(m, hypotheses, cplex);
  }

  ////
  //// add count constraints. if none are present,
  //// var_state_coeff_constraints_ will be empty
  ////
  add_count_hard_constraints( cplex );
}


void ModelBuilder::set_cplex_timeout( double seconds ) {
  cplex_timeout_ = seconds;
}


inline void ModelBuilder::add_detection_vars( const MultiHypothesesGraph& hypotheses, Model& m ) const {
  if(!has_detection_vars()) {
    throw std::runtime_error("multihypotheses::ModelBuilder::add_detection_vars(): called without has_detection_vars()");
  }
  LOG(logDEBUG) << "ModelBuilder::add_detection_vars: entered";
  MultiHypothesesGraph::TraxelMap& traxels = hypotheses.get(node_traxel());
  for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
    const Traxel& trax = traxels[n];
    if (timestep_range_specified() &&
        (trax.Timestep < first_timestep() || trax.Timestep > last_timestep())) {
      continue;
    }
    m.opengm_model->addVariable(2);
    m.node_var_.left.insert(Model::node_var_map::value_type(n, m.opengm_model->numberOfVariables() - 1));
    LOG(logDEBUG4) << trax.Timestep << ',' << first_timestep() << ',' << last_timestep();
    LOG(logDEBUG4) << "ModelBuilder::add_detection_vars: added var " << m.opengm_model->numberOfVariables() - 1
                     << " for " << trax;
  }
}






inline void ModelBuilder::add_assignment_vars( const MultiHypothesesGraph& hypotheses, Model& m ) const {
  LOG(logINFO) << "ModelBuilder::add_assignment_vars() ";
  for (MultiHypothesesGraph::ArcIt a(hypotheses); a != lemon::INVALID; ++a) {
    m.opengm_model->addVariable(2);
    m.arc_var_.left.insert(Model::arc_var_map::value_type(a, m.opengm_model->numberOfVariables() - 1));
  }
}


std::vector<OpengmModel::IndexType> ModelBuilder::vars_for_outgoing_factor( const MultiHypothesesGraph& hypotheses,
                                                                            const Model& m,
                                                                            const MultiHypothesesGraph::Node& node) const {
  LOG(logDEBUG2) << "ModelBuilder::vars_for_outgoing_factor() -- entered";
  std::vector<OpengmModel::IndexType> vi; // opengm variable indices; may be empty if no det vars
  if (has_detection_vars()) {
    vi.push_back(m.var_of_node(node)); // first detection var, the others will be assignment vars
  }
  for (MultiHypothesesGraph::OutArcIt a(hypotheses, node); a != lemon::INVALID; ++a) {
    vi.push_back(m.var_of_arc(a));
  }
  return vi;
}


std::vector<OpengmModel::IndexType> ModelBuilder::vars_for_incoming_factor( const MultiHypothesesGraph& hypotheses,
                                                                            const Model& m,
                                                                            const MultiHypothesesGraph::Node& node) const {
  LOG(logDEBUG2) << "ModelBuilder::vars_for_incoming_factor() -- entered";
  std::vector<OpengmModel::IndexType> vi; // opengm variable indices; may be empty if no det vars
  if (has_detection_vars()) {
    vi.push_back(m.var_of_node(node)); // first detection var, the others will be assignment vars
  }
  for (MultiHypothesesGraph::InArcIt a(hypotheses, node); a != lemon::INVALID; ++a) {
    vi.push_back(m.var_of_arc(a));
  }
  return vi;
}


void ModelBuilder::couple_outgoing( const std::vector<size_t>& cplex_idxs, OpengmLPCplex& cplex ) {
  LOG(logDEBUG1) << "MultiHypotheses::couple_outgoing()";
  MAX_OUTGOING_ARCS max_outgoing_arcs = DIVISION;
  
  // maximum division level constraint by soft constraint from feature?
  // if (s->Level > max_division_level_) {
  //   max_outgoing_arcs = TRANSITION;
  // }
  if (cplex_idxs.size() > 0) {
    std::vector<int> coeffs(cplex_idxs.size(), 1);
    const size_t max_on = has_divisions() ? max_outgoing_arcs : TRANSITION;
    // 0 <= 1*transition + ... + 1*transition <= 2 [div and level <= max_level] or 1 [no div or level > max_level]
    cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, max_on);
  }
}


void ModelBuilder::couple_incoming( const std::vector<size_t>& cplex_idxs, OpengmLPCplex& cplex ) {
  LOG(logDEBUG1) << "MultiHypotheses::couple_incoming() entering";
  if (cplex_idxs.size() > 0) {
    std::vector<int> coeffs(cplex_idxs.size(), 1);
    // 0 <= 1*transition + ... + 1*transition <= 1
    cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 1);
  }
  LOG(logDEBUG1) << "MultiHypotheses::couple_incoming() done";
}


void ModelBuilder::couple( const multihypotheses::Model& m,
                           const MultiHypothesesGraph::Node& n,
                           const MultiHypothesesGraph::Arc& a,
                           OpengmLPCplex& cplex ) {
  std::vector<size_t> cplex_idxs;
  cplex_idxs.push_back(cplex_id(cplex, m.var_of_node(n)));
  cplex_idxs.push_back(cplex_id(cplex, m.var_of_arc(a)));
  std::vector<int> coeffs;
  coeffs.push_back(1);
  coeffs.push_back(-1);
  // 0 <= 1*detection - 1*transition <= 1
  cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 1);
}


void ModelBuilder::couple_conflicts( const multihypotheses::Model& m,
                                     const MultiHypothesesGraph& hypotheses,
                                     OpengmLPCplex& cplex ) {
  LOG(logDEBUG1) << "ModelBuilder::couple_conflicts(): entered";
  const ConflictMap& conflicts = hypotheses.get_conflicts();
  const MultiHypothesesGraph::TraxelMap& traxels = hypotheses.get(node_traxel());
  for (ConflictMap::const_iterator timestep = conflicts.begin(); timestep != conflicts.end(); ++timestep) {
    const std::vector<std::vector<unsigned> >& conflict_vector = timestep->second;
    for (std::vector<std::vector<unsigned> >::const_iterator conflict = conflict_vector.begin();
         conflict != conflict_vector.end();
         ++conflict) {
      std::vector<size_t> cplex_idxs;
      stringstream name;
      name << "Conflict constraint: ";
      for (std::vector<unsigned>::const_iterator id = conflict->begin(); id != conflict->end(); ++id) {
        LOG(logDEBUG4) <<  "ModelBuilder::couple_conflicts(): adding id " << *id << " at time " << timestep->first;
        MultiHypothesesGraph::TraxelMap::ItemIt n(traxels, Traxel(*id, timestep->first));
        cplex_idxs.push_back(cplex_id(cplex, m.var_of_node(n)));
        name << Traxel(*id, timestep->first) << "(" << *(cplex_idxs.rbegin()) <<"),";
      }
      name << "\b";
      couple_conflict(cplex_idxs, cplex, name.str());
    }
  }
  LOG(logDEBUG1) << "ModelBuilder::couple_conflicts(): done";
}


void ModelBuilder::couple_conflict( const std::vector<size_t>& cplex_idxs, OpengmLPCplex& cplex, const std::string& name ) {
  LOG(logDEBUG3) << "ModelBuilder::couple_conflict(): entered";
  if (cplex_idxs.size() > 0) {
    std::vector<int> coeffs(cplex_idxs.size(), 1);
    // 0 <= 1*detection + ... + 1*detection <= 1
    cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 1, name.c_str());
  }
  LOG(logDEBUG3) << "ModelBuilder::couple_conflict(): done";
}


void ModelBuilder::add_count_hard_constraints( OpengmLPCplex& cplex ) const {
  for(std::vector<std::vector<std::pair<std::pair<size_t, size_t>, int> > >::const_iterator constraint_it =
          var_state_coeff_constraints_.begin(); constraint_it != var_state_coeff_constraints_.end();
      ++constraint_it) {
    std::vector<size_t> cplex_idxs;
    std::vector<int> coeffs;
    for (std::vector<std::pair<std::pair<size_t, size_t>, int> >::const_iterator it = constraint_it->begin();
         it != constraint_it->end(); ++it) {
      cplex_idxs.push_back(cplex_id(cplex, it->first.first, it->first.second));
      coeffs.push_back(it->second);
    }
    cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 0);
  }
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



////
//// class CVPR2014ModelBuilder
////

boost::shared_ptr<ModelBuilder> CVPR2014ModelBuilder::clone() const {
  return boost::shared_ptr<ModelBuilder>(new CVPR2014ModelBuilder(*this));
}


boost::shared_ptr<Model> CVPR2014ModelBuilder::build(const MultiHypothesesGraph& hypotheses) {
  LOG(logDEBUG) << "CVPR2014ModelBuilder::build() -- entered";
  
  if( !has_detection_vars() ) {
    throw std::runtime_error("CVPR2014ModelBuilder::build(): option without detection vars not yet implemented");
  }

  if (timestep_range_specified() == false) {
    ModelBuilder::first_timestep_ = hypotheses.earliest_timestep();
    ModelBuilder::last_timestep_ = hypotheses.latest_timestep();
  }

  boost::shared_ptr<Model> model(new Model);

  
  if (has_detection_vars()) {
    add_detection_vars( hypotheses, *model );
  }
  add_assignment_vars( hypotheses, *model );

  MultiHypothesesGraph::node_timestep_map& timesteps = hypotheses.get(node_timestep());
  if ( has_detection_vars() ) {
    if (has_maximal_conflict_cliques() && has_conflict_factors()) {
      add_conflict_factors( hypotheses, *model );
    } else {
      for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
        if (timestep_range_specified() &&
            (timesteps[n] < first_timestep() || timesteps[n] > last_timestep())) {
          continue;
        }
        add_detection_factor( hypotheses, *model, n );
      }
    }
  }

  add_count_factors(hypotheses, *model);

  for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
    if (timestep_range_specified() &&
        (timesteps[n] < first_timestep() || timesteps[n] > last_timestep())) {
      continue;
    }
    add_outgoing_factor( hypotheses, *model, n );
    add_incoming_factor( hypotheses, *model, n );
  }

  return model;
}


void CVPR2014ModelBuilder::add_conflict_factors( const MultiHypothesesGraph& hypotheses, Model& m ) const {
  LOG(logDEBUG3) << "CVPR2014ModelBuilder::add_conflict_factors() -- add factors for conflict sets";
  const ConflictMap& conflict_map = hypotheses.get_conflicts();
  const MultiHypothesesGraph::TraxelMap& traxel_map = hypotheses.get(node_traxel());
  for (ConflictMap::const_iterator timestep = conflict_map.begin(); timestep != conflict_map.end(); ++timestep) {
    if (timestep_range_specified() &&
        (timestep->first < first_timestep() || timestep->first > last_timestep())) {
      continue;
    }
    for (ConflictSetVector::const_iterator conflict = timestep->second.begin();
         conflict != timestep->second.end();
         ++conflict) {
      size_t table_dim = conflict->size();
      std::vector<const Traxel*> traxels;
      std::vector<size_t> vi;
      LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_conflict_factors(): get variable indices (vi) and traxels at time t="
                     << timestep->first;
      for (ConflictSet::const_iterator det = conflict->begin();
           det != conflict->end();
           ++det) {
        MultiHypothesesGraph::TraxelMap::ItemIt n(traxel_map, Traxel(*det, timestep->first));
        traxels.push_back(&traxel_map[n]);
        vi.push_back(m.var_of_node(n));
      }
    
      // create empty table from vi and coordinate vector
      OpengmExplicitFactor<double> table(vi);
      std::vector<size_t> coords(table_dim, 0);

      // set table values for one active detection in conflict set and get maximum energy
      assert(traxels.size() == table_dim && "There should be a traxel for each detection in the conflict.");
    
      double deactivated_energy_max = 0.;
      for (size_t coord_index = 0; coord_index < table_dim; ++coord_index) {
        LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_conflict_factor() -- current deactivated_energy_max=" << deactivated_energy_max
                       << "," << *(traxels[coord_index]);
        coords[coord_index] = 1;
        double deactivated_energy_curr = detection()(*(traxels[coord_index]), 0);
        if (deactivated_energy_curr > deactivated_energy_max) {
          deactivated_energy_max = deactivated_energy_curr;
        }
        table.set_value( coords, detection()(*(traxels[coord_index]), 1));
        coords[coord_index] = 0;
      }
      table.set_value( coords, deactivated_energy_max + opportunity_cost());
      LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_conflict_factor() -- maximum deactivation energy: "
                     << deactivated_energy_max;
    }
  }
}


void CVPR2014ModelBuilder::add_detection_factor( const MultiHypothesesGraph& hypotheses,
                                                 Model& m,
                                                 const MultiHypothesesGraph::Node& node) const {
  LOG(logDEBUG3) << "CVPR2014ModelBuilder::add_detection_factor() -- entered";
  size_t vi[] = {m.var_of_node(node)};
  std::vector<size_t> coords(1, 0);
  OpengmExplicitFactor<double> table(vi, vi+1);

  const Traxel& trax = hypotheses.get(node_traxel())[node];

  coords[0] = 0;
  table.set_value( coords, non_detection()(trax, 0) +opportunity_cost() );

  coords[0] = 1;
  table.set_value( coords, detection()(trax, 1) );

  LOG(logDEBUG3) << "CVPR2014ModelBuilder::add_detection_factor: for "
                 << trax << ": detection=" << table.get_value(std::vector<size_t>(1,1))
                 << ", non_detection=" << table.get_value(std::vector<size_t>(1,0));
  table.add_to( *(m.opengm_model) );

}

void CVPR2014ModelBuilder::add_count_factors( const MultiHypothesesGraph& hypotheses, Model& m) {
  const MultiHypothesesGraph::ConnectedComponentMap& components = hypotheses.get(node_connected_component());
  const MultiHypothesesGraph::TraxelMap& traxel_map = hypotheses.get(node_traxel());
  for (MultiHypothesesGraph::ConnectedComponentMap::ValueIt value = components.beginValue();
       value != components.endValue();
       ++value) {
    std::vector<size_t> vi;
    const Traxel* component_traxel = NULL;
    for (MultiHypothesesGraph::ConnectedComponentMap::ItemIt n(components, *value);
         n != lemon::INVALID;
         ++n) {
      const Traxel* current_traxel = &traxel_map[n];
      if (current_traxel->Id == *value) {
        component_traxel = current_traxel;
      }
      vi.push_back(m.var_of_node(n));
    }
    assert(component_traxel != NULL && "There must be a traxel for the connected component.");
  }
  if (has_hierarchical_counting_factor()) {
    
  } else {
    
  }
}

void CVPR2014ModelBuilder::add_hierarchical_count_factor( Model& m,
                                                          const std::vector<size_t>& vi,
                                                          const Traxel& component_traxel ) {
  LOG(logDEBUG3) << "CVPR2014ModelBuilder::add_hierachical_count_factor: entered";
  assert(has_hierarchical_counting_factor());
  size_t maximum_active_regions = component_traxel.features.find("cardinality")->second[0];
  // maybe use feature_array& and a non-const Traxel& as function argument here
  feature_array probabilities = component_traxel.features.find("count_prediction")->second;
  fill_probabilities(probabilities, maximum_active_regions);
  // add counter variable
  size_t num_states = std::min(vi.size() + 1, maximum_active_regions + 1);
  size_t var_id;
  if (num_states > 2) {
    m.opengm_model->addVariable(num_states);
    var_id = m.opengm_model->numberOfVariables() - 1;
    assert(m.opengm_model->numberOfLabels(var_id) == num_states);
    
    // add helper hard constraints to be added later
    var_state_coeff_constraints_.push_back(std::vector<std::pair<std::pair<size_t, size_t>, int> >());
    std::vector<std::pair<std::pair<size_t, size_t>, int> >& constraint = *(var_state_coeff_constraints_.rbegin());
    std::vector<std::pair<size_t,size_t> > vars;
    for (std::vector<size_t>::const_iterator index = vi.begin(); index != vi.end(); ++index) {
      std::pair<size_t, size_t> var_state = std::pair<size_t, size_t>(*index, 1);
      std::pair<std::pair<size_t, size_t>, int> var_coeff(var_state, 1);
      constraint.push_back(var_coeff);
    }
    for (size_t state = 0; state < num_states; ++state) {
      std::pair<size_t, size_t> var_state = std::pair<size_t, size_t>(var_id, state);
      std::pair<std::pair<size_t, size_t>, int> var_coeff(var_state, -state);
      constraint.push_back(var_coeff);
    }
    
  } else {
    // if there is only two possible states, we have only one region, aka the connected component
    var_id = vi[0];
  }

  // add the count prior
  size_t var[] = {var_id};
  assert(m.opengm_model->numberOfLabels(var[0]) == probabilities.size());
  std::vector<size_t> coords(1, 0);
  OpengmExplicitFactor<double> table(var, var+1, 0, probabilities.size());
  for (size_t state = 0; state < probabilities.size(); ++state) {
    coords[0] = state;
    table.set_value( coords, count()(probabilities[state]) );
  }
  LOG(logDEBUG3) << "CVPR2014ModelBuilder::add_hierarchical_count_factor: "
                 << "helper count factor added";
  table.add_to( *(m.opengm_model) );
}


void CVPR2014ModelBuilder::add_explicit_count_factor( Model& m,
                                                      const std::vector<size_t>& vi,
                                                      const Traxel& component_traxel ) const {
  LOG(logDEBUG) << "CVPR2014ModelBuilder::add_explicit_count_factor: entered";
  assert(!has_hierarchical_counting_factor());
  assert(vi.size() > 0);
  assert(component_traxel.features.find("count_prediction") != component_traxel.features.end());
  // maybe use feature_array& and a non-const Traxel& as function argument here
  feature_array probabilities = component_traxel.features.find("count_prediction")->second;
  size_t table_dim = vi.size();
  size_t maximum_active_regions = component_traxel.features.find("cardinality")->second[0];
  
  assert(probabilities.size() > 0);

  // probabilities must have one state more than maximum_active_regions to
  // include the state for all regions deactivated
  fill_probabilities(probabilities, maximum_active_regions);

    std::vector<size_t> coords(table_dim, 0);
    OpengmExplicitFactor<double> table( vi, forbidden_cost() );

    for (size_t i = 0; i < static_cast<size_t>(std::pow(2., static_cast<int>(table_dim))); ++i) {
      std::vector<size_t> state = DecToBin(i);
      assert(state.size() <= coords.size());
      std::copy(state.begin(), state.end(), coords.begin());
      LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_count_factor: " << coords.size()
                   << ", " << vi.size();

      size_t active_count = std::accumulate(coords.begin(), coords.end(), 0);
      if (active_count <= maximum_active_regions) {
        table.set_value(coords, count()(probabilities[active_count]));
      }
      coords = std::vector<size_t>(table_dim, 0);
    }
    table.add_to( *m.opengm_model );
    LOG(logDEBUG) << "CVPR2014ModelBuilder::add_explicit_count_factor: done";
}


void CVPR2014ModelBuilder::fill_probabilities(feature_array& probabilities, size_t maximum_active_regions) const {
  if (probabilities.size() < maximum_active_regions + 1) {
    // fill the missing states with copies of the last state
    probabilities.insert(probabilities.end(),
                         maximum_active_regions + 1 - probabilities.size(),
                         *(probabilities.rbegin()));
  } else {
    // merge states n, n+1, ..., end into state n and divide by the
    // number of excess states
    // n is the last valid state, i.e. n == maximum_active_regions
    size_t overhead = probabilities.size() - maximum_active_regions;
    assert(probabilities.begin() + maximum_active_regions < probabilities.end());
    probabilities[maximum_active_regions] =
        std::accumulate(probabilities.begin() + maximum_active_regions,
                        probabilities.end(),
                        0.) / overhead;
    probabilities.resize(maximum_active_regions + 1);
  }
  
  assert(probabilities.size() == maximum_active_regions + 1);

  // renormalize
  feature_type sum = std::accumulate(probabilities.begin(), probabilities.end(), 0.);
  if (sum > 0) {
    for(feature_array::iterator p = probabilities.begin(); p != probabilities.end(); ++p) {
      *p /= sum;
    }
  } else {
    feature_type constant = 1/probabilities.size();
    std::fill(probabilities.begin(), probabilities.end(), constant);
  }
}


namespace {
double get_transition_prob(double distance, size_t state, double alpha) {
  double prob = exp(-distance / alpha);
  if (state == 0) {
    return 1 - prob;
  }
  return prob;
}
}


void CVPR2014ModelBuilder::add_outgoing_factor( const MultiHypothesesGraph& hypotheses,
                                                Model& m,
                                                const MultiHypothesesGraph::Node& node) const {
  const MultiHypothesesGraph::TraxelMap& traxel_map = hypotheses.get(node_traxel());
  const MultiHypothesesGraph::ConnectedComponentMap& component_map = hypotheses.get(node_connected_component());
  const Traxel& trax = traxel_map[node];
  const Traxel& connected_component_traxel =
      traxel_map[MultiHypothesesGraph::ConnectedComponentMap::ItemIt(component_map, std::make_pair(trax.Timestep, trax.Component))];
  feature_type maximum_cardinality = connected_component_traxel.features.find("cardinality")->second[0];
  
  const vector<size_t> vi = vars_for_outgoing_factor(hypotheses, m, node);
  LOG(logDEBUG3) << "CVPR2014ModelBuilder::add_outgoing_factor(): entered for " << trax
                 << " at " << trax.features.find("com")->second[0] << "," << trax.features.find("com")->second[1] << ","
                 << trax.features.find("com")->second[2] << " lvl: " << trax.Level << " - factor order: " << vi.size();
  if (vi.size() == 0) {
    // nothing to do here
    // happens in case of no det vars and no outgoing arcs
    return;
  }

  // collect arcs for use in feature functions
  // and target traxels for potential transitions
  std::vector<Model::arc_t> arcs;
  std::vector<const Traxel*> target_traxels;
  for (std::vector<size_t>::const_iterator v = vi.begin()+1;
       v != vi.end();
       ++v) {
    arcs.push_back(m.arc_of_var(*v));
    target_traxels.push_back(&traxel_map[hypotheses.target(*(arcs.rbegin()))]);
  }

  size_t table_dim = vi.size();
  assert(table_dim == target_traxels.size() + 1);
  LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_outgoing_factor(): table_dim=" << table_dim;

    
  std::vector<size_t> coords(table_dim, 0);
  OpengmExplicitFactor<double> table( vi, forbidden_cost() );
  
  // opportunity?
  table.set_value( coords, 0 );
  
  
  LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_outgoing_factor(): initializing minimum energies";

  // if no arcs, then initialize energy to zero
  double maximum_non_move_energy = 0.;
  double maximum_non_division_energy = 0.;

  
  // move configuration
  coords[0] = 1;
  for (size_t i = 1; i < table_dim; ++i) {
    coords[i] = 1;
    feature_type probability = 0.;
    if (has_classifiers() && transition_parameter() == 0) {
      probability = hypotheses.get(node_move_features())[node]
         .find(target_traxels[i-1]->Id)->second[1];
      LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_outgoing_factor: move using classifier, prob: " << probability
                     << ", energy: " << move()(trax, *(target_traxels[i-1]), probability);
    } else if (transition_parameter() != 0) {
      probability = (double) get_transition_prob(trax.distance_to(*(target_traxels[i-1])), /* state */ 1, transition_parameter());
    }


    double move_energy = move()(trax, *(target_traxels[i-1]), probability);
    double non_move_energy = move()(trax, *(target_traxels[i-1]), 1.-probability);

    if (non_move_energy > maximum_non_move_energy) {
      maximum_non_move_energy = non_move_energy;
    }
    LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_outgoing_factor: move energy: " << move_energy;
    table.set_value( coords,  (trax.features.find("cardinality")->second[0]/maximum_cardinality)*move_energy);
    LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_outgoing_factor: probability = "
                   << probability << ", move=" << table.get_value( coords );
    LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_outgoing_factor: cardinality =" << trax.features.find("cardinality")->second[0]
                   << ", maximum_cardinality=" <<  maximum_cardinality;
    coords[i] = 0;
  }
  coords[0] = 0;

  // division configuration
  if (has_divisions()) {
    coords[0] = 1;
    for (unsigned int i = 1; i < table_dim - 1; ++i) {
      coords[i] = 1;
      for (unsigned int j = i + 1; j < table_dim; ++j) {
        coords[j] = 1;
        feature_type probability = 0.;
        if (has_classifiers()) {
          probability = hypotheses.get(node_division_features())[node]
              .find(std::make_pair(std::min(target_traxels[i-1]->Id, target_traxels[j-1]->Id),
                                   std::max(target_traxels[i-1]->Id, target_traxels[j-1]->Id)))->second[1];
        }
        double division_energy = division()(trax, *(target_traxels[i-1]), *(target_traxels[j-1]), probability);
        double non_division_energy = division()(trax, *(target_traxels[i-1]), *(target_traxels[j-1]), 1.-probability);
        if (non_division_energy > maximum_non_division_energy) {
          maximum_non_division_energy = non_division_energy;
        }
        assert(trax.features.find("cardinality") != trax.features.end() && "Cardinality must be present");
        table.set_value(coords, (trax.features.find("cardinality")->second[0]/maximum_cardinality)*division_energy);
        LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_outgoing_factor: division="
                       << table.get_value( coords );
        coords[j] = 0;
      }
      coords[i] = 0;
    }
    coords[0] = 0;
  }

  // disappearance configuration
  if (trax.Timestep < hypotheses.latest_timestep()) {
    coords[0] = 1;
    table.set_value( coords,
                     (trax.features.find("cardinality")->second[0]/maximum_cardinality)*(disappearance()(trax)) +
                     std::max(maximum_non_move_energy, maximum_non_division_energy)
                     );
    LOG(logDEBUG3) << "CVPR2014ModelBuilder::add_outgoing_factor: at least two outgoing arcs: "
                   << "forbidden=" << forbidden_cost() << ", disappearance=" << table.get_value(coords);
    coords[0] = 0;
  }

    
  table.add_to( *m.opengm_model );

  LOG(logDEBUG2) << "CVPR2014ModelBuilder::add_outgoing_factor(): leaving";
}


void CVPR2014ModelBuilder::add_incoming_factor( const MultiHypothesesGraph& hypotheses,
                                                Model& m,
                                                const MultiHypothesesGraph::Node& node ) const {
  const MultiHypothesesGraph::TraxelMap& traxel_map = hypotheses.get(node_traxel());
  const MultiHypothesesGraph::ConnectedComponentMap& component_map = hypotheses.get(node_connected_component());
  const Traxel& trax = traxel_map[node];
  const Traxel& connected_component_traxel =
      traxel_map[MultiHypothesesGraph::ConnectedComponentMap::ItemIt(component_map, std::make_pair(trax.Timestep, trax.Component))];
  feature_type maximum_cardinality = connected_component_traxel.features.find("cardinality")->second[0];
  LOG(logDEBUG2) << "CVPR2014ModelBuilder::add_incoming_factor(): entered for " << trax;
  const std::vector<size_t> vi = vars_for_incoming_factor(hypotheses, m, node);
  if (vi.size() == 0) {
    // nothing to be done here
    // no det vars and no incoming arcs
    return;
  }

  // construct factor
  LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_incoming_factor: " << vi.size();

  
  LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_incoming_factor: constructing factor for "
                 << trax << ": " << std::pow(2., static_cast<int>(vi.size())) << " entries (2^" << vi.size() << ")";
  const size_t table_dim = vi.size();
  std::vector<size_t> coords(table_dim, 0);
  OpengmExplicitFactor<double> table( vi, forbidden_cost() );

  // appearance
  if (trax.Timestep > hypotheses.earliest_timestep()) {
    coords[0] = 1;
    table.set_value( coords, (trax.features.find("cardinality")->second[0]/maximum_cardinality)*appearance()(trax) );
    // assert(table.get_value( coords ) == (trax.features.find("cardinality")->second[0]/maximum_cardinality+1.5)*appearance()(trax) );
    LOG(logDEBUG2) << "CVPR2014ModelBuilder::add_incoming_factor: appearance="
                   << table.get_value( coords );
    coords[0] = 0;
  }

  table.add_to( *m.opengm_model );
  LOG(logDEBUG2) << "CVPR2014ModelBuilder::add_incoming_factor: done";
}


    

} // namespace multihypotheses
} // namespace pgm
} // namespace pgmlink 
