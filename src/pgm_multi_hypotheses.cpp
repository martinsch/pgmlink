// stl
#include <algorithm>
#include <cassert>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>
#include <iterator>
#include <functional>

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
             const trax_var_map& trax_var,
             const arc_var_map& arc_var
             )
    : opengm_model(m) {
  trax_var_.left = trax_var;
  arc_var_.left = arc_var;
  init();
}

const Model::trax_var_map& Model::var_of_trax() const {
  return trax_var_.left;
}

const Model::var_trax_map& Model::trax_of_var() const {
  return trax_var_.right;
}

const Model::arc_var_map& Model::var_of_arc() const {
  return arc_var_.left;
}

const Model::var_arc_map& Model::arc_of_var() const
{
  return arc_var_.right;      
}

Model::var_t Model::var_of_trax(const Traxel& e) const {
  trax_var_map::const_iterator it = var_of_trax().find(e);
  if(it!=var_of_trax().end()) {
    return it->second;
  } else {
    throw std::out_of_range("chaingraph::Model::var_of_trax(): key does not exist");
  }
}

Model::var_t Model::var_of_arc(const TraxelArc& e) const {
  arc_var_map::const_iterator it = var_of_arc().find(e);
  if(it!=var_of_arc().end()) {
    return it->second;
  } else {
    throw std::out_of_range("MultiHypothesesModel::var_of_arc(): key does not exist");
  }
}

Traxel Model::trax_of_var(var_t e) const {
  var_trax_map::const_iterator it = trax_of_var().find(e);
  if(it!=trax_of_var().end()) {
    return it->second;
  } else {
    throw std::out_of_range("MultiHypothesesModel::trax_of_var(): key does not exist");
  }
}

Model::TraxelArc Model::arc_of_var(var_t e) const {
  var_arc_map::const_iterator it = arc_of_var().find(e);
  if(it!=arc_of_var().end()) {
    return it->second;
  } else {
    throw std::out_of_range("ChaingraphModel::arc_of_var(): key does not exist");
  }
}

Model::VarCategory Model::var_category(var_t e) const {
  if(arc_of_var().count(e)) {
    return Model::arc_var;
  } else if(trax_of_var().count(e)) {
    return Model::trax_var;
  } else {
    throw std::out_of_range("ChaingraphModel::var_category(): key does not exist");
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
}


size_t ModelBuilder::cplex_id(OpengmLPCplex& cplex, const size_t opengm_id, const size_t state) const {
    return cplex.lpNodeVi(opengm_id, state);
}

size_t ModelBuilder::cplex_id(OpengmLPCplex& cplex, const size_t opengm_id) const {
    return cplex_id(cplex, opengm_id, 1);
}


void ModelBuilder::add_hard_constraints(const Model& m, const MultiHypothesesGraph& hypotheses, OpengmLPCplex& cplex) {
  LOG(logDEBUG) << "MultiHypotheses::add_hard_constraints: entered";
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  MultiHypothesesGraph::node_timestep_map& timesteps = hypotheses.get(node_timestep());

  
  for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
    if (timestep_range_specified() &&
        (timesteps[n] < first_timestep() || timesteps[n] > last_timestep())) {
      continue;
    }
    LOG(logDEBUG1) << "MultiHypotheses::add_hard_constraints: outgoing transitions";
    const std::vector<Traxel>& traxels = regions[n];
    std::vector<Traxel> traxels_dest;
    for (MultiHypothesesGraph::OutArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
      const std::vector<Traxel>& traxels_at = regions[hypotheses.target(a)];
      traxels_dest.insert(traxels_dest.end(), traxels_at.begin(), traxels_at.end());      
    }
    couple_outgoing(m, traxels, traxels_dest, cplex);

    LOG(logDEBUG1) << "MultiHypotheses::add_hard_constraints: incoming transitions";
    std::vector<Traxel> traxels_src;
    for (MultiHypothesesGraph::InArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
      const std::vector<Traxel>& traxels_at = regions[hypotheses.source(a)];
      traxels_src.insert(traxels_src.end(), traxels_at.begin(), traxels_at.end());
    }
    couple_incoming(m, traxels_src, traxels, cplex);

    if (has_detection_vars()) {
      LOG(logDEBUG1) << "MultiHypotheses::add_hard_constraints: coupling conflicts";
      couple_conflicts(m, hypotheses, n, traxels, cplex);
      
      // LOG(logDEBUG1) << "MultiHypotheses::add_hard_constraints: coupling count";
      // couple_count(m, traxels, cplex);
    }
  }
  LOG(logDEBUG) << "MultiHypotheses::add_hard_constraints: exited";
}


void ModelBuilder::set_cplex_timeout( double seconds ) {
  cplex_timeout_ = seconds;
}


inline void ModelBuilder::add_detection_vars( const MultiHypothesesGraph& hypotheses, Model& m ) const {
  if(!has_detection_vars()) {
    throw std::runtime_error("multihypotheses::ModelBuilder::add_detection_vars(): called without has_detection_vars()");
  }
  LOG(logDEBUG) << "ModelBuilder::add_detection_vars: entered";
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  MultiHypothesesGraph::node_timestep_map& timesteps = hypotheses.get(node_timestep());
  for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
    std::vector<Traxel>& traxels = regions.get_value(n);
    traxels[0].Level = 1;
    const MultiHypothesesGraph::node_timestep_map::Value& timestep = timesteps[n];
    for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
      if (timestep_range_specified() &&
          (timestep < first_timestep() || timestep > last_timestep())) {
        continue;
      }
      m.opengm_model->addVariable(2);
      m.trax_var_.left.insert(Model::trax_var_map::value_type(*t, m.opengm_model->numberOfVariables() - 1));
      LOG(logDEBUG4) << timestep << ',' << first_timestep() << ',' << last_timestep();
      LOG(logDEBUG4) << "ModelBuilder::add_detection_vars: added var " << m.opengm_model->numberOfVariables() - 1
                     << " for " << m.trax_of_var(m.opengm_model->numberOfVariables() - 1);
    }
  }
}


inline void ModelBuilder::add_assignment_vars( const MultiHypothesesGraph& hypotheses, Model& m ) const {
  LOG(logINFO) << "add_assignment_vars() -- all neighbors";
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  MultiHypothesesGraph::node_timestep_map& timesteps = hypotheses.get(node_timestep());
  for (MultiHypothesesGraph::ArcIt a(hypotheses); a != lemon::INVALID; ++a) {
    const MultiHypothesesGraph::Node& source_node = hypotheses.source(a);
    const MultiHypothesesGraph::Node& dest_node = hypotheses.target(a);
    if (timestep_range_specified() &&
        (timesteps[source_node] < first_timestep() || timesteps[dest_node] > last_timestep())) {
      continue;
    }
    const std::vector<Traxel>& source = regions[source_node];
    const std::vector<Traxel>& dest = regions[dest_node];
    for (std::vector<Traxel>::const_iterator s = source.begin(); s != source.end(); ++s) {
      for (std::vector<Traxel>::const_iterator d = dest.begin(); d != dest.end(); ++d) {
        m.opengm_model->addVariable(2);
        m.arc_var_.left.insert(Model::arc_var_map::value_type(Model::TraxelArc(*s, *d), m.opengm_model->numberOfVariables() - 1));
      }
    }
  }
}


namespace {
class TraxelProbabilityPairLessThan {
 public:
  bool operator()(const std::pair<Traxel, feature_type>& p1, const std::pair<Traxel, feature_type>& p2) {
    return p1.second < p2.second;
  }
};


class TraxelProbabilityPairGreaterThan {
 public:
  bool operator()(const std::pair<Traxel, feature_type>& p1, const std::pair<Traxel, feature_type>& p2) {
    return p1.second > p2.second;
  }
};
} /* namespace */


inline void ModelBuilder::add_assignment_vars( const MultiHypothesesGraph& hypotheses, Model& m, const MultiHypothesesGraph::MoveFeatureMap& moves ) const {
  LOG(logINFO) << "add_assignment_vars() -- " << maximum_outgoing_arcs_ << " best neighbors";
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  MultiHypothesesGraph::node_timestep_map& timesteps = hypotheses.get(node_timestep());
  for (MultiHypothesesGraph::ArcIt a(hypotheses); a != lemon::INVALID; ++a) {
    const MultiHypothesesGraph::Node& source_node = hypotheses.source(a);
    const MultiHypothesesGraph::Node& dest_node = hypotheses.target(a);
    LOG(logDEBUG4) << "add_assignment_vars() -- " << timesteps[source_node] << ',' << first_timestep() << ',' << timesteps[dest_node] << ',' << last_timestep();
    if (timestep_range_specified() &&
        (timesteps[source_node] < first_timestep() || timesteps[dest_node] > last_timestep())) {
      continue;
    }
    const std::vector<Traxel>& source = regions[source_node];
    const std::vector<Traxel>& dest = regions[dest_node];
    for (std::vector<Traxel>::const_iterator s = source.begin(); s != source.end(); ++s) {
      const std::map<Traxel, feature_array>& probabilities = moves[hypotheses.source(a)].find(*s)->second;
      std::vector<std::pair<Traxel, feature_type> > k_best_probabilities;
      for (std::map<Traxel, feature_array>::const_iterator it = probabilities.begin();
           it != probabilities.end();
           ++it) {
        k_best_probabilities.push_back(std::make_pair(it->first, it->second[1]));
        std::sort(k_best_probabilities.begin(), k_best_probabilities.end(), TraxelProbabilityPairGreaterThan());
        k_best_probabilities.resize(std::min(k_best_probabilities.size(), static_cast<size_t>(maximum_outgoing_arcs_)));
      }
      for (std::vector<Traxel>::const_iterator d = dest.begin(); d != dest.end(); ++d) {
        if (std::find(k_best_probabilities.begin(), k_best_probabilities.end(), *d) != k_best_probabilities.end()) {
          m.opengm_model->addVariable(2);
          m.arc_var_.left.insert(Model::arc_var_map::value_type(Model::TraxelArc(*s, *d), m.opengm_model->numberOfVariables() - 1));
        }
      }
    }
  }
  LOG(logINFO) << "add_assignment_vars() -- finished";
}


std::vector<OpengmModel::IndexType> ModelBuilder::vars_for_outgoing_factor( const MultiHypothesesGraph& hypotheses,
                                                                            const Model& m,
                                                                            const MultiHypothesesGraph::Node& node,
                                                                            const Traxel& trax) const {
  LOG(logDEBUG2) << "ModelBuilder::vars_for_outgoing_factor() -- entered";
  std::vector<OpengmModel::IndexType> vi; // opengm variable indices; may be empty if no det vars
  if (has_detection_vars()) {
    vi.push_back(m.var_of_trax(trax)); // first detection var, the others will be assignment vars
  }
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  for (MultiHypothesesGraph::OutArcIt a(hypotheses, node); a != lemon::INVALID; ++a) {
    const std::vector<Traxel>& traxels = regions[hypotheses.target(a)];
    for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
      Model::arc_var_map::const_iterator it = m.var_of_arc().find(Model::TraxelArc(trax, *t));
      if (it != m.var_of_arc().end()) {
          vi.push_back(it->second);
      }
      // vi.push_back(m.var_of_arc(Model::TraxelArc(trax, *t)));
    }
  }
  return vi;
}


std::vector<OpengmModel::IndexType> ModelBuilder::vars_for_incoming_factor( const MultiHypothesesGraph& hypotheses,
                                                                            const Model& m,
                                                                            const MultiHypothesesGraph::Node& node,
                                                                            const Traxel& trax) const {
  LOG(logDEBUG2) << "ModelBuilder::vars_for_incoming_factor() -- entered";
  std::vector<OpengmModel::IndexType> vi; // opengm variable indices; may be empty if no det vars
  if (has_detection_vars()) {
    vi.push_back(m.var_of_trax(trax)); // first detection var, the others will be assignment vars
  }
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  for (MultiHypothesesGraph::InArcIt a(hypotheses, node); a != lemon::INVALID; ++a) {
    const std::vector<Traxel>& traxels = regions[hypotheses.source(a)];
    for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
      Model::arc_var_map::const_iterator it = m.var_of_arc().find(Model::TraxelArc(*t, trax));
      if (it != m.var_of_arc().end()) {
          vi.push_back(it->second);
      }
      // vi.push_back(m.var_of_arc(Model::TraxelArc(*t, trax)));
    }
  }

  
  // std::reverse(vi.begin(), vi.end()); // det var should be the first index to be consistent with vars_for_outgoing_factor() WHY? -> ASK BERNHARD!!
  // it seems to be the last one here
  return vi;
}


void ModelBuilder::couple_outgoing(const Model& m, const std::vector<Traxel>& source, const std::vector<Traxel>& dest, OpengmLPCplex& cplex) {
  LOG(logDEBUG1) << "MultiHypotheses::couple_outgoing()";
  if (has_detection_vars()) {
    couple_detections_assignments(m, source, dest, cplex);
  }
  couple_outgoing_assignments(m, source, dest, cplex);
}


void ModelBuilder::couple_incoming(const Model& m, const std::vector<Traxel>& source, const std::vector<Traxel>& dest, OpengmLPCplex& cplex) {
  LOG(logDEBUG1) << "MultiHypotheses::couple_incoming()";
  if (has_detection_vars()) {
    couple_detections_assignments_incoming(m, source, dest, cplex);
  }
  couple_incoming_assignments(m, source, dest, cplex);
}


void ModelBuilder::couple_detections_assignments(const Model& m, const std::vector<Traxel>& source, const std::vector<Traxel>& dest, OpengmLPCplex& cplex) {
  // LOG(logDEBUG1) << "MultiHypotheses::couple_detection_assignments()";
  LOG(logDEBUG2) << "MultiHypotheses::couple_detection_assignments: "
                 << source.size() << " source(s) and "
                 << dest.size() << " target(s)";
  for (std::vector<Traxel>::const_iterator s = source.begin(); s != source.end(); ++s) {
    for (std::vector<Traxel>::const_iterator d = dest.begin(); d != dest.end(); ++d) { 
      Model::arc_var_map::const_iterator it = m.var_of_arc().find(Model::TraxelArc(*s, *d));
      if (it != m.var_of_arc().end()) {
        std::vector<size_t> cplex_idxs;
        LOG(logDEBUG4) << "MultiHypotheses::couple_detecion_assignments: "
                       << *s << "," << *d;
        cplex_idxs.push_back(cplex_id(cplex, m.var_of_trax(*s)));
        cplex_idxs.push_back(cplex_id(cplex, it->second));
        std::vector<int> coeffs;
        coeffs.push_back(1);
        coeffs.push_back(-1);
        // 0 <= 1*detection - 1*transition <= 1
        cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 1);
      }
    }
  }
}

void ModelBuilder::couple_detections_assignments_incoming(const Model& m, const std::vector<Traxel>& source, const std::vector<Traxel>& dest, OpengmLPCplex& cplex) {
  LOG(logDEBUG1) << "MultiHypotheses::couple_detection_assignments_incoming()";
  LOG(logDEBUG2) << "MultiHypotheses::couple_detection_assignments_incoming: "
                 << source.size() << " source(s) and "
                 << dest.size() << " target(s)";
  for (std::vector<Traxel>::const_iterator d = dest.begin(); d != dest.end(); ++d) {
    for (std::vector<Traxel>::const_iterator s = source.begin(); s != source.end(); ++s) {
      Model::arc_var_map::const_iterator it = m.var_of_arc().find(Model::TraxelArc(*s, *d));
      if (it != m.var_of_arc().end()) {
        std::vector<size_t> cplex_idxs;
        LOG(logDEBUG4) << "MultiHypotheses::couple_detecion_assignments_incoming: "
                       << *s << "," << *d;
        cplex_idxs.push_back(cplex_id(cplex, m.var_of_trax(*d)));
        cplex_idxs.push_back(cplex_id(cplex, it->second));
        std::vector<int> coeffs;
        coeffs.push_back(1);
        coeffs.push_back(-1);
        // 0 <= 1*detection - 1*transition <= 1
        cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 1);
      }
    }
  }
}


void ModelBuilder::couple_outgoing_assignments(const Model& m, const std::vector<Traxel>& source, const std::vector<Traxel>& dest, OpengmLPCplex& cplex) {
  LOG(logDEBUG1) << "MultiHypotheses::couple_outgoing_assignments()";
  for (std::vector<Traxel>::const_iterator s = source.begin(); s != source.end(); ++s) {
    MAX_OUTGOING_ARCS max_outgoing_arcs = DIVISION;
    //if (s->Level > max_division_level_) {
    //  max_outgoing_arcs = TRANSITION;
    //}
    std::vector<size_t> cplex_idxs;
    for (std::vector<Traxel>::const_iterator d = dest.begin(); d != dest.end(); ++d) {
      Model::arc_var_map::const_iterator it = m.var_of_arc().find(Model::TraxelArc(*s, *d));
      if (it != m.var_of_arc().end()) {
        cplex_idxs.push_back(cplex_id(cplex, it->second));
      }
    }
    if (cplex_idxs.size() > 0) {
      std::vector<int> coeffs(cplex_idxs.size(), 1);
      const size_t max_on = has_divisions() ? max_outgoing_arcs : TRANSITION;
      // 0 <= 1*transition + ... + 1*transition <= 2 [div and level <= max_level] or 1 [no div or level > max_level]
      cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, max_on);
    }
  }
}


void ModelBuilder::couple_incoming_assignments(const Model& m, const std::vector<Traxel>& source, const std::vector<Traxel>& dest, OpengmLPCplex& cplex) {
  LOG(logDEBUG1) << "MultiHypotheses::couple_incoming_assignments()";
  for (std::vector<Traxel>::const_iterator d = dest.begin(); d != dest.end(); ++d) {
    std::vector<size_t> cplex_idxs;
    for (std::vector<Traxel>::const_iterator s = source.begin(); s != source.end(); ++s) {
      Model::arc_var_map::const_iterator it = m.var_of_arc().find(Model::TraxelArc(*s, *d));
      if (it != m.var_of_arc().end()) {
        cplex_idxs.push_back(cplex_id(cplex, it->second));
      }
    }
    if (cplex_idxs.size() > 0) {
      std::vector<int> coeffs(cplex_idxs.size(), 1);
      // 0 <= 1*transition + ... + 1*transition <= 1
      cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 1);
    }
  }
}


void ModelBuilder::couple_count( const Model& m, const std::vector<Traxel>& traxels, OpengmLPCplex& cplex) {
  LOG(logDEBUG1) << "MultiHypotheses::couple_count()";
  std::vector<size_t> cplex_idxs;
  for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
    cplex_idxs.push_back(cplex_id(cplex, m.var_of_trax(*t)));
  }
  if (cplex_idxs.size() > 0) {
    std::vector<int> coeffs(cplex_idxs.size(), 1);
    // 0 <= 1*detection + ... + 1*detection <= max_count_
    cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, max_count_);
  }
}


void ModelBuilder::couple_conflicts( const Model& m,
                                     const MultiHypothesesGraph& hypotheses,
                                     const MultiHypothesesGraph::Node& node,
                                     const std::vector<Traxel>& traxels,
                                     OpengmLPCplex& cplex) {
  LOG(logDEBUG4) << "MultiHypotheses::couple_conflicts()";
  if (has_maximal_conflict_cliques()) {
    couple_conflicts_maximal_cliques( m,
                                      traxels[0].Timestep,
                                      hypotheses.get(node_conflict_sets())[node],
                                      cplex );
  }
  else {
    couple_conflicts_pairwise( m, traxels, cplex);
  }
}


void ModelBuilder::couple_conflicts_pairwise( const Model& m,
                                              const std::vector<Traxel>& traxels,
                                              OpengmLPCplex& cplex) {
  for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
    FeatureMap::const_iterator feature = t->features.find("conflicts");
    if (feature == t->features.end()) {
      throw std::runtime_error("MultiHypotheses: couple_conflicts_pairwise - Feature conflicts not found in traxel!");
    }
    for (feature_array::const_iterator conflict = feature->second.begin(); conflict != feature->second.end(); ++conflict) {
      std::vector<size_t> cplex_idxs;
      cplex_idxs.push_back(cplex_id(cplex, m.var_of_trax(*t)));
      LOG(logDEBUG4) << "MultiHypotheses::couple_conflicts_pairwise: Adding cplex ids for " << *t
                     << " conflicting with " << *conflict;
      assert(std::find(traxels.begin(), traxels.end(), Traxel(*conflict, t->Timestep)) != traxels.end());
      cplex_idxs.push_back(cplex_id(cplex, m.var_of_trax(Traxel(*conflict, t->Timestep))));
    
      
      std::vector<int> coeffs(cplex_idxs.size(), 1);
      // 0 <= 1*detection + 1*detection <= 1
      cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), int(has_one_active_per_component_constraint()), 1);
      
    }
  }  
}


void ModelBuilder::couple_conflicts_maximal_cliques( const Model& m,
                                                     int timestep,
                                                     const std::vector<std::vector<unsigned> >& conflict_sets,
                                                     OpengmLPCplex& cplex) {
  LOG(logDEBUG) << "MultiHypotheses::couple_conflicts_maximal_cliques: entered";
  for (std::vector<std::vector<unsigned> >::const_iterator set = conflict_sets.begin();
       set != conflict_sets.end();
       ++set) {
    if (set->size() <= 1) {
      continue;
    }
    std::vector<size_t> cplex_idxs;
    for (std::vector<unsigned>::const_iterator object_id = set->begin();
         object_id != set->end();
         ++object_id) {
      cplex_idxs.push_back(cplex_id(cplex, m.var_of_trax(Traxel(*object_id, timestep))));
    }
    std::vector<int> coeffs(cplex_idxs.size(), 1);
    // 0 <= 1*detection + ... + 1*detection <= 1
    cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 1);
  }
}


////
//// class TrainableModelBuilder
////
boost::shared_ptr<ModelBuilder> TrainableModelBuilder::clone() const {
  return boost::shared_ptr<ModelBuilder>(new TrainableModelBuilder(*this));
}


boost::shared_ptr<Model> TrainableModelBuilder::build(const MultiHypothesesGraph& hypotheses) {

  if ( timestep_range_specified() ) {
    throw std::runtime_error("TrainableModelBuilder not compatible with restricted timestep range!");
  }

  boost::shared_ptr<Model> model(new Model);
  assert(model->opengm_model->numberOfWeights() == 0);
  model->opengm_model->increaseNumberOfWeights(3);
  model->weight_map[Model::mov_weight].push_back(0);
  model->weight_map[Model::app_weight].push_back(1);
  model->weight_map[Model::dis_weight].push_back(2);

  if (has_divisions()) {
    model->opengm_model->increaseNumberOfWeights(1);
    model->weight_map[Model::div_weight].push_back(3);
  }

  if (has_detection_vars()) {
    model->opengm_model->increaseNumberOfWeights(1);
    model->weight_map[Model::det_weight].push_back(4);
  }

  if (has_detection_vars()) {
    add_detection_vars( hypotheses, *model );
  }
  add_assignment_vars( hypotheses, *model );

  if ( has_detection_vars() ) {
    for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
      add_detection_factors( hypotheses, *model, n );
    }
  }

  for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
    add_outgoing_factors( hypotheses, *model, n );
    add_incoming_factors( hypotheses, *model, n );
  }

  return model;
}


void TrainableModelBuilder::add_detection_factors( const MultiHypothesesGraph& hypotheses,
                                                   Model& m,
                                                   const MultiHypothesesGraph::Node& n) const {
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  const std::vector<Traxel>& traxels = regions[n];
  for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
    add_detection_factor(m, *t);
  }
}


void TrainableModelBuilder::add_outgoing_factors( const MultiHypothesesGraph& hypotheses,
                                                  Model& m,
                                                  const MultiHypothesesGraph::Node& n ) const {
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  const std::vector<Traxel>& traxels = regions[n];
  for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
    std::vector<Traxel> neighbors;
    for (MultiHypothesesGraph::OutArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
      const std::vector<Traxel>& neighbors_at = regions[hypotheses.target(a)];
      neighbors.insert(neighbors.end(), neighbors_at.begin(), neighbors_at.end());      
    }
    add_outgoing_factor(hypotheses, m, n, *t, neighbors);
  }

}


void TrainableModelBuilder::add_incoming_factors( const MultiHypothesesGraph& hypotheses,
                                                   Model& m,
                                                   const MultiHypothesesGraph::Node& n ) const {
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  const std::vector<Traxel>& traxels = regions[n];
  for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
    std::vector<Traxel> neighbors;
    for (MultiHypothesesGraph::InArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
      const std::vector<Traxel>& neighbors_at = regions[hypotheses.target(a)];
      neighbors.insert(neighbors.end(), neighbors_at.begin(), neighbors_at.end());
    }
    add_incoming_factor(hypotheses, m, n, *t, neighbors);
  }

}


void TrainableModelBuilder::add_detection_factor( Model& m,
                                                  const Traxel& trax) const {
  
  std::vector<size_t> var_indices;
  var_indices.push_back(m.var_of_trax(trax));
  size_t shape[] = {2};
  size_t indicate[] = {0};

  LOG(logDEBUG2) << "TrainableModelBuilder::add_detection_factor(): entered for " << trax
                 << " - non_detection: " << non_detection()(trax, 0)
                 << " detection: " << detection()(trax, 1);
  OpengmWeightedFeature<OpengmModel::ValueType>(var_indices, shape, shape+1, indicate, non_detection()(trax, 0))
      .add_as_feature_to( *(m.opengm_model), m.weight_map[Model::det_weight].front() );

  indicate[0] = 1;
  OpengmWeightedFeature<OpengmModel::ValueType>(var_indices, shape, shape+1, indicate, detection()(trax, 1))
      .add_as_feature_to( *(m.opengm_model), m.weight_map[Model::det_weight].front() );
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


void TrainableModelBuilder::add_outgoing_factor(const MultiHypothesesGraph& hypotheses,
                                                Model& m,
                                                const MultiHypothesesGraph::Node& node,
                                                const Traxel& trax,
                                                const std::vector<Traxel>& neighbors) const {
  LOG(logDEBUG2) << "TrainableModelBuilder::add_outgoing_factor(): entered for " << trax;
  const vector<size_t> vi = vars_for_outgoing_factor(hypotheses, m, node, trax);
  if (vi.size() == 0) {
    // nothing to do here
    // happens in case of no det vars and no outgoing arcs
    return;
  }

  // collect TraxelArcs for use in feature functions
  std::vector<Model::TraxelArc> arcs;
  for (std::vector<Traxel>::const_iterator neighbor = neighbors.begin();
       neighbor != neighbors.end();
       ++neighbor) {
    arcs.push_back(Model::TraxelArc(trax, *neighbor));
  }

  // construct factor
  const size_t table_dim = vi.size();
  assert(table_dim > 0);
  const std::vector<size_t> shape(table_dim, 2);
  std::vector<size_t> coords;
  

  

  std::set<size_t> entries;
  for (size_t i = 0; i < static_cast<size_t>(std::pow(2., static_cast<int>(table_dim))); ++i) {
    entries.insert(entries.end(), i);
  }

  // #if #FILELOG_MAX_LEVEL == pgmlink::DEBUG4
  /* if (FILELOG_MAX_LEVEL >= pgmlink::logDEBUG4) {
    std::ostringstream os;
    std::ostream_iterator<size_t> ostream_it(os, ", ");
    std::copy(entries.begin(), entries.end(), ostream_it);
    LOG(logDEBUG4) << "TrainableModelBuilder::add_outgoing_factor(): coords size="
                   << entries.size() << " and contents: " << os.str();
  } */
  // #endif /* logDEBUG4 */

  // opportunity configuration??

  // disappearance configuration
  {
    coords = std::vector<size_t>(table_dim, 0);
    if(has_detection_vars()) {
      coords[0] = 1;
    }
    size_t check = entries.erase(BinToDec(coords));
    assert (check == 1);
    double weight = disappearance()(trax);
    if (trax.Timestep == hypotheses.latest_timestep()) {
      weight = 0.;
    }
    LOG(logDEBUG2) << "TrainableModelBuilder::add_outgoing_factor: disappearance cost for "
                   << trax << ": " << weight
                   << " " << disappearance()(trax);
    OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape.begin(), shape.end(), coords.begin(), weight)
        .add_as_feature_to( *(m.opengm_model), m.weight_map[Model::dis_weight].front() );
  }

  // move configuration
  {
    coords = std::vector<size_t>(table_dim, 0);
    size_t assignment_begin = 0;
    if (has_detection_vars()) {
      coords[0] = 1;
      assignment_begin = 1;
    }
    // (1, 0, 0, ..., 1, ..., 0)
  
    for (size_t i = assignment_begin; i < table_dim; ++i) {
      coords[i] = 1;
      // LOG(logDEBUG4) << "TrainableModelBuilder::add_outgoing_factor(): moves: "
      // << "entry to be erased: " << BinToDec(coords);
      size_t check = entries.erase(BinToDec(coords));
      // LOG(logDEBUG4) << "TrainableModelBuilder::add_outgoing_factor(): moves: "
      // << "successfully erased? " << check;
      assert(check == 1);
      feature_type probability  = 0.;
      if (has_classifiers()) {
      // FIXME: NEEDS TO BE IMPLENENTED
        throw std::runtime_error("TrainableModelBuilder does not support classifier priors yet!");
      }
      OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape.begin(), shape.end(), coords.begin(), move()(trax, arcs[i - assignment_begin].second, probability) )
          .add_as_feature_to( *(m.opengm_model), m.weight_map[Model::mov_weight].front() );
      coords[i] = 0;
    }
  }


  // division configurations
  if (has_divisions()) {
    coords = std::vector<size_t>(table_dim, 0);
    size_t assignment_begin = 0;
    if(has_detection_vars()) {
      coords[0] = 1;
      assignment_begin = 1;
    }
    

    // (1, 0, 0, ..., 1, ..., 1, ..., 0)
    for (unsigned i = assignment_begin; i < table_dim - 1; ++i) {
      coords[i] = 1;
      for (unsigned j = i+1; j < table_dim; ++j) {
        coords[j] = 1;
        // LOG(logDEBUG4) << "TrainableModelBuilder::add_outgoing_factor(): divisions: "
        // << "entry to be erased: " << BinToDec(coords);
        size_t check = entries.erase(BinToDec(coords));
        // LOG(logDEBUG4) << "TrainableModelBuilder::add_outgoing_factor(): divisions: "
        // << "successfully erased? " << check;
        assert(check == 1);
        feature_type probability = 0.;
        if (has_classifiers()) {
          // FIXME: NEEDS TO BE IMPLENENTED
          throw std::runtime_error("TrainableModelBuilder does not support classifier priors yet!");
        }
        OpengmModel::ValueType value = division()(trax, arcs[i-assignment_begin].second, arcs[j-assignment_begin].second, probability);
        OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape.begin(), shape.end(), coords.begin(), value)
            .add_as_feature_to( *(m.opengm_model), m.weight_map[Model::div_weight].front() );
        coords[j] = 0;
      }
      coords[i] = 0;
    }
  }

  // forbidden configurations
  {
    for (std::set<size_t>::iterator it = entries.begin(); it != entries.end(); ++it) {
      coords = DecToBin(*it);
      if (coords.size() < table_dim) {
        coords.insert(coords.begin(), table_dim-coords.size(), 0);
      }
      assert (coords.size() == table_dim);
      OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape.begin(), shape.end(), coords.begin(), forbidden_cost())
          .add_to( *(m.opengm_model) );
    }
  }

  // LOG(logDEBUG1) << "TrainableModelBuilder::add_outgoing_factor(): leaving";
}


void TrainableModelBuilder::add_incoming_factor(const MultiHypothesesGraph& hypotheses,
                                                Model& m,
                                                const MultiHypothesesGraph::Node& node,
                                                const Traxel& trax,
                                                const std::vector<Traxel>& neighbors) const {
  // LOG(logDEBUG2) << "TrainableModelBuilder::add_incoming_factor(): entered for " << trax;
  // LOG(logDEBUG1) << "TrainableModelBuilder::add_incoming_factor(): entered";
  const std::vector<size_t> vi = vars_for_incoming_factor(hypotheses, m, node, trax);
  if (vi.size() == 0) {
    // nothing to be done here
    // no det vars and no incoming arcs
    return;
  }

  // construct factor
  LOG(logDEBUG2) << "TrainableModelBuilder::add_incoming_factor: constructing factor for "
  << trax << ": " << std::pow(2., static_cast<int>(vi.size())) << " entries (2^" << vi.size() << ")";
  const size_t table_dim = vi.size();
  const std::vector<size_t> shape(table_dim, 2);
  std::vector<size_t> coords;

  std::set<size_t> entries;
  for (size_t i = 0; i < static_cast<size_t>(std::pow(2., static_cast<int>(table_dim))); ++i) {
    entries.insert( entries.end(), i );
  }

  // opportunity ??

  // appearance
  coords = std::vector<size_t>(table_dim, 0);
  if (has_detection_vars()) {
    // *(coords.rbegin()) = 1; reverse or not reverse in vars_for_incoming_factor?
    coords[0] = 1; // (1, 0, ..., 0)
  }
  size_t check = entries.erase(BinToDec(coords));
  assert(check == 1);
  double weight = appearance()(trax);
  if (trax.Timestep == hypotheses.earliest_timestep()) {
    weight = 0.;
  }
  LOG(logDEBUG2) << "TrainableModelBuilder::add_incoming_factor: appearance cost for "
                   << trax << ": " << weight;
  OpengmWeightedFeature<OpengmModel::ValueType>( vi, shape.begin(), shape.end(), coords.begin(), weight )
      .add_as_feature_to( *(m.opengm_model), m.weight_map[Model::app_weight].front() );

  // moves
  // (1, 0, 0, ..., 1, ..., 0)
  LOG(logDEBUG2) << "TrainableModelBuilder::add_incoming_factor: moves for " << trax;
  coords = std::vector<size_t>(table_dim, 0);
  size_t assignment_begin = 0;
  if (has_detection_vars()) {
    coords[0] = 1;
    assignment_begin = 1;
  }
  for (size_t i = assignment_begin; i < table_dim; ++i) {
    coords[i] = 1;
    size_t check = entries.erase(BinToDec(coords));
    assert(check == 1);
    coords[i] = 0;
  }

  // forbidden configurations
  LOG(logDEBUG2) << "TrainableModelBuilder::add_incoming_factor: forbidden configurations for " << trax;
  for (std::set<size_t>::iterator it = entries.begin(); it != entries.end(); ++it) {
    coords = DecToBin(*it);
    // zero padding up to table dim
    if(coords.size() < table_dim) {
      coords.insert(coords.begin(), table_dim-coords.size(), 0);
    }
    assert(coords.size() == table_dim);
    OpengmWeightedFeature<OpengmModel::ValueType>(vi, shape.begin(), shape.end(), coords.begin(), forbidden_cost())
        .add_to( *(m.opengm_model) );
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

  boost::shared_ptr<Model> model(new Model);


  if (has_detection_vars()) {
    add_detection_vars( hypotheses, *model );
  }
  if (has_maximum_arcs()) {
    if (!has_classifiers()) {
      // throw std::runtime_error("maximum arcs required without classifiers present!");
    }
    const MultiHypothesesGraph::MoveFeatureMap& moves = hypotheses.get(node_move_features());
    add_assignment_vars( hypotheses, *model, moves );
  } else {
    add_assignment_vars( hypotheses, *model );
  }

  MultiHypothesesGraph::node_timestep_map& timesteps = hypotheses.get(node_timestep());
  if ( has_detection_vars() ) {
    for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
      if (timestep_range_specified() &&
          (timesteps[n] < first_timestep() || timesteps[n] > last_timestep())) {
        continue;
      }
      if (has_maximal_conflict_cliques()) {
        add_conflict_factors( hypotheses, *model, n );
      } else {
        add_detection_factors( hypotheses, *model, n );
      }
    }
  }

  add_count_factors(hypotheses, *model);

  for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
    if (timestep_range_specified() &&
        (timesteps[n] < first_timestep() || timesteps[n] > last_timestep())) {
      continue;
    }
    add_outgoing_factors( hypotheses, *model, n );
    add_incoming_factors( hypotheses, *model, n );
  }

  return model;
}


void CVPR2014ModelBuilder::add_detection_factors( const MultiHypothesesGraph& hypotheses,
                                                   Model& m,
                                                   const MultiHypothesesGraph::Node& n) const {
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  const std::vector<Traxel>& traxels = regions[n];
  for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
    add_detection_factor(m, *t);
  }
}

void CVPR2014ModelBuilder::add_count_factors( const MultiHypothesesGraph& hypotheses, Model& m) {
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  MultiHypothesesGraph::node_timestep_map& timesteps = hypotheses.get(node_timestep());
  for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
    if (timestep_range_specified() &&
        (timesteps[n] < first_timestep() || timesteps[n] > last_timestep())) {
      continue;
    }
    const std::vector<Traxel>& traxels = regions[n];
    if (has_maximal_conflict_cliques()) {
      MultiHypothesesGraph::ConflictSetMap& conflicts = hypotheses.get(node_conflict_sets());
      add_count_factor(m, traxels, (conflicts[n]).size());
    } else {
      add_count_factor(m, traxels, traxels.size());
    }
  }
}


void CVPR2014ModelBuilder::add_conflict_factors( const MultiHypothesesGraph& hypotheses, Model& m, const MultiHypothesesGraph::Node& n ) const {
  LOG(logDEBUG3) << "CVPR2014ModelBuilder::add_conflict_factor() -- add factors for conflict sets";
  const ConflictSets& conflicts = hypotheses.get(node_conflict_sets())[n];
  int timestep = hypotheses.get(node_timestep())[n];
  for (ConflictSets::const_iterator conflict = conflicts.begin(); conflict != conflicts.end(); ++conflict) {
    add_conflict_factor( hypotheses, m, n, *conflict, timestep);
  }
}

void CVPR2014ModelBuilder::add_conflict_factor( const MultiHypothesesGraph& hypotheses,
                                                Model& m,
                                                const MultiHypothesesGraph::Node& n,
                                                const ConflictSet& conflict,
                                                int timestep ) const {
  // MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  const std::vector<Traxel>& traxels_in_component = hypotheses.get(node_regions_in_component())[n];
  std::vector<size_t> vi;
  std::vector<Traxel> traxels;
  for (ConflictSet::const_iterator id = conflict.begin(); id != conflict.end(); ++id) {
    traxels.push_back(*std::find(traxels_in_component.begin(), traxels_in_component.end(), Traxel(*id, timestep)));
    vi.push_back(m.var_of_trax(*traxels.rbegin()));
  }
  size_t table_dim = vi.size();
  std::vector<size_t> coords(table_dim, 0);
  OpengmExplicitFactor<double> table(vi);
  double deactivated_sum = 0.0;
  for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
    deactivated_sum += detection()(*t, 0);
  }
  table.set_value( coords, deactivated_sum );

  for (size_t i = 0; i < table_dim; ++i) {
    coords[i] = 1;
    table.set_value( coords, detection()(traxels[i], 1 ));
    coords[i] = 0;
  }
}

void CVPR2014ModelBuilder::add_outgoing_factors( const MultiHypothesesGraph& hypotheses,
                                                  Model& m,
                                                  const MultiHypothesesGraph::Node& n ) const {
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  const std::vector<Traxel>& traxels = regions[n];
  for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
    std::vector<Traxel> neighbors;
    for (MultiHypothesesGraph::OutArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
      const std::vector<Traxel>& neighbors_at = regions[hypotheses.target(a)];
      neighbors.insert(neighbors.end(), neighbors_at.begin(), neighbors_at.end());      
    }
    add_outgoing_factor(hypotheses, m, n, *t, neighbors);
  }

}


void CVPR2014ModelBuilder::add_incoming_factors( const MultiHypothesesGraph& hypotheses,
                                                   Model& m,
                                                   const MultiHypothesesGraph::Node& n ) {
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  const std::vector<Traxel>& traxels = regions[n];
  std::vector<Traxel> neighbors;
  for (MultiHypothesesGraph::InArcIt a(hypotheses, n); a != lemon::INVALID; ++a) {
    const std::vector<Traxel>& neighbors_at = regions[hypotheses.source(a)];
    neighbors.insert(neighbors.end(), neighbors_at.begin(), neighbors_at.end());
  }
  for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
    add_incoming_factor(hypotheses, m, n, *t, neighbors);
  }

}


void CVPR2014ModelBuilder::add_detection_factor( Model& m,
                                                 const Traxel& trax) const {
  LOG(logDEBUG3) << "CVPR2014ModelBuilder::add_detection_factor() -- entered";
  size_t vi[] = {m.var_of_trax(trax)};
  std::vector<size_t> coords(1, 0);
  OpengmExplicitFactor<double> table(vi, vi+1);

  coords[0] = 0;
  table.set_value( coords, non_detection()(trax, 0) );
  // table.set_value( coords, 0. );


  coords[0] = 1;
  table.set_value( coords, detection()(trax, 1) );

  // table = OpengmExplicitFactor<double>(vi, vi+1, 0);
  LOG(logDEBUG3) << "CVPR2014ModelBuilder::add_detection_factor: for "
                 << trax << ": detection=" << table.get_value(std::vector<size_t>(1,1))
                 << ", non_detection=" << table.get_value(std::vector<size_t>(1,0));
  table.add_to( *(m.opengm_model) );

}


void CVPR2014ModelBuilder::add_count_factor( Model& m,
                                     const std::vector<Traxel>& traxels,
                                     size_t maximum_active_regions) {
	LOG(logDEBUG) << "ModelBuilder::add_count_factor: entered";

	if (has_hierarchical_counting_factor()) {
		add_hierarchical_count_factor(m,traxels,maximum_active_regions);
	} else {
		add_explicit_count_factor(m,traxels,maximum_active_regions);
	}
}

void CVPR2014ModelBuilder::add_hierarchical_count_factor( Model& m,
                                                  const std::vector<Traxel>& traxels,
                                                  size_t maximum_active_regions) {
	LOG(logDEBUG3) << "CVPR2014ModelBuilder::add_hierarchical_count_factor: entered";
	assert(has_hierarchical_counting_factor());
	assert(traxels.size() > 0);
  	assert(traxels[0].features.find("count_prediction") != traxels[0].features.end());

  	feature_array probabilities = traxels[0].features.find("count_prediction")->second;
  	assert(probabilities.size() > 0);

  	// probabilities must have one state more than maximum_active_regions to
  	// include the state for all regions deactivated
  	fill_probabilities(probabilities, maximum_active_regions);

  	// add counter variable
  	size_t num_states = std::min(traxels.size() + 1, maximum_active_regions + 1);
  	size_t var_id;
  	if (num_states > 2) {
		m.opengm_model->addVariable(num_states);
		var_id = m.opengm_model->numberOfVariables() - 1;
		assert(m.opengm_model->numberOfLabels(var_id) == num_states);

		// add helper hard constraints to be added later
		std::vector<std::pair<std::pair<size_t, size_t>, int> > constraint;
		std::vector<std::pair<size_t,size_t> > vars;
		for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
			std::pair<size_t, size_t> var_state = std::pair<size_t, size_t>(m.var_of_trax(*t), 1);
			std::pair<std::pair<size_t, size_t>, int> var_coeff(var_state, 1);
			constraint.push_back(var_coeff);
		}
		for (size_t state = 0; state < num_states; ++state) {
			std::pair<size_t, size_t> var_state = std::pair<size_t, size_t>(var_id, state);
			std::pair<std::pair<size_t, size_t>, int> var_coeff(var_state, -state);
			constraint.push_back(var_coeff);
		}
		var_state_coeff_constraints_.push_back(constraint);
  	} else {
  		var_id = m.var_of_trax(traxels[0]);
  	}


  	// add the count prior
	size_t var[] = {var_id};
	assert(m.opengm_model->numberOfLabels(var[0]) == probabilities.size());

	std::vector<size_t> coords(1, 0);
	OpengmExplicitFactor<double> table(var, var+1, 0, probabilities.size());

	for(size_t state = 0; state < probabilities.size(); ++state) {
		coords[0] = state;
                if (state == 0) {
                  table.set_value( coords, count()(probabilities[state]) + opportunity_cost() );
                } else {
                  table.set_value( coords, count()(probabilities[state]) );
                }
	}

	LOG(logDEBUG3) << "CVPR2014ModelBuilder::add_hierarchical_count_factor: "
                      << "helper count factor added";
	table.add_to( *(m.opengm_model) );
}


//void ModelBuilder::add_hierarchical_count_factor( Model& m,
//                                               const std::vector<Traxel>& traxels,
//                                               size_t maximum_active_regions) {
//	LOG(logDEBUG) << "ModelBuilder::add_hierarchical_count_factor: entered";
//	assert(has_hierarchical_counting_factor());
//	assert(traxels.size() > 0);
//  	assert(traxels[0].features.find("count_prediction") != traxels[0].features.end());
//
//  	feature_array probabilities = traxels[0].features.find("count_prediction")->second;
//  	assert(probabilities.size() > 0);
//
//  	// probabilities must have one state more than maximum_active_regions to
//  	// include the state for all regions deactivated
//  	fill_probabilities(probabilities, maximum_active_regions);
//
//  	std::vector<std::vector<std::pair<size_t, size_t> > > lower_level;
//	for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
//		std::vector<std::pair<size_t,size_t> > vars;
//		vars.push_back(std::pair<size_t, size_t>(m.var_of_trax(*t), (size_t) 0));
//		vars.push_back(std::pair<size_t, size_t>(m.var_of_trax(*t), (size_t) 1));
//		lower_level.push_back(vars);
//	}
//
//	std::vector<std::vector<std::pair<size_t, size_t> > > higher_level;
//	do {
//		add_count_helper(m, lower_level, maximum_active_regions, higher_level);
//	} while(higher_level.size() > 1);
//
//	assert(higher_level.size() == 1);
//
//	// add the count prior
//	size_t var_id[] = {m.opengm_model->numberOfVariables() - 1};
//	assert(m.opengm_model->numberOfLabels(var_id[0]) == probabilities.size());
//
//	std::vector<size_t> coords(1, 0);
//	OpengmExplicitFactor<double> table(var_id, var_id+1);
//
//	for(size_t state = 0; state < probabilities.size(); ++state) {
//		coords[0] = state;
//		table.set_value( coords, probabilities[state] );
//	}
//
//	LOG(logDEBUG2) << "ModelBuilder::add_hierarchical_count_factor: "
//	                 << "helper count factor added";
//	table.add_to( *(m.opengm_model) );
//}


void CVPR2014ModelBuilder::add_count_helper( Model& m,
                                     std::vector<std::vector<std::pair<size_t, size_t> > >& lower_level, const size_t maximum_active_regions,
                                     std::vector<std::vector<std::pair<size_t, size_t> > >& higher_level) {

	LOG(logDEBUG1) << "add_count_helper: entered";

	// add helper variables
	for (size_t i = 0; i < lower_level.size() - 1;  i += 2) { // note the integer division!
		std::vector<std::pair<size_t, size_t> > vars;
		size_t num_states = std::min(lower_level[i].size() + lower_level[i+1].size() - 1, maximum_active_regions + 1);
		m.opengm_model->addVariable(num_states);
		size_t var_id = m.opengm_model->numberOfVariables() - 1;
		assert(m.opengm_model->numberOfLabels(var_id) == num_states);
		for (size_t j = 0; j < num_states; ++j) {
			vars.push_back(std::pair<size_t, size_t>(var_id, j));
		}
		higher_level.push_back(vars);
	}
	if (lower_level.size() % 2 == 1) {
		higher_level.push_back(lower_level.back());
	}

	// add helper hard constraints to be added later
	std::vector<std::pair<std::pair<size_t, size_t>, int> > constraint;
	for (size_t i = 0; i < lower_level.size() - 1;  i += 2) { // note the integer division!
		for (size_t nu = 0; nu < lower_level[i].size(); ++nu) {
			std::pair<std::pair<size_t, size_t>, int> var_coeff(lower_level[i][nu], nu);
			constraint.push_back(var_coeff);
		}
		for (size_t nu = 0; nu < lower_level[i+1].size(); ++nu) {
			std::pair<std::pair<size_t, size_t>, int> var_coeff(lower_level[i+1][nu], nu);
			constraint.push_back(var_coeff);
		}
		for (size_t nu = 0; nu < higher_level[i/2].size(); ++nu) {
			std::pair<std::pair<size_t, size_t>, int> var_coeff(higher_level[i/2][nu], -nu);
			constraint.push_back(var_coeff);

		}
      var_state_coeff_constraints_.push_back(constraint);
	}
	return;
}

void ModelBuilder::add_count_hard_constraints(const Model& m,
                                                      const MultiHypothesesGraph& hypotheses, OpengmLPCplex& cplex) const {
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



void CVPR2014ModelBuilder::add_explicit_count_factor( Model& m,
                                              const std::vector<Traxel>& traxels,
                                              size_t maximum_active_regions ) const {
  LOG(logDEBUG) << "CVPR2014ModelBuilder::add_explicit_count_factor: entered";
  assert(not has_hierarchical_counting_factor());
  assert(traxels.size() > 0);
  assert(traxels[0].features.find("count_prediction") != traxels[0].features.end());
  feature_array probabilities = traxels[0].features.find("count_prediction")->second;
  size_t table_dim = traxels.size();
  
  assert(probabilities.size() > 0);

  // probabilities must have one state more than maximum_active_regions to
  // include the state for all regions deactivated
  fill_probabilities(probabilities, maximum_active_regions);

  vector<size_t> vi;
    for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
      vi.push_back(m.var_of_trax(*t));
    }

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
        if (active_count == 0) {
          table.set_value(coords, count()(probabilities[active_count]) + opportunity_cost());
        } else {
          table.set_value(coords, count()(probabilities[active_count]));
        }
      }
      coords = std::vector<size_t>(table_dim, 0);
    }
    table.add_to( *m.opengm_model );
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



void CVPR2014ModelBuilder::add_outgoing_factor( const MultiHypothesesGraph& hypotheses,
                                                Model& m,
                                                const MultiHypothesesGraph::Node& node,
                                                const Traxel& trax,
                                                const std::vector<Traxel>& neighbors ) const {
  
  const vector<size_t> vi = vars_for_outgoing_factor(hypotheses, m, node, trax);
  LOG(logDEBUG2) << "CVPR2014ModelBuilder::add_outgoing_factor(): entered for " << trax
                 << " at " << trax.features.find("com")->second[0] << "," << trax.features.find("com")->second[1] << ","
                 << trax.features.find("com")->second[2] << " lvl: " << trax.Level << " - factor order: " << vi.size();
  if (vi.size() == 0) {
    // nothing to do here
    // happens in case of no det vars and no outgoing arcs
    return;
  }

  // collect TraxelArcs for use in feature functions
  std::vector<Model::TraxelArc> arcs;
  for (std::vector<size_t>::const_iterator v = vi.begin()+1;
       v != vi.end();
       ++v) {
    arcs.push_back(m.arc_of_var(*v));
  }

  size_t table_dim = vi.size();
  assert(table_dim <= neighbors.size() + 1);

  // construct factor
  // only one detection var no outgoing arcs
  /* if (table_dim == 1) {
    std::vector<size_t> coords(table_dim, 0);
    OpengmExplicitFactor<double> table( vi );

    // oppportunity?
    table.set_value( coords, 0 );
    
    // disappearance
    if (trax.Timestep < hypotheses.latest_timestep()) {
      coords[0] = 1;
      table.set_value( coords, disappearance()(trax) );
      coords[0] = 0;
    }

    table.add_to( *m.opengm_model );

  } else if( table_dim == 2) {
    // no division possible
    std::vector<size_t> coords(table_dim, 0);
    OpengmExplicitFactor<double> table( vi, forbidden_cost() );

    // opportunity?
    table.set_value( coords, 0 );

    // disappearance configuration
    if (trax.Timestep < hypotheses.latest_timestep()) {
      coords[0] = 1;
      table.set_value( coords, disappearance()(trax) );
      coords[0] = 0;
    }

    // move configuration
    feature_type probability = 0.;
    if (has_classifiers()) {
      // messy! needs better implementation
      probability = hypotheses.get(node_move_features())[node]
          .find(trax)->second
          .find(neighbors[0])->second[0];
    }
    coords[0] = 1; coords[1] = 1;
    table.set_value( coords, move()(trax, neighbors[0], probability) );
    coords[0] = 0; coords[1] = 0;

    table.add_to( *m.opengm_model ); */

  // } else {
    
    std::vector<size_t> coords(table_dim, 0);
    OpengmExplicitFactor<double> table( vi, forbidden_cost() );

    // opportunity?
    table.set_value( coords, 0 );

    // disappearance configuration
    if (trax.Timestep < hypotheses.latest_timestep()) {
      coords[0] = 1;
      table.set_value( coords, disappearance()(trax) );
      LOG(logDEBUG3) << "CVPR2014ModelBuilder::add_outgoing_factor: at least two outgoing arcs: "
                     << "forbidden=" << forbidden_cost() << ", disappearance=" << table.get_value(coords);
      coords[0] = 0;
    }



    // move configuration

    coords[0] = 1;
    for (size_t i = 1; i < table_dim; ++i) {
      coords[i] = 1;
      feature_type probability = 0.;
      if (has_classifiers()) {
        probability = hypotheses.get(node_move_features())[node]
            .find(trax)->second
            .find(arcs[i-1].second)->second[1];
        LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_outgoing_factor: move using classifier, prob: " << probability
                       << ", energy: " << move()(trax, arcs[i-1].second, probability);
      }
     
      table.set_value( coords, move()(trax, arcs[i-1].second, probability) );
      LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_outgoing_factor: move="
                     << table.get_value( coords );
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
                .find(trax)->second
                .find(std::make_pair(arcs[i-1].second, arcs[j-1].second))->second[0];
          }     
          table.set_value(coords, division()(trax,
                                             arcs[i-1].second,
                                             arcs[j-1].second,
                                             probability
                                             ));
          LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_outgoing_factor: division="
                         << table.get_value( coords );
          coords[j] = 0;
        }
        coords[i] = 0;
      }
      coords[0] = 0;
    }
    // table = OpengmExplicitFactor<double>( vi, 999999 );
    table.add_to( *m.opengm_model );

    // }
  // LOG(logDEBUG2) << "CVPR2014ModelBuilder::add_outgoing_factor(): leaving";
}


void CVPR2014ModelBuilder::add_incoming_factor(const MultiHypothesesGraph& hypotheses,
                                                Model& m,
                                                const MultiHypothesesGraph::Node& node,
                                                const Traxel& trax,
                                                const std::vector<Traxel>& neighbors) {
  LOG(logDEBUG2) << "CVPR2014ModelBuilder::add_incoming_factor(): entered for " << trax;
  const std::vector<size_t> vi = vars_for_incoming_factor(hypotheses, m, node, trax);
  if (vi.size() == 0) {
    // nothing to be done here
    // no det vars and no incoming arcs
    return;
  }

  // construct factor
  LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_incoming_factor: " << vi.size()
                 << " - " <<  neighbors.size();

  assert(vi.size() <= neighbors.size() + has_detection_vars() ? 1 : 0);

  if (has_counting_incoming_factor()) {
	  // create one variable representing the sum of all incoming variables, then
	  // the pairwise factor between this counter and the detection variable is much cheaper
	  LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_incoming_factor: constructing count-helper";

	  vector<size_t> vi_pairwise;
	  vi_pairwise.push_back(vi[0]);
	  if (vi.size() > 2) {
		  m.opengm_model->addVariable(2); // add a binary counting variable (count can only be 0 or 1 for incomings)
		  size_t vi_from = m.opengm_model->numberOfVariables() - 1;

		  // add hard constraints assuring the counting (to be added later)
	      std::vector<std::pair<std::pair<size_t, size_t>, int> > constraint;
		  std::vector<std::pair<size_t,size_t> > vars;
		  for (size_t i = 1; i < vi.size(); ++i) {
				std::pair<size_t, size_t> var_state = std::pair<size_t, size_t>(vi[i], 1);
				std::pair<std::pair<size_t, size_t>, int> var_coeff(var_state, 1);
				constraint.push_back(var_coeff);
		  }
		  std::pair<size_t, size_t> var_state = std::pair<size_t, size_t>(vi_from, 1);
		  std::pair<std::pair<size_t, size_t>, int> var_coeff(var_state, -1);
		  constraint.push_back(var_coeff);
		  var_state_coeff_constraints_.push_back(constraint);
		  vi_pairwise.push_back(vi_from);
	  } else if (vi.size() == 2) {
		  size_t vi_from = vi[1]; // the only incoming variable
		  vi_pairwise.push_back(vi_from);
	  } else if (vi.size() == 1) {
		  // do nothing
	  } else { // vi.size() == 0
		  assert(false); // cannot happen, it is already returned earlier in this function
	  }

	  // create a table for a pairwise factor
	  assert(vi_pairwise.size() <= 2);
	  assert(vi_pairwise.size() > 0);
	  const size_t table_dim = vi_pairwise.size();
	  std::vector<size_t> coords(table_dim, 0);
	  OpengmExplicitFactor<double> table( vi_pairwise, 0 );

	  // appearance
	  if (trax.Timestep > hypotheses.earliest_timestep()) {
		  coords[0] = 1;
		  table.set_value( coords, appearance()(trax) );
		  assert(table.get_value( coords ) == appearance()(trax) );
		  LOG(logDEBUG2) << "CVPR2014ModelBuilder::add_incoming_factor: appearance="
							   << table.get_value( coords );

		  table.add_to( *m.opengm_model );
		  LOG(logDEBUG2) << "CVPR2014ModelBuilder::add_incoming_factor: done";
	  }
  } else { // create a table with 2^n rows
	  LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_incoming_factor: constructing factor for "
	  << trax << ": " << std::pow(2., static_cast<int>(vi.size())) << " entries (2^" << vi.size() << ")";
	  const size_t table_dim = vi.size();
	  std::vector<size_t> coords(table_dim, 0);
	  OpengmExplicitFactor<double> table( vi, forbidden_cost() );

	  // opportunity?
	  table.set_value( coords, 0 );

	  // appearance
	  if (trax.Timestep > hypotheses.earliest_timestep()) {
		coords[0] = 1;
		table.set_value( coords, appearance()(trax));
		assert(table.get_value( coords ) == appearance()(trax) );
		LOG(logDEBUG2) << "CVPR2014ModelBuilder::add_incoming_factor: appearance="
					   << table.get_value( coords );
		coords[0] = 0;
	  }

	  // move
	  // move cost already included in outgoing factor
	  coords[0] = 1;
	  for (size_t i = 1; i < table_dim; ++i) {
		coords[i] = 1;
		table.set_value( coords, 0 );
		LOG(logDEBUG4) << "CVPR2014ModelBuilder::add_incoming_factor: move="
					   << table.get_value( coords );
		coords[i] = 0;
	  }
	  coords[0] = 0;

	  // table = OpengmExplicitFactor<double>( vi, 999999 );
	  table.add_to( *m.opengm_model );
	  LOG(logDEBUG2) << "CVPR2014ModelBuilder::add_incoming_factor: done";
  }
}


    

} /* namespace multihypotheses */
} /* namespace pgm */
} /* namespace pgmlink */ 
