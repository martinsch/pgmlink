#include <algorithm>
#include <cassert>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>
#include <boost/scoped_ptr.hpp>
#include <opengm/inference/lpcplex.hxx>
#include <opengm/datastructures/marray/marray.hxx>

#include "pgmlink/pgm.h"
#include "pgmlink/log.h"
#include "pgmlink/pgm_multi_hypotheses.h"
#include "pgmlink/traxels.h"
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
    throw std::out_of_range("ChaingraphModel::var_of_arc(): key does not exist");
  }
}

Traxel Model::trax_of_var(var_t e) const {
  var_trax_map::const_iterator it = trax_of_var().find(e);
  if(it!=trax_of_var().end()) {
    return it->second;
  } else {
    throw std::out_of_range("ChaingraphModel::trax_of_var(): key does not exist");
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
ModelBuilder& ModelBuilder::with_detection_vars(  function<double (const Traxel&)> detection,
                                                  function<double (const Traxel&)> non_detection) {
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


ModelBuilder& ModelBuilder::with_divisions( function<double (const Traxel&, const Traxel&, const Traxel&)> division ) {
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

namespace {
inline size_t cplex_id(size_t opengm_id) {
  return 2*opengm_id+1;
}
}


void ModelBuilder::add_hard_constraints(const Model& m, const MultiHypothesesGraph& hypotheses, OpengmLPCplex& cplex) {
  LOG(logDEBUG) << "MultiHypotheses::add_hard_constraints: entered";
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());

  
  for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
    LOG(logDEBUG1) << "MultiHypotheses::add_hard_constraints: outgoing transitions";
    const std::vector<Traxel>& traxels = regions[n];
    for (MultiHypothesesGraph::OutArcIt a(hypotheses, n); a != lemon::INVALID; ++n) {
      const std::vector<Traxel>& traxels_dest = regions[hypotheses.target(a)];
      couple_outgoing(m, traxels, traxels_dest, cplex);
    }

    LOG(logDEBUG1) << "MultiHypotheses::add_hard_constraints: incoming transitions";
    for (MultiHypothesesGraph::InArcIt a(hypotheses, n); a != lemon::INVALID; ++n) {
      const std::vector<Traxel>& traxels_src = regions[hypotheses.source(a)];
      couple_incoming(m, traxels_src, traxels, cplex);
    }

    if (has_detection_vars()) {
      LOG(logDEBUG1) << "MultiHypotheses::add_hard_constraints: coupling conflicts";
      couple_count(m, traxels, cplex);
      
      LOG(logDEBUG1) << "MultiHypotheses::add_hard_constraints: coupling count";
      couple_conflicts(m, traxels, cplex);
    }
  }
}


void ModelBuilder::set_cplex_timeout( double seconds ) {
  cplex_timeout_ = seconds;
}


inline void ModelBuilder::add_detection_vars( const MultiHypothesesGraph& hypotheses, Model& m ) const {
  if(!has_detection_vars()) {
    throw std::runtime_error("multihypotheses::ModelBuilder::add_detection_vars(): called without has_detection_vars()");
  }
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  for (MultiHypothesesGraph::NodeIt n(hypotheses); n != lemon::INVALID; ++n) {
    const std::vector<Traxel>& traxels = regions[n];
    for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
      m.opengm_model->addVariable(2);
      m.trax_var_.left.insert(Model::trax_var_map::value_type(*t, m.opengm_model->numberOfVariables() - 1));
    }
  }
}


inline void ModelBuilder::add_assignment_vars( const MultiHypothesesGraph& hypotheses, Model& m ) const {
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  for (MultiHypothesesGraph::ArcIt a(hypotheses); a != lemon::INVALID; ++a) {
    const std::vector<Traxel>& source = regions[hypotheses.source(a)];
    const std::vector<Traxel>& dest = regions[hypotheses.target(a)];
    for (std::vector<Traxel>::const_iterator s = source.begin(); s != source.end(); ++s) {
      for (std::vector<Traxel>::const_iterator d = dest.begin(); d != dest.end(); ++d) {
        m.opengm_model->addVariable(2);
        m.arc_var_.left.insert(Model::arc_var_map::value_type(Model::TraxelArc(*s, *d), m.opengm_model->numberOfVariables() - 1));
      }
    }
  }
}


std::vector<OpengmModel::IndexType> ModelBuilder::vars_for_outgoing_factor( const MultiHypothesesGraph& hypotheses,
                                                                            const Model& m,
                                                                            const MultiHypothesesGraph::Node& node,
                                                                            const Traxel& trax) const {
  std::vector<OpengmModel::IndexType> vi; // opengm variable indices; may be empty if no det vars
  if (has_detection_vars()) {
    vi.push_back(m.var_of_trax(trax)); // first detection var, the others will be assignment vars
  }
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  for (MultiHypothesesGraph::OutArcIt a(hypotheses, node); a != lemon::INVALID; ++a) {
    const std::vector<Traxel>& traxels = regions[hypotheses.target(a)];
    for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
      vi.push_back(m.var_of_arc(Model::TraxelArc(trax, *t)));
    }
  }
  return vi;
}


std::vector<OpengmModel::IndexType> ModelBuilder::vars_for_incoming_factor( const MultiHypothesesGraph& hypotheses,
                                                                            const Model& m,
                                                                            const MultiHypothesesGraph::Node& node,
                                                                            const Traxel& trax) const {
  std::vector<OpengmModel::IndexType> vi; // opengm variable indices; may be empty if no det vars
  if (has_detection_vars()) {
    vi.push_back(m.var_of_trax(trax)); // first detection var, the others will be assignment vars
  }
  MultiHypothesesGraph::ContainedRegionsMap& regions = hypotheses.get(node_regions_in_component());
  for (MultiHypothesesGraph::InArcIt a(hypotheses, node); a != lemon::INVALID; ++a) {
    const std::vector<Traxel>& traxels = regions[hypotheses.source(a)];
    for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
      vi.push_back(m.var_of_arc(Model::TraxelArc(*t, trax)));
    }
  }
  std::reverse(vi.begin(), vi.end()); // det var should be the first index to be consistent with vars_for_outgoing_factor() WHY? -> ASK BERNHARD!!
  return vi;
}


void ModelBuilder::couple_outgoing(const Model& m, const std::vector<Traxel>& source, const std::vector<Traxel>& dest, OpengmLPCplex& cplex) {
  if (has_detection_vars()) {
    couple_detections_assignments(m, source, dest, cplex);
  }
  couple_outgoing_assignments(m, source, dest, cplex);
}


void ModelBuilder::couple_incoming(const Model& m, const std::vector<Traxel>& source, const std::vector<Traxel>& dest, OpengmLPCplex& cplex) {
  if (has_detection_vars()) {
    couple_detections_assignments(m, source, dest, cplex);
  }
  couple_incoming_assignments(m, source, dest, cplex);
}


void ModelBuilder::couple_detections_assignments(const Model& m, const std::vector<Traxel>& source, const std::vector<Traxel>& dest, OpengmLPCplex& cplex) {
  for (std::vector<Traxel>::const_iterator s = source.begin(); s != source.end(); ++s) {
    for (std::vector<Traxel>::const_iterator d = dest.begin(); d != dest.end(); ++s) {
      std::vector<size_t> cplex_idxs;
      cplex_idxs.push_back(cplex_id(m.var_of_trax(*s)));
      cplex_idxs.push_back(cplex_id(m.var_of_arc(Model::TraxelArc(*s, *d))));
      std::vector<int> coeffs;
      coeffs.push_back(1);
      coeffs.push_back(-1);
      // 0 <= 1*detection - 1*transition <= 1
      cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 1);
    }
  }
}


void ModelBuilder::couple_outgoing_assignments(const Model& m, const std::vector<Traxel>& source, const std::vector<Traxel>& dest, OpengmLPCplex& cplex) {
  for (std::vector<Traxel>::const_iterator s = source.begin(); s != source.end(); ++s) {
    MAX_OUTGOING_ARCS max_outgoing_arcs = DIVISION;
    FeatureMap::const_iterator feature = s->features.find("level");
    if (feature == s->features.end()) {
      throw std::runtime_error("MultiHypotheses: couple_outgoing_assignments - Feature level not found in traxel!");
    }
    if (feature->second[0] > max_division_level_) {
      max_outgoing_arcs = TRANSITION;
    }
    std::vector<size_t> cplex_idxs;
    for (std::vector<Traxel>::const_iterator d = dest.begin(); d != dest.end(); ++d) {
      cplex_idxs.push_back(cplex_id(m.var_of_arc(Model::TraxelArc(*s, *d))));
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
  for (std::vector<Traxel>::const_iterator d = dest.begin(); d != dest.end(); ++d) {
    std::vector<size_t> cplex_idxs;
    for (std::vector<Traxel>::const_iterator s = source.begin(); s != source.end(); ++s) {
      cplex_idxs.push_back(cplex_id(m.var_of_arc(Model::TraxelArc(*s, *d))));
    }
    if (cplex_idxs.size() > 0) {
      std::vector<int> coeffs(cplex_idxs.size(), 1);
      // 0 <= 1*transition + ... + 1*transition <= 1
      cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 1);
    }
  }
}


void ModelBuilder::couple_count( const Model& m, const std::vector<Traxel>& traxels, OpengmLPCplex& cplex) {
  std::vector<size_t> cplex_idxs;
  for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
    cplex_idxs.push_back(cplex_id(m.var_of_trax(*t)));
  }
  if (cplex_idxs.size() > 0) {
    std::vector<int> coeffs(cplex_idxs.size(), 1);
    // 0 <= 1*detection + ... + 1*detection <= max_count_
    cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, max_count_);
  }
}


void ModelBuilder::couple_conflicts( const Model& m, const std::vector<Traxel>& traxels, OpengmLPCplex& cplex) {
  for (std::vector<Traxel>::const_iterator t = traxels.begin(); t != traxels.end(); ++t) {
    std::vector<size_t> cplex_idxs;
    FeatureMap::const_iterator feature = t->features.find("conflicts");
    if (feature == t->features.end()) {
      throw std::runtime_error("MultiHypotheses: couple_conflicts - Feature conflicts not found in traxel!");
    }
    cplex_idxs.push_back(cplex_id(m.var_of_trax(*t)));
    for (feature_array::const_iterator conflict = feature->second.begin(); conflict != feature->second.end(); ++conflict) {
      cplex_idxs.push_back(cplex_id(m.var_of_trax(*std::find(traxels.begin(), traxels.end(), Traxel(*conflict, t->Timestep)))));
    }
    if (cplex_idxs.size() > 1) {
      std::vector<int> coeffs(cplex_idxs.size(), 1);
      // 0 <= 1*detection + ... + 1*detection <= 1
      cplex.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 1);
    }
  }
}


////
//// class TrinableChaingraphModelBuilder
////
/* boost::shared_ptr<ModelBuilder> TrainableModelBuilder::clone() const {
  return boost::shared_ptr<ModelBuilder(new TrainableModelBuilder(*this))>;
} */


boost::shared_ptr<Model> TrainableModelBuilder::build(const MultiHypothesesGraph& hypotheses) const {
  boost::shared_ptr<Model> model(new Model);
  assert(model->opengm_model->numberOfWeights() == 0);
  model->opengm_model->increaseNumberOfWeights(3);
  model->weight_map[Model::mov_weight].push_back(0);
  model->weight_map[Model::app_weight].push_back(1);
  model->weight_map[Model::dis_weight].push_back(2);

  return model;
}
    

} /* namespace multihypotheses */
} /* namespace pgm */
} /* namespace pgmlink */ 
