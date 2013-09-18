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
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/pgm_multi_hypotheses.h"
#include "pgmlink/traxels.h"
#include "pgmlink/multi_hypotheses_graph.h"

//#include <ostream>

using namespace std;

namespace pgmlink {
namespace pgm {
namespace multihypotheses {
////
//// class chaingraph::Model
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

    

} /* namespace chaingraph */
} /* namespace pgm */
} /* namespace pgmlink */ 
