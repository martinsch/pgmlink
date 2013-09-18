/**
   @file
   @ingroup pgm
   @brief Graphical Model for Multi Hypotheses Segmentation
*/

#ifndef PGMLINK_PGM_MULTI_HYPOTHESES_H
#define PGMLINK_PGM_MULTI_HYPOTHESES_H


// stl
#include <algorithm>
#include <iterator>
#include <map>
#include <vector>
#include <utility>

// boost
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bimap.hpp>

// opengm
#include <opengm/inference/inference.hxx>
#include <opengm/inference/lpcplex.hxx>

// pgmlink
#include "pgmlink/pgm.h"
#include "pgmlink/multi_hypotheses_graph.h"
#include "pgmlink/feature.h"
#include "pgmlink/traxels.h"


namespace pgmlink {
namespace pgm {
namespace multihypotheses {
typedef opengm::LPCplex<OpengmModel, opengm::Minimizer> OpengmLPCplex;
using boost::function;
using std::map;
using std::vector;

class ModelBuilder;
/**
   @brief Opengm formulation for the multi hypotheses tracking
	 
   Represents an opengm model to solve the matching problem in a
   Multi HypothesesGraph. Use a multihypotheses::ModelBuilder to construct
   the model. During construction of the multihypotheses::Model a random
   variable is added to the graphical model for every... 
	 
   @see multihypotheses::ModelBuilder
   @see MultiHypothesesGraph
*/
class Model {
 public:

  
  typedef std::pair<Traxel, Traxel> TraxelArc;
  typedef OpengmModel::IndexType var_t;
  typedef boost::bimap<Traxel, var_t>::left_map trax_var_map;
  typedef boost::bimap<Traxel, var_t>::right_map var_trax_map;
  typedef boost::bimap<TraxelArc, var_t>::left_map arc_var_map;
  typedef boost::bimap<TraxelArc, var_t>::right_map var_arc_map;

  Model();
  Model( shared_ptr<OpengmModel>,
         const trax_var_map&,
         const arc_var_map&
         );

  shared_ptr<OpengmModel> opengm_model; // opengm model usually constructed by multihypotheses::ModelBuilder

  const trax_var_map& var_of_trax() const; // maps traxels to random variables representing detections
  const var_trax_map& trax_of_var() const;
  const arc_var_map& var_of_arc() const; // maps pairs of traxels to random variables representing links
  const var_arc_map& arc_of_var() const;

  var_t var_of_trax(const Traxel&) const;
  var_t var_of_arc(const TraxelArc&) const;
  Traxel trax_of_var(var_t) const;
  TraxelArc arc_of_var(var_t) const;

  enum VarCategory {trax_var, arc_var};
  VarCategory var_category(var_t) const;

  enum WeightType {det_weight, mov_weight, div_weight, app_weight, dis_weight};
  map<WeightType, vector<OpengmModel::IndexType> > weight_map; // associate events with weight ids

 private:
  friend class ModelBuilder;

  void init();

  boost::bimap<Traxel, var_t> trax_var_;
  boost::bimap<TraxelArc, var_t> arc_var_;
  
};

} /* namespace multihypotheses */
} /* namespace pgm */
} /* namespace pgmlink */



#endif /* PGMLINK_PGM_MULTI_HYPOTHESES_H */
