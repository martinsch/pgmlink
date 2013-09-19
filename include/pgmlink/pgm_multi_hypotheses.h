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



class ModelBuilder {
 public:
  enum MAX_OUTGOING_ARCS {TRANSITION=1, DIVISION=2};
  ModelBuilder(boost::function<double (const Traxel&)> appearance = ConstantFeature(1000),
               boost::function<double (const Traxel&)> disappearance = ConstantFeature(1000),
               boost::function<double (const Traxel&, const Traxel&)> move = SquaredDistance(),
               double forbidden_cost = 1000000,
               unsigned max_division_level=0,
               unsigned max_count=2)
      : with_detection_vars_(false),
        with_divisions_(false),
        appearance_(appearance),
        disappearance_(disappearance),
    move_(move),
    forbidden_cost_(forbidden_cost),
    cplex_timeout_(1e+75),
    max_division_level_(max_division_level),
    max_count_(max_count) {}

  virtual multihypotheses::ModelBuilder* clone() const = 0;
  virtual ~ModelBuilder() {}

  // mandatory parameters
  function<double (const Traxel&)> appearance() const { return appearance_; }
  ModelBuilder& appearance( function<double (const Traxel&)> );

  function<double (const Traxel&)> disappearance() const { return disappearance_; }
  ModelBuilder& disappearance( function<double (const Traxel&)> );

  function<double (const Traxel&, const Traxel&)> move() const { return move_; }
  ModelBuilder& move(function<double (const Traxel&, const Traxel&)> );

  double forbidden_cost() const { return forbidden_cost_; }
  ModelBuilder& forbidden_cost( double c ) { forbidden_cost_ = c; return *this; }

  // optional parameters
  // detection vars
  ModelBuilder& with_detection_vars( function<double (const Traxel&)> detection=ConstantFeature(10),
                                     function<double (const Traxel&)> non_detection=ConstantFeature(200));
  ModelBuilder& without_detection_vars();
  bool has_detection_vars() const { return with_detection_vars_; }
  function<double (const Traxel&)> detection() const { return detection_; }
  function<double (const Traxel&)> non_detection() const { return non_detection_; }

  // divisions
  ModelBuilder& with_divisions( function<double (const Traxel&, const Traxel&, const Traxel&)> div = KasterDivision(10) );
  ModelBuilder& without_divisions();
  bool has_divisions() const { return with_divisions_; }
  function<double (const Traxel&, const Traxel&, const Traxel&)> division() const { return division_; }

  // build
  virtual boost::shared_ptr<multihypotheses::Model> build( const MultiHypothesesGraph& ) const = 0;

  // refinement
  void add_hard_constraints( const Model&, const MultiHypothesesGraph&, OpengmLPCplex& );

  // cplex parameters
  void set_cplex_timeout( double seconds );


 protected:
  void add_detection_vars( const MultiHypothesesGraph&, Model& ) const;
  void add_assignment_vars( const MultiHypothesesGraph&, Model& ) const;

  vector<OpengmModel::IndexType> vars_for_outgoing_factor( const MultiHypothesesGraph&,
                                                           const Model&,
                                                           const MultiHypothesesGraph::Node&) const;
  vector<OpengmModel::IndexType> vars_for_incoming_factor( const MultiHypothesesGraph&,
                                                           const Model&,
                                                           const MultiHypothesesGraph::Node&) const;

 private:
  void couple_outgoing( const Model&, const std::vector<Traxel>&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_incoming( const Model&, const std::vector<Traxel>&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_detections_assignments( const Model&, const std::vector<Traxel>&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_outgoing_assignments( const Model&, const std::vector<Traxel>&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_incoming_assignments( const Model&, const std::vector<Traxel>&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_count( const Model&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_conflicts( const Model&, const std::vector<Traxel>&, OpengmLPCplex& );

  bool with_detection_vars_;
  bool with_divisions_;

  function<double (const Traxel&)> detection_;
  function<double (const Traxel&)> non_detection_;
  function<double (const Traxel&)> appearance_;
  function<double (const Traxel&)> disappearance_;
  function<double (const Traxel&, const Traxel&)> move_;
  function<double (const Traxel&, const Traxel&, const Traxel&)> division_;
  double forbidden_cost_;
  double cplex_timeout_;
  unsigned max_division_level_;
  unsigned max_count_;
  
};


} /* namespace multihypotheses */
} /* namespace pgm */
} /* namespace pgmlink */



#endif /* PGMLINK_PGM_MULTI_HYPOTHESES_H */
