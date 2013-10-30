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
bool operator==(const std::pair<Traxel, feature_type>& p, const Traxel& t);
bool operator==(const Traxel& t, const std::pair<Traxel, feature_type>& p);
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
               boost::function<double (const Traxel&, const Traxel&, feature_type)> move = SquaredDistance(),
               boost::function<double (feature_type)> count = ConstantFeature(100),
               double forbidden_cost = 1000000,
               double opportunity_cost = 1000,
               unsigned max_division_level=0,
               unsigned max_count=2)
      : with_detection_vars_(false),
    with_divisions_(false),
    with_classifier_priors_(false),
    with_maximal_conflict_cliques_(false),
    with_hierarchical_counting_factor_(false),
    with_counting_incoming_factor_(false),
    with_maximum_arcs_(false),
    with_timestep_range_(false),
    appearance_(appearance),
    disappearance_(disappearance),
    move_(move),
    count_(count),
    forbidden_cost_(forbidden_cost),
    opportunity_cost_(opportunity_cost),
    cplex_timeout_(1e+75),
    max_division_level_(max_division_level),
    max_count_(max_count),
    maximum_outgoing_arcs_(0) {}

  virtual boost::shared_ptr<ModelBuilder> clone() const = 0;
  virtual ~ModelBuilder() {}

  // mandatory parameters
  function<double (const Traxel&)> appearance() const { return appearance_; }
  ModelBuilder& appearance( function<double (const Traxel&)> );

  function<double (const Traxel&)> disappearance() const { return disappearance_; }
  ModelBuilder& disappearance( function<double (const Traxel&)> );

  function<double (const Traxel&, const Traxel&, feature_type)> move() const { return move_; }
  ModelBuilder& move(function<double (const Traxel&, const Traxel&, feature_type)> );

  function<double (feature_type)> count() const { return count_; }
  ModelBuilder& count( function<double (feature_type)> );

  double forbidden_cost() const { return forbidden_cost_; }
  ModelBuilder& forbidden_cost( double c ) { forbidden_cost_ = c; return *this; }

  double opportunity_cost() const { return opportunity_cost_; }
  ModelBuilder& opportunity_cost( double c ) { opportunity_cost_ = c; return *this; }

  // optional parameters
  // detection vars
  ModelBuilder& with_detection_vars( function<double (const Traxel&, size_t)> detection=ConstantFeature(10),
                                     function<double (const Traxel&, size_t)> non_detection=ConstantFeature(200));
  ModelBuilder& without_detection_vars();
  bool has_detection_vars() const { return with_detection_vars_; }
  function<double (const Traxel&, size_t)> detection() const { return detection_; }
  function<double (const Traxel&, size_t)> non_detection() const { return non_detection_; }

  // divisions
  ModelBuilder& with_divisions( function<double (const Traxel&, const Traxel&, const Traxel&, feature_type)> div = KasterDivision(10) );
  ModelBuilder& without_divisions();
  bool has_divisions() const { return with_divisions_; }
  function<double (const Traxel&, const Traxel&, const Traxel&, feature_type)> division() const { return division_; }

  // classifier priors
  ModelBuilder& with_classifier_priors( function<double (const Traxel&, const Traxel&, feature_type)> move,
                                        function<double (feature_type)> count );
  ModelBuilder& without_classifier_priors();
  bool has_classifiers() const { return with_classifier_priors_; }

  // maximal conflict cliques
  ModelBuilder& with_maximal_conflict_cliques(bool);
  bool has_maximal_conflict_cliques() const { return with_maximal_conflict_cliques_; }

  // build
  virtual boost::shared_ptr<Model> build( const MultiHypothesesGraph& ) = 0;

  // refinement
  void add_hard_constraints( const Model&, const MultiHypothesesGraph&, OpengmLPCplex& );

  // cplex parameters
  void set_cplex_timeout( double seconds );

  // hierarchical counting factor
  ModelBuilder& with_hierarchical_counting_factor(bool);
  bool has_hierarchical_counting_factor() const { return with_hierarchical_counting_factor_; }
  void add_count_hard_constraints(const Model& m,
                                  const MultiHypothesesGraph& hypotheses, OpengmLPCplex& cplex) const;

  // counter variable for incoming assignments
  ModelBuilder& with_counting_incoming_factor(bool);
  bool has_counting_incoming_factor() const { return with_counting_incoming_factor_; }

  // maximum outgoing, incoming arcs
  ModelBuilder& with_maximum_arcs(unsigned);
  ModelBuilder& without_maximum_arcs();
  bool has_maximum_arcs() const { return with_maximum_arcs_; }

  // with timestep range specified
  ModelBuilder& with_timestep_range(int, int);
  ModelBuilder& without_timestep_range();
  bool timestep_range_specified() const { return with_timestep_range_; }
  int first_timestep() const { return first_timestep_; }
  int last_timestep() const { return last_timestep_; }

 protected:
  size_t cplex_id(OpengmLPCplex& cplex, const size_t opengm_id, const size_t state) const;
  size_t cplex_id(OpengmLPCplex& cplex, const size_t opengm_id) const;
  void add_detection_vars( const MultiHypothesesGraph&, Model& ) const;
  void add_assignment_vars( const MultiHypothesesGraph&, Model& ) const;
  void add_assignment_vars( const MultiHypothesesGraph&, Model&, const MultiHypothesesGraph::MoveFeatureMap& ) const;

  vector<OpengmModel::IndexType> vars_for_outgoing_factor( const MultiHypothesesGraph&,
                                                           const Model&,
                                                           const MultiHypothesesGraph::Node&) const;
  vector<OpengmModel::IndexType> vars_for_incoming_factor( const MultiHypothesesGraph&,
                                                           const Model&,
                                                           const MultiHypothesesGraph::Node& ) const;

  std::vector<OpengmModel::IndexType> vars_for_outgoing_factor( const MultiHypothesesGraph&,
                                                                const Model&,
                                                                const MultiHypothesesGraph::Node&,
                                                                const Traxel& ) const;
  std::vector<OpengmModel::IndexType> vars_for_incoming_factor( const MultiHypothesesGraph&,
                                                                const Model&,
                                                                const MultiHypothesesGraph::Node&,
                                                                const Traxel& ) const;
  
  
  std::vector<std::vector<std::pair<std::pair<size_t, size_t>, int> > > var_state_coeff_constraints_;

 private:
  void couple_outgoing( const Model&, const std::vector<Traxel>&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_incoming( const Model&, const std::vector<Traxel>&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_detections_assignments( const Model&, const std::vector<Traxel>&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_detections_assignments_incoming( const Model&, const std::vector<Traxel>&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_outgoing_assignments( const Model&, const std::vector<Traxel>&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_incoming_assignments( const Model&, const std::vector<Traxel>&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_count( const Model&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_conflicts( const Model&, const MultiHypothesesGraph&, const MultiHypothesesGraph::Node&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_conflicts_pairwise( const Model&, const std::vector<Traxel>&, OpengmLPCplex& );
  void couple_conflicts_maximal_cliques( const Model&, int, const std::vector<std::vector<unsigned> >&, OpengmLPCplex& );

  bool with_detection_vars_;
  bool with_divisions_;
  bool with_classifier_priors_;
  bool with_maximal_conflict_cliques_;
  bool with_hierarchical_counting_factor_;
  bool with_counting_incoming_factor_;
  bool with_maximum_arcs_;
  bool with_timestep_range_;

  function<double (const Traxel&, size_t)> detection_;
  function<double (const Traxel&, size_t)> non_detection_;
  function<double (const Traxel&)> appearance_;
  function<double (const Traxel&)> disappearance_;
  function<double (const Traxel&, const Traxel&, feature_type)> move_;
  function<double (const Traxel&, const Traxel&, const Traxel&, feature_type)> division_;
  function<double (feature_type)> count_;
  double forbidden_cost_;
  double opportunity_cost_;
  double cplex_timeout_;
  unsigned max_division_level_;
  unsigned max_count_;  
  int maximum_outgoing_arcs_;
  int first_timestep_;
  int last_timestep_;
};


class TrainableModelBuilder : public ModelBuilder {
 public:
  TrainableModelBuilder(boost::function<double (const Traxel&)> appearance = ConstantFeature(1000),
                        boost::function<double (const Traxel&)> disappearance = ConstantFeature(1000),
                        boost::function<double (const Traxel&, const Traxel&, feature_type)> move = SquaredDistance(),
                        boost::function<double (feature_type)> count = ConstantFeature(100),
                        double forbidden_cost = 1000000,
                        double opportunity_cost = 1000,
                        unsigned max_division_level=0,
                        unsigned max_count=2) :
  ModelBuilder(appearance, disappearance, move, count, forbidden_cost, opportunity_cost, max_division_level, max_count) {}
  virtual boost::shared_ptr<ModelBuilder> clone() const;
  virtual ~TrainableModelBuilder() {}

  virtual boost::shared_ptr<Model> build( const MultiHypothesesGraph& );

 private:
  void add_detection_factors( const MultiHypothesesGraph&, Model&, const MultiHypothesesGraph::Node& ) const;
  void add_outgoing_factors( const MultiHypothesesGraph&, Model&, const MultiHypothesesGraph::Node& ) const;
  void add_incoming_factors( const MultiHypothesesGraph&, Model&, const MultiHypothesesGraph::Node& ) const;

  void add_detection_factor( Model&, const Traxel& ) const;
  void add_outgoing_factor( const MultiHypothesesGraph&,
                            Model&,
                            const MultiHypothesesGraph::Node&,
                            const Traxel&,
                            const std::vector<Traxel>& ) const;
  void add_incoming_factor( const MultiHypothesesGraph&,
                            Model&,
                            const MultiHypothesesGraph::Node&,
                            const Traxel&,
                            const std::vector<Traxel>& ) const;
  

};


class CVPR2014ModelBuilder : public ModelBuilder {
 public:
  CVPR2014ModelBuilder(boost::function<double (const Traxel&)> appearance = ConstantFeature(1000),
                       boost::function<double (const Traxel&)> disappearance = ConstantFeature(1000),
                       boost::function<double (const Traxel&, const Traxel&, feature_type)> move = SquaredDistance(),
                       boost::function<double (feature_type)> count = ConstantFeature(100),
                       double forbidden_cost = 1000000,
                       double opportunity_cost = 1000,
                       unsigned max_division_level=0,
                       unsigned max_count=2) :
  ModelBuilder(appearance, disappearance, move, count, forbidden_cost, opportunity_cost, max_division_level, max_count) {}
  virtual boost::shared_ptr<ModelBuilder> clone() const;
  virtual ~CVPR2014ModelBuilder() {}

  virtual boost::shared_ptr<Model> build( const MultiHypothesesGraph& );

  void add_outgoing_factors( const MultiHypothesesGraph&, Model&, const MultiHypothesesGraph::Node& ) const;
  void add_count_factors( const MultiHypothesesGraph&, Model&);
  void add_count_factors( const MultiHypothesesGraph& hypotheses, Model& m, OpengmLPCplex& cplex) const;

 private:
  void add_detection_factors( const MultiHypothesesGraph&, Model&, const MultiHypothesesGraph::Node& ) const;
  void add_incoming_factors( const MultiHypothesesGraph&, Model&, const MultiHypothesesGraph::Node& );
  void add_detection_factor( Model&, const Traxel& ) const;
  void add_outgoing_factor( const MultiHypothesesGraph&,
                            Model&,
                            const MultiHypothesesGraph::Node&,
                            const Traxel&,
                            const std::vector<Traxel>& ) const;
  void add_incoming_factor( const MultiHypothesesGraph&,
                            Model&,
                            const MultiHypothesesGraph::Node&,
                            const Traxel&,
                            const std::vector<Traxel>& );
  void fill_probabilities(feature_array& probabilities, size_t maximum_active_regions ) const;
  void add_count_factor( Model& m,
                         const std::vector<Traxel>& traxels,
                         size_t maximum_active_regions );
  void add_explicit_count_factor( Model& m,
                                  const std::vector<Traxel>& traxels,
                                  size_t maximum_active_regions ) const;
  void add_hierarchical_count_factor( Model& m,
                                      const std::vector<Traxel>& traxels,
                                      size_t maximum_active_regions );
  void add_count_helper( Model& m,
                         std::vector<std::vector<std::pair<size_t, size_t> > >& lower_level, const size_t maximum_active_regions,
                         std::vector<std::vector<std::pair<size_t, size_t> > >& out );

};


} /* namespace multihypotheses */
} /* namespace pgm */
} /* namespace pgmlink */



#endif /* PGMLINK_PGM_MULTI_HYPOTHESES_H */
