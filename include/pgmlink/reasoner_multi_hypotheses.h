/**
   @file
   @ingroup pgm
   @brief graphical model-based reasoner
*/

#ifndef REASONER_MULTI_HYPOTHESES_H
#define REASONER_PGM_H

#include <algorithm>
#include <iterator>
#include <map>
#include <vector>
#include <utility>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bimap.hpp>
#include <opengm/inference/inference.hxx>
#include <opengm/inference/lpcplex.hxx>

#include "pgmlink/event.h"
#include "pgmlink/feature.h"
#include "pgmlink/pgm.h"
#include "pgmlink/multi_hypotheses_graph.h"
#include "pgmlink/reasoner.h"
#include "pgmlink/pgm_multi_hypotheses.h"

namespace pgmlink {
class Traxel;
namespace pgm {
typedef opengm::LPCplex<OpengmModel, opengm::Minimizer> OpengmLPCplex;
} /* namespace pgm */

class MultiHypotheses { // : public Reasoner {
 public:
  typedef pgm::multihypotheses::Model::trax_var_map trax_var_map;
  typedef pgm::multihypotheses::Model::arc_var_map arc_var_map;

  MultiHypotheses(bool with_constraints = true,
                  double ep_gap = 0.01,
                  bool fixed_detections = false,
                  double cplex_timeout = 1e+75
                  )
      : optimizer_(NULL),
        with_constraints_(with_constraints),
        fixed_detections_(fixed_detections),
        ep_gap_(ep_gap),
        cplex_timeout_(cplex_timeout),
        builder_(new pgm::multihypotheses::TrainableModelBuilder())
  {(*builder_).with_detection_vars().with_divisions(); }
    

  MultiHypotheses(const pgm::multihypotheses::ModelBuilder& builder,
                  bool with_constraints = true,
                  double ep_gap = 0.01,
                  bool fixed_detections = false,
                  double cplex_timeout = 1e+75
                  ) 
      : optimizer_(NULL),
        with_constraints_(with_constraints),
        fixed_detections_(fixed_detections),
        ep_gap_(ep_gap),
        cplex_timeout_(cplex_timeout),
        builder_(builder.clone())
  {};
  ~MultiHypotheses();

  virtual void formulate( const MultiHypothesesGraph& );
  virtual void infer();
  virtual void conclude( MultiHypothesesGraph& );

  double forbidden_cost() const;
  bool with_constraints() const;
  const pgm::multihypotheses::ModelBuilder& builder() { return *builder_; }
  void builder(const pgm::multihypotheses::ModelBuilder& builder) {
    builder_ = builder.clone(); }

  /** Return current state of graphical model
   *
   * The returned pointer may be NULL before formulate() is called
   * the first time.
   **/
  const pgm::OpengmModel* get_graphical_model() const;

  /** Return mapping from HypothesesGraph nodes to graphical model variable ids
   *
   * The map is populated after the first call to formulate().
   */
  const trax_var_map& get_trax_map() const;

  /** Return mapping from HypothesesGraph arcs to graphical model variable ids
   *
   * The map is populated after the first call to formulate().
   */
  const arc_var_map& get_arc_map() const;
    

 private:
  // copy and assingment have to be implemented, yet
  MultiHypotheses(const MultiHypotheses&) {};
  MultiHypotheses& operator=(const MultiHypotheses&) { return *this;};
  void reset();
    
  pgm::OpengmLPCplex* optimizer_;
  shared_ptr<pgm::multihypotheses::Model> linking_model_;

  bool with_constraints_;
  bool fixed_detections_;

  double ep_gap_;
  double cplex_timeout_;
  boost::shared_ptr<pgm::multihypotheses::ModelBuilder> builder_;
};

} /* namespace pgmlink */
#endif /* REASONER_MULTI_HYPOTHESES_H */
