/**
   @file
   @ingroup pgm
   @brief graphical model-based reasoner
*/

#ifndef MRF_REASONER_H
#define MRF_REASONER_H

#include <map>
#include <vector>
#include <utility>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/inference/lpcplex.hxx>
#include <opengm/operations/adder.hxx>
#include <opengm/graphicalmodel/loglinearmodel.hxx>

#include "pgmlink/event.h"
#include "pgmlink/graphical_model.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/reasoner.h"

namespace pgmlink {
  class Traxel;

  namespace pgm {
    typedef opengm::GraphicalModel<double, opengm::Adder> OpengmModel;
    typedef opengm::LPCplex<OpengmModel, opengm::Adder> OpengmLPCplex;
    using boost::shared_ptr;
    using boost::function;
    using std::map;
    using std::vector;

    /**
       \brief Accessing entries of a Factor/Function that was already added to a graphical model.

       Manages a pointer to an element of an array-like opengm function (usually opengm::ExplicitFunction).
       Validity of the pointer is ensured by owning a smart pointer to the full model.

       Use this class to modify factor elements of an already instantiated opengm graphical model.
     */
    class FactorEntry {
    public:
    FactorEntry() : entry_(NULL) {}
    FactorEntry( shared_ptr<OpengmModel> m, /**< has to be valid */
		 OpengmModel::ValueType* entry /**< has to point into the opengm model to ensure the same lifetime */
		 ) :
      m_(m), entry_(entry) {}
      
      void set( OpengmModel::ValueType );
      OpengmModel::ValueType get() const;

      shared_ptr<OpengmModel> model() const { return m_; }

    private:
      shared_ptr<OpengmModel> m_;
      OpengmModel::ValueType* entry_;
    };



    struct TrackingModel {
    TrackingModel()
    : opengm_model( new OpengmModel() ),
	hypotheses( new HypothesesGraph() ) {}
    TrackingModel( shared_ptr<OpengmModel> m,
		   shared_ptr<HypothesesGraph> g,
		   vector<std::pair<Event, FactorEntry> > events = vector<std::pair<Event, FactorEntry> >(),
		   map<HypothesesGraph::Node, OpengmModel::IndexType> node_var = map<HypothesesGraph::Node, OpengmModel::IndexType>(),
		   map<HypothesesGraph::Arc, OpengmModel::IndexType> arc_var = map<HypothesesGraph::Arc, OpengmModel::IndexType>()
		   )
    : opengm_model(m), hypotheses(g), events(events), node_var(node_var), arc_var(arc_var) {}
      
      shared_ptr<OpengmModel> opengm_model;
      shared_ptr<HypothesesGraph> hypotheses;

      vector<std::pair<Event, FactorEntry> > events;
      map<HypothesesGraph::Node, OpengmModel::IndexType> node_var;
      map<HypothesesGraph::Arc, OpengmModel::IndexType> arc_var;
    };



    class ChaingraphModelBuilder {
    public:
      ChaingraphModelBuilder(shared_ptr<HypothesesGraph> g,
    			     boost::function<double (const Traxel&)> appearance,
    			     boost::function<double (const Traxel&)> disappearance,
    			     boost::function<double (const Traxel&, const Traxel&)> move,
    			     double opportunity_cost = 0,
    			     double forbidden_cost = 100000)
    	: hypotheses_(g),
	with_detection_vars_(false),
	with_divisions_(false),
	with_hard_constraints_(false),
    	appearance_(appearance),
    	disappearance_(disappearance),
    	move_(move),
    	opportunity_cost_(opportunity_cost),
    	forbidden_cost_(forbidden_cost) {}

      // mandatory parameters
      shared_ptr<HypothesesGraph> hypotheses() { return hypotheses_; }
      ChaingraphModelBuilder& hypotheses( shared_ptr<HypothesesGraph> );

      function<double (const Traxel&)> appearance() { return appearance_; }
      ChaingraphModelBuilder& appearance( function<double (const Traxel&)> );

      function<double (const Traxel&)> disappearance() { return disappearance_; }
      ChaingraphModelBuilder& disappearance( function<double (const Traxel&)> );

      function<double (const Traxel&,const Traxel&)> move() { return move_; }
      ChaingraphModelBuilder& move( function<double (const Traxel&,const Traxel&)> );

      double opportunity_cost() { return opportunity_cost_; }
      ChaingraphModelBuilder& opportunity_cost( double c ) { opportunity_cost_ = c; return *this; }

      double forbidden_cost() { return forbidden_cost_; }
      ChaingraphModelBuilder& forbidden_cost( double c ) { forbidden_cost_ = c; return *this; } 

      //// optional parameters
      // detection vars
      ChaingraphModelBuilder& with_detection_vars( function<double (const Traxel&)> detection,
						   function<double (const Traxel&)> non_detection );
      ChaingraphModelBuilder& without_detection_vars();
      bool with_detection_vars() { return with_detection_vars_; }
      function<double (const Traxel&)> detection() { return detection_; }
      function<double (const Traxel&)> non_detection() { return non_detection_; }

      // divisions
      ChaingraphModelBuilder& with_divisions( function<double (const Traxel&,const Traxel&,const Traxel&)> );
      ChaingraphModelBuilder& without_divisions();
      bool with_divisions() { return with_divisions_; }
      function<double (const Traxel&,const Traxel&,const Traxel&)> division() { return division_; }

      // hard constraints
      ChaingraphModelBuilder& with_hard_constraints( shared_ptr<OpengmLPCplex> );
      ChaingraphModelBuilder& without_hard_constraints();
      bool with_hard_constraints() { return with_hard_constraints_; }
      shared_ptr<OpengmLPCplex> cplex() { return cplex_; }

      // build
      shared_ptr<TrackingModel> build();

    private:
      void add_detection_vars( TrackingModel& ) const;
      void add_assignment_vars( TrackingModel& ) const;
      void add_detection_factor( TrackingModel&, const HypothesesGraph::Node& ) const;
      void add_outgoing_factor( TrackingModel&,
				const HypothesesGraph::Node&,
				boost::function<double (const Traxel&)> disappearance,
				boost::function<double (const Traxel&, const Traxel&)> move,
				boost::function<double (const Traxel&, const Traxel&, const Traxel&)> division,
				double opportunity_cost = 0,
				double forbidden_cost = 10000000) const;

      void add_incoming_factor( TrackingModel&,
				const HypothesesGraph::Node&,
				boost::function<double (const Traxel&)> appearance,
				double forbidden_cost = 10000000) const;

      void add_constraints( TrackingModel& ) const;
      void couple( TrackingModel&, HypothesesGraph::Node&, HypothesesGraph::Arc& ) const;

      shared_ptr<HypothesesGraph> hypotheses_;
      shared_ptr<OpengmLPCplex> cplex_;
      
      bool with_detection_vars_;
      bool with_divisions_;
      bool with_hard_constraints_;

      function<double (const Traxel&)> detection_;
      function<double (const Traxel&)> non_detection_;
      function<double (const Traxel&)> appearance_;
      function<double (const Traxel&)> disappearance_;
      function<double (const Traxel&,const Traxel&)> move_;
      function<double (const Traxel&,const Traxel&,const Traxel&)> division_;
      double opportunity_cost_;
      double forbidden_cost_;
    };
  } /* namespace gm */



  class Chaingraph : public Reasoner {
    public:
    Chaingraph(boost::function<double (const Traxel&)> detection,
	       boost::function<double (const Traxel&)> non_detection,
	       boost::function<double (const Traxel&)> appearance,
	       boost::function<double (const Traxel&)> disappearance,
	       boost::function<double (const Traxel&, const Traxel&)> move,
	       boost::function<double (const Traxel&, const Traxel&, const Traxel&)> division,
	       double opportunity_cost = 0,
	       double forbidden_cost = 0,
	       bool with_constraints = true,
	       bool fixed_detections = false,
	       double ep_gap = 0.01
    ) 
    : detection_(detection), 
    non_detection_(non_detection),
    appearance_(appearance),
    disappearance_(disappearance),
    move_(move),
    division_(division),
    opportunity_cost_(opportunity_cost),
    forbidden_cost_(forbidden_cost),
    mrf_(NULL),
    optimizer_(NULL),
    with_constraints_(with_constraints),
    fixed_detections_(fixed_detections),
    ep_gap_(ep_gap)
    { };
    ~Chaingraph();

    virtual void formulate( const HypothesesGraph& );
    virtual void infer();
    virtual void conclude( HypothesesGraph& );

    double forbidden_cost() const;
    bool with_constraints() const;

    /** Return current state of graphical model
     *
     * The returned pointer may be NULL before formulate() is called
     * the first time.
     **/
    const OpengmModel<>::ogmGraphicalModel* get_graphical_model() const;

    /** Return mapping from HypothesesGraph nodes to graphical model variable ids
     *
     * The map is populated after the first call to formulate().
     */
    const std::map<HypothesesGraph::Node, size_t>& get_node_map() const;

    /** Return mapping from HypothesesGraph arcs to graphical model variable ids
     *
     * The map is populated after the first call to formulate().
     */
    const std::map<HypothesesGraph::Arc, size_t>& get_arc_map() const;
    

    private:
    // copy and assingment have to be implemented, yet
    Chaingraph(const Chaingraph&) {};
    Chaingraph& operator=(const Chaingraph&) { return *this;};

    void reset();
    void add_constraints( const HypothesesGraph& );
    void add_detection_nodes( const HypothesesGraph& );
    void add_transition_nodes( const HypothesesGraph& );
    void add_finite_factors( const HypothesesGraph& );

    // helper
    void couple( HypothesesGraph::Node&, HypothesesGraph::Arc& );
    void fix_detections( const HypothesesGraph&, size_t value );
    void add_outgoing_factor( const HypothesesGraph&, const HypothesesGraph::Node& );
    void add_incoming_factor( const HypothesesGraph&, const HypothesesGraph::Node& );

    // energy functions
    boost::function<double (const Traxel&)> detection_;
    boost::function<double (const Traxel&)> non_detection_;
    boost::function<double (const Traxel&)> appearance_;
    boost::function<double (const Traxel&)> disappearance_;
    boost::function<double (const Traxel&, const Traxel&)> move_;
    boost::function<double (const Traxel&, const Traxel&, const Traxel&)> division_;
    double opportunity_cost_;
    double forbidden_cost_;
    
    OpengmModel<>::ogmGraphicalModel* mrf_;
    OpengmModel<>::ogmInference* optimizer_;

    std::map<HypothesesGraph::Node, size_t> node_map_;
    std::map<HypothesesGraph::Arc, size_t> arc_map_;

    bool with_constraints_;
    bool fixed_detections_;

    double ep_gap_;
};

} /* namespace pgmlink */
#endif /* MRF_REASONER_H */
