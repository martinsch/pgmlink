/**
   @file
   @ingroup pgm
   @brief graphical model-based reasoner
*/

#ifndef REASONER_OPENGM_H
#define REASONER_OPENGM_H

#include <map>
#include <vector>
#include <utility>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <opengm/inference/inference.hxx>
#include <opengm/inference/lpcplex.hxx>

#include "pgmlink/event.h"
#include "pgmlink/feature.h"
#include "pgmlink/graphical_model.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/reasoner.h"

namespace pgmlink {
  class Traxel;

  namespace pgm {
    typedef opengm::LPCplex<OpengmModel, opengm::Minimizer> OpengmLPCplex;
    using boost::shared_ptr;
    using boost::function;
    using std::map;
    using std::vector;

    /**
       @brief Accessing entries of a Factor/Function that was already added to a graphical model.

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


    /**
       @brief Chaingraph model formulated as an Opengm graphical model.

       Represents an opengm model to solve the matching problem in a
       HypothesesGraph. Use a ChaingraphModelBuilder to construct the
       model. During construction of the ChaingraphModel a random
       variable is added to the graphical model for every node and
       every arc in the HypothesesGraph. The mapping between nodes
       resp. arcs and random variables is stored in the fields
       node_var and arc_var.

       A node in the HypothesesGraph describes a detection in the link
       model. The corresponding random variable determines wether it
       is an actual object or a misdetection. Similarly, an arc is
       interpreted as a possible link between two objects whose state is
       determined by the corresponding random variable.

       @see ChaingraphModelBuilder
       @see HypothesesGraph
     */
    struct ChaingraphModel {
    ChaingraphModel() : opengm_model( new OpengmModel() ) {}
    ChaingraphModel( shared_ptr<OpengmModel> m,
		   map<HypothesesGraph::Node, OpengmModel::IndexType> node_var = map<HypothesesGraph::Node, OpengmModel::IndexType>(),
		   map<HypothesesGraph::Arc, OpengmModel::IndexType> arc_var = map<HypothesesGraph::Arc, OpengmModel::IndexType>()
		   )
    : opengm_model(m), node_var(node_var), arc_var(arc_var) {
      weight_map[det_weight] = vector<OpengmModel::IndexType>();
      weight_map[mov_weight] = vector<OpengmModel::IndexType>();
      weight_map[div_weight] = vector<OpengmModel::IndexType>();
      weight_map[app_weight] = vector<OpengmModel::IndexType>();
      weight_map[dis_weight] = vector<OpengmModel::IndexType>();
      weight_map[opp_weight] = vector<OpengmModel::IndexType>();
    }
      
      shared_ptr<OpengmModel> opengm_model; ///< opengm model usually constructed by ChaingraphModelBuilder
      map<HypothesesGraph::Node, OpengmModel::IndexType> node_var; ///< maps nodes to random variables representing detections
      map<HypothesesGraph::Arc, OpengmModel::IndexType> arc_var; ///< maps arcs to random variables representing links
      
      enum WeightType {det_weight, mov_weight, div_weight, app_weight, dis_weight, opp_weight};
      map<WeightType, vector<OpengmModel::IndexType> > weight_map; ///< associates events with their corresponding weight ids
      
      //void set_weights( WeightType, vector<OpengmModel::ValueType> );
      //const vector<OpengmModel::ValueType>& get_weights( WeightType );
    };

    class ChaingraphModelBuilder {
    public:
    ChaingraphModelBuilder(boost::function<double (const Traxel&)> appearance = ConstantFeature(1000),
    			     boost::function<double (const Traxel&)> disappearance = ConstantFeature(1000),
    			     boost::function<double (const Traxel&, const Traxel&)> move = SquaredDistance(),
    			     double opportunity_cost = 0,
    			     double forbidden_cost = 100000)
	: with_detection_vars_(false),
	with_divisions_(false),
    	appearance_(appearance),
    	disappearance_(disappearance),
    	move_(move),
    	opportunity_cost_(opportunity_cost),
    	forbidden_cost_(forbidden_cost) {}

      virtual ChaingraphModelBuilder* clone() const = 0;

      // mandatory parameters
      function<double (const Traxel&)> appearance() const { return appearance_; }
      ChaingraphModelBuilder& appearance( function<double (const Traxel&)> );

      function<double (const Traxel&)> disappearance() const { return disappearance_; }
      ChaingraphModelBuilder& disappearance( function<double (const Traxel&)> );

      function<double (const Traxel&,const Traxel&)> move() const { return move_; }
      ChaingraphModelBuilder& move( function<double (const Traxel&,const Traxel&)> );

      double opportunity_cost() const { return opportunity_cost_; }
      ChaingraphModelBuilder& opportunity_cost( double c ) { opportunity_cost_ = c; return *this; }

      double forbidden_cost() const { return forbidden_cost_; }
      ChaingraphModelBuilder& forbidden_cost( double c ) { forbidden_cost_ = c; return *this; } 

      //// optional parameters
      // detection vars
      ChaingraphModelBuilder& with_detection_vars( function<double (const Traxel&)> detection=ConstantFeature(10),
						   function<double (const Traxel&)> non_detection=ConstantFeature(200));
      ChaingraphModelBuilder& without_detection_vars();
      bool has_detection_vars() const { return with_detection_vars_; }
      function<double (const Traxel&)> detection() const { return detection_; }
      function<double (const Traxel&)> non_detection() const { return non_detection_; }

      // divisions
      ChaingraphModelBuilder& with_divisions( function<double (const Traxel&,const Traxel&,const Traxel&)> div = KasterDivision(10) );
      ChaingraphModelBuilder& without_divisions();
      bool has_divisions() const { return with_divisions_; }
      function<double (const Traxel&,const Traxel&,const Traxel&)> division() const { return division_; }

      // build
      virtual shared_ptr<ChaingraphModel> build( const HypothesesGraph& ) const = 0;      

      // refinement
      static void add_hard_constraints( const ChaingraphModel&, const HypothesesGraph&, OpengmLPCplex& );
      static void fix_detections( const ChaingraphModel&, const HypothesesGraph&, OpengmLPCplex& );

    protected:
      void add_detection_vars( const HypothesesGraph&, ChaingraphModel& ) const;
      void add_assignment_vars( const HypothesesGraph&, ChaingraphModel& ) const;

    private:
      static void couple( const ChaingraphModel&, const HypothesesGraph::Node&, const HypothesesGraph::Arc&, OpengmLPCplex& );
      
      bool with_detection_vars_;
      bool with_divisions_;

      function<double (const Traxel&)> detection_;
      function<double (const Traxel&)> non_detection_;
      function<double (const Traxel&)> appearance_;
      function<double (const Traxel&)> disappearance_;
      function<double (const Traxel&,const Traxel&)> move_;
      function<double (const Traxel&,const Traxel&,const Traxel&)> division_;
      double opportunity_cost_;
      double forbidden_cost_;
    };

    class TrainableChaingraphModelBuilder : public ChaingraphModelBuilder {
    public:
      TrainableChaingraphModelBuilder(boost::function<double (const Traxel&)> appearance = ConstantFeature(1000),
				      boost::function<double (const Traxel&)> disappearance = ConstantFeature(1000),
				      boost::function<double (const Traxel&, const Traxel&)> move = SquaredDistance(),
				      double opportunity_cost = 0,
				      double forbidden_cost = 100000)
    	: ChaingraphModelBuilder(appearance, disappearance, move, opportunity_cost, forbidden_cost) {}
      virtual TrainableChaingraphModelBuilder* clone() const;

      // build
      virtual shared_ptr<ChaingraphModel> build( const HypothesesGraph& ) const;

    private:
      void add_detection_factor( const HypothesesGraph&, ChaingraphModel&, const HypothesesGraph::Node& ) const;
      void add_outgoing_factor( const HypothesesGraph&, ChaingraphModel&, const HypothesesGraph::Node& ) const;
      void add_incoming_factor( const HypothesesGraph&, ChaingraphModel&, const HypothesesGraph::Node& ) const;
    };

    class ChaingraphModelBuilderECCV12 : public ChaingraphModelBuilder {
    public:
      ChaingraphModelBuilderECCV12(boost::function<double (const Traxel&)> appearance = ConstantFeature(1000),
				   boost::function<double (const Traxel&)> disappearance = ConstantFeature(1000),
				   boost::function<double (const Traxel&, const Traxel&)> move = SquaredDistance(),
				   double opportunity_cost = 0,
				   double forbidden_cost = 100000)
    	: ChaingraphModelBuilder(appearance, disappearance, move, opportunity_cost, forbidden_cost) {}
      virtual ChaingraphModelBuilderECCV12* clone() const;

      // build
      virtual shared_ptr<ChaingraphModel> build( const HypothesesGraph& ) const;

    private:
      void add_detection_factor( const HypothesesGraph&, ChaingraphModel&, const HypothesesGraph::Node& ) const;
      void add_outgoing_factor( const HypothesesGraph&, ChaingraphModel&, const HypothesesGraph::Node& ) const;
      void add_incoming_factor( const HypothesesGraph&, ChaingraphModel&, const HypothesesGraph::Node& ) const;
    };

  } /* namespace pgm */


  class Chaingraph : public Reasoner {
    public:

    Chaingraph(bool with_constraints = true,
	       double ep_gap = 0.01,
	       bool fixed_detections = false
	       )
      : optimizer_(NULL),
      with_constraints_(with_constraints),
      fixed_detections_(fixed_detections),
      ep_gap_(ep_gap),
      builder_(NULL)
	{ builder_ = new pgm::ChaingraphModelBuilderECCV12(); (*builder_).with_detection_vars().with_divisions(); }
    

  Chaingraph(const pgm::ChaingraphModelBuilder& builder,
	     bool with_constraints = true,
	     double ep_gap = 0.01,
	     bool fixed_detections = false
    ) 
    : optimizer_(NULL),
    with_constraints_(with_constraints),
    fixed_detections_(fixed_detections),
    ep_gap_(ep_gap),
    builder_(builder.clone())
    {};
    ~Chaingraph();

    virtual void formulate( const HypothesesGraph& );
    virtual void infer();
    virtual void conclude( HypothesesGraph& );

    double forbidden_cost() const;
    bool with_constraints() const;
    const pgm::ChaingraphModelBuilder& builder() { return *builder_; }
    void builder(const pgm::ChaingraphModelBuilder& builder) {
      if(builder_) delete builder_; builder_ = builder.clone(); }

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
    
    pgm::OpengmLPCplex* optimizer_;
    shared_ptr<pgm::ChaingraphModel> linking_model_;

    bool with_constraints_;
    bool fixed_detections_;

    double ep_gap_;
    pgm::ChaingraphModelBuilder* builder_;
};

} /* namespace pgmlink */
#endif /* REASONER_OPENGM_H */
