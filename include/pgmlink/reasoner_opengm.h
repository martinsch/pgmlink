/**
   @file
   @ingroup pgm
   @brief graphical model-based reasoner
*/

#ifndef REASONER_OPENGM_H
#define REASONER_OPENGM_H

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

    class ChaingraphModelBuilder;
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
    class ChaingraphModel {
    public:
      typedef HypothesesGraph::Node node_t;
      typedef HypothesesGraph::Arc arc_t;
      typedef OpengmModel::IndexType var_t;
      typedef boost::bimap<node_t, var_t>::left_map node_var_map;
      typedef boost::bimap<node_t, var_t>::right_map var_node_map;
      typedef boost::bimap<arc_t, var_t>::left_map arc_var_map;
      typedef boost::bimap<arc_t, var_t>::right_map var_arc_map;

      ChaingraphModel();
      ChaingraphModel( shared_ptr<OpengmModel>,
		       const node_var_map&,
		       const arc_var_map&
		       );
      
      shared_ptr<OpengmModel> opengm_model; ///< opengm model usually constructed by ChaingraphModelBuilder

      const node_var_map& var_of_node() const; ///< maps nodes to random variables representing detections
      const var_node_map& node_of_var() const;
      const arc_var_map& var_of_arc() const; ///< maps arcs to random variables representing links
      const var_arc_map& arc_of_var() const;

      var_t var_of_node(node_t) const;
      var_t var_of_arc(arc_t) const;
      node_t node_of_var(var_t) const;
      arc_t arc_of_var(var_t) const;

      enum VarCategory {node_var, arc_var};
      VarCategory var_category(var_t) const;
      
      enum WeightType {det_weight, mov_weight, div_weight, app_weight, dis_weight, opp_weight};
      map<WeightType, vector<OpengmModel::IndexType> > weight_map; ///< associates events with their corresponding weight ids
      
      //void set_weights( WeightType, vector<OpengmModel::ValueType> );
      //const vector<OpengmModel::ValueType>& get_weights( WeightType );
    private:
      friend class ChaingraphModelBuilder;

      void init();

      boost::bimap<node_t, var_t> node_var_;
      boost::bimap<arc_t, var_t> arc_var_;
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
      virtual ChaingraphModel* build( const HypothesesGraph& ) const = 0;      

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
      virtual ChaingraphModel* build( const HypothesesGraph& ) const;

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
      virtual ChaingraphModel* build( const HypothesesGraph& ) const;

    private:
      void add_detection_factor( const HypothesesGraph&, ChaingraphModel&, const HypothesesGraph::Node& ) const;
      void add_outgoing_factor( const HypothesesGraph&, ChaingraphModel&, const HypothesesGraph::Node& ) const;
      void add_incoming_factor( const HypothesesGraph&, ChaingraphModel&, const HypothesesGraph::Node& ) const;
    };

    class ChaingraphModelTrainer {
    public:
      template<class IT1, class IT2, class IT3>
	std::vector<OpengmModel::ValueType> train(IT1 samples_begin, IT1 samples_end, IT2 node_labels, IT3 arc_labels) const;
    };

  } /* namespace pgm */


  class Chaingraph : public Reasoner {
    public:
    typedef pgm::ChaingraphModel::node_var_map node_var_map;
    typedef pgm::ChaingraphModel::arc_var_map arc_var_map;

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
    const node_var_map& get_node_map() const;

    /** Return mapping from HypothesesGraph arcs to graphical model variable ids
     *
     * The map is populated after the first call to formulate().
     */
    const arc_var_map& get_arc_map() const;
    

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



  /**/
  /* implementation */
  /**/
  
  template<class IT1, class IT2, class IT3>
    std::vector<pgm::OpengmModel::ValueType> pgm::ChaingraphModelTrainer::train(IT1 samples_begin, IT1 samples_end, IT2 node_labels, IT3 arc_labels) const {
    // for each sample: build chaingraph model
    

    // convert HypothesesGraph labels to OpengmModel labels
    
    return std::vector<pgm::OpengmModel::ValueType>();
  }


} /* namespace pgmlink */
#endif /* REASONER_OPENGM_H */
