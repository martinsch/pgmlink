/**
   @file
   @ingroup pgm
   @brief graphical model-based reasoner
*/

#ifndef MRF_REASONER_H
#define MRF_REASONER_H

#include <map>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include "pgmlink/graphical_model.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/reasoner.h"

namespace pgmlink {
class Traxel;

 /* class ChaingraphModelBuilder { */
 /*   ChaingraphModelBuilder( */
 /* 	       double opportunity_cost = 0, */
 /* 	       double forbidden_cost = 0 */
 /* 			  ) */

 /*   ChaingraphModelBuilder& opportunity_cost(double); */
 /*   ChaingraphModelBuilder& forbidden_cost(double); */
 /*   shared_ptr<GraphicalModelType> build(); */

 /* private: */
 /*   double opportunity_cost_; */
 /*   double forbidden_cost_; */
 /* }; */

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
    template<typename table_t, typename const_iter>
      void add_factor( const table_t& table, const_iter first_idx, const_iter last_idx );

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



/******************/
/* Implementation */
/******************/
 
 template< typename table_t, typename const_iter >
   void Chaingraph::add_factor( const table_t& table, const_iter first_idx, const_iter last_idx ){
   OpengmModel<>::FunctionIdentifier id=mrf_->addFunction(table);
   mrf_->addFactor(id, first_idx, last_idx);
 }
 
} /* namespace pgmlink */
#endif /* MRF_REASONER_H */
