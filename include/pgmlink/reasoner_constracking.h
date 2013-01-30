#ifndef CONSTRACKING_REASONER_H
#define CONSTRACKING_REASONER_H

#include <map>
#include <boost/function.hpp>
#include <opengm/inference/inference.hxx>
#include <opengm/inference/lpcplex.hxx>

#include "pgmlink/graphical_model.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/reasoner.h"

namespace pgmlink {
class Traxel;


class ConservationTracking : public Reasoner {
    public:
	ConservationTracking(
				unsigned int max_number_objects,
				boost::function<double (const Traxel&, const size_t)> detection,
			    boost::function<double (const Traxel&, const size_t)> division,
			    boost::function<double (const double)> transition,
                double forbidden_cost = 0,
			    bool with_constraints = true,
			    bool fixed_detections = false,
			    double ep_gap = 0.01,
			    bool with_appearance = false,
			    bool with_disappearance = false,
			    bool with_tracklets = false
    )
    : max_number_objects_(max_number_objects),
    detection_(detection),
    division_(division),
    transition_(transition),
    forbidden_cost_(forbidden_cost),
//    pgm_(NULL),
    optimizer_(NULL),
    with_constraints_(with_constraints),
    fixed_detections_(fixed_detections),
    ep_gap_(ep_gap),
    with_appearance_(with_appearance),
    with_disappearance_(with_disappearance),
    with_tracklets_(with_tracklets),
    number_of_appearance_nodes_(0),
    number_of_disappearance_nodes_(0)
    { };
    ~ConservationTracking();

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
//    const pgm::OpengmModelDeprecated* get_graphical_model() const;

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
    ConservationTracking(const ConservationTracking&) {};
    ConservationTracking& operator=(const ConservationTracking&) { return *this;};

    void reset();
    void add_constraints( const HypothesesGraph& );
    void add_detection_nodes( const HypothesesGraph& );
    void add_appearance_nodes( const HypothesesGraph& );
    void add_disappearance_nodes( const HypothesesGraph& );
    void add_transition_nodes( const HypothesesGraph& );
    void add_division_nodes(const HypothesesGraph& );
    void add_finite_factors( const HypothesesGraph& );

    // helper
    void fix_detections( const HypothesesGraph& );
    size_t cplex_id(size_t opengm_id, size_t state);
//    template<typename table_t, typename const_iter>
//      void add_factor( const table_t& table, const_iter first_idx, const_iter last_idx );


    unsigned int max_number_objects_;

    // energy functions
    boost::function<double (const Traxel&, const size_t)> detection_;
    boost::function<double (const Traxel&, const size_t)> division_;
//    boost::function<double (const Traxel&, const Traxel&, const size_t)> transition_;
    boost::function<double (const double)> transition_;
    double forbidden_cost_;
    
//    pgm::OpengmModelDeprecated* pgm_;
    shared_ptr<pgm::OpengmModelDeprecated> pgm_;
//    pgmlink::pgm::OpengmModelDeprecated::ogmLPCplex* optimizer_;
    opengm::LPCplex<pgm::OpengmModelDeprecated::ogmGraphicalModel, pgm::OpengmModelDeprecated::ogmAccumulator>* optimizer_;

    std::map<HypothesesGraph::Node, size_t> node_map_;
    std::map<HypothesesGraph::Node, size_t> div_node_map_;
    std::map<HypothesesGraph::Node, size_t> app_node_map_;
    std::map<HypothesesGraph::Node, size_t> dis_node_map_;
    std::map<HypothesesGraph::Arc, size_t> arc_map_;

    bool with_constraints_;
    bool fixed_detections_;
    double ep_gap_;

    bool with_appearance_;
    bool with_disappearance_;

    bool with_tracklets_;

    unsigned int number_of_detection_nodes_, number_of_transition_nodes_, number_of_division_nodes_;
    unsigned int number_of_appearance_nodes_, number_of_disappearance_nodes_;

    HypothesesGraph tracklet_graph_;
    std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > tracklet2traxel_node_map_;
};



/******************/
/* Implementation */
/******************/
 
// template< typename table_t, typename const_iter >
//   void ConservationTracking::add_factor( const table_t& table, const_iter first_idx, const_iter last_idx ){
//   OpengmModelDeprecated::FunctionIdentifier id=pgm_->Model()->addFunction(table);
//   pgm_->Model()->addFactor(id, first_idx, last_idx);
// }
 
} /* namespace pgmlink */
#endif /* MRF_REASONER_H */
  
