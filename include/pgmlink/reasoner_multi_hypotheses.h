#ifndef REASONER_MULTI_HYPOTHESES_H
#define REASONER_MULTI_HYPOTHESES_H

// stl

// boost

// opengm

// pgmlink
#include <pgmlink/pgm.h>
#include <pgmlink/reasoner.h>
#include <pgmlink/multi_hyptheses_graph.h>


namespace pgmlink {
class MultiHypothesesTracking : public Reasoner {
      public:
	MultiHypothesesTracking(
            unsigned int max_number_objects,
            boost::function<double (const Traxel&, const size_t)> detection,
            boost::function<double (const Traxel&, const size_t)> division,
            boost::function<double (const double)> transition,
            double forbidden_cost = 0,
            double ep_gap = 0.01,
            bool with_tracklets = false,
            bool with_divisions = true,
            boost::function<double (const Traxel&)> disappearance_cost_fn = ConstantFeature(500.0),
            boost::function<double (const Traxel&)> appearance_cost_fn = ConstantFeature(500.0),
            bool with_misdetections_allowed = true,
            bool with_appearance = true,
            bool with_disappearance = true,
            double transition_parameter = 5,
            bool with_constraints = true
                             )
        : max_number_objects_(max_number_objects),
          detection_(detection),
          division_(division),
          transition_(transition),
          forbidden_cost_(forbidden_cost),
          optimizer_(NULL),
          ep_gap_(ep_gap),
          with_tracklets_(with_tracklets),
          with_divisions_(with_divisions),
          disappearance_cost_(disappearance_cost_fn),
          appearance_cost_(appearance_cost_fn),
          number_of_transition_nodes_(0), 
          number_of_division_nodes_(0),
          number_of_appearance_nodes_(0),
          number_of_disappearance_nodes_(0),
          with_misdetections_allowed_(with_misdetections_allowed),
          with_appearance_(with_appearance),
          with_disappearance_(with_disappearance),
          transition_parameter_(transition_parameter),
          with_constraints_(with_constraints)
    { };
    ~MultiHypothesesTracking();

    virtual void formulate( const MultiHypothesesGraph& );
    virtual void infer();
    virtual void conclude( MultiHypothesesGraph& );

    double forbidden_cost() const;
    bool with_constraints() const;

    /** Return current state of graphical model
     *
     * The returned pointer may be NULL before formulate() is called
     * the first time.
     **/
//    const pgm::OpengmModelDeprecated* get_graphical_model() const;

    /** Return mapping from MultiHypothesesGraph nodes to graphical model variable ids
     *
     * The map is populated after the first call to formulate().
     */
//    const std::map<MultiHypothesesGraph::Node, size_t>& get_node_map() const;

    /** Return mapping from MultiHypothesesGraph arcs to graphical model variable ids
     *
     * The map is populated after the first call to formulate().
     */
    const std::map<MultiHypothesesGraph::Arc, size_t>& get_arc_map() const;
    

    private:
    // copy and assingment have to be implemented, yet
    MultiHypothesesTracking(const MultiHypothesesTracking&) {};
    MultiHypothesesTracking& operator=(const MultiHypothesesTracking&) { return *this;};

    void reset();
    void add_constraints( const MultiHypothesesGraph& );
    void add_detection_nodes( const MultiHypothesesGraph& );
    void add_appearance_nodes( const MultiHypothesesGraph& );
    void add_disappearance_nodes( const MultiHypothesesGraph& );
    void add_transition_nodes( const MultiHypothesesGraph& );
    void add_division_nodes(const MultiHypothesesGraph& );
    void add_finite_factors( const MultiHypothesesGraph& );

    // helper
    size_t cplex_id(size_t opengm_id, size_t state);


    unsigned int max_number_objects_;

    // energy functions
    boost::function<double (const Traxel&, const size_t)> detection_;
    boost::function<double (const Traxel&, const size_t)> division_;
    boost::function<double (const double)> transition_;

    double forbidden_cost_;
    
    shared_ptr<pgm::OpengmModelDeprecated> pgm_;
    opengm::LPCplex<pgm::OpengmModelDeprecated::ogmGraphicalModel, pgm::OpengmModelDeprecated::ogmAccumulator>* optimizer_;

    std::map<MultiHypothesesGraph::Node, size_t> div_node_map_;
    std::map<MultiHypothesesGraph::Node, size_t> app_node_map_;
    std::map<MultiHypothesesGraph::Node, size_t> dis_node_map_;
    std::map<MultiHypothesesGraph::Arc, size_t> arc_map_;

    double ep_gap_;

    bool with_tracklets_, with_divisions_;

    boost::function<double (const Traxel&)> disappearance_cost_;
    boost::function<double (const Traxel&)> appearance_cost_;

    unsigned int number_of_transition_nodes_, number_of_division_nodes_;
    unsigned int number_of_appearance_nodes_, number_of_disappearance_nodes_;

    bool with_misdetections_allowed_;
    bool with_appearance_;
    bool with_disappearance_;

    double transition_parameter_;

    bool with_constraints_;

    MultiHypothesesGraph tracklet_graph_;
    std::map<MultiHypothesesGraph::Node, std::vector<MultiHypothesesGraph::Node> > tracklet2traxel_node_map_;
};
}
}

#endif /* REASONER_MULTI_HYPOTHESES_H */
