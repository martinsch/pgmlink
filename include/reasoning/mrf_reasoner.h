#ifndef MRF_REASONER_H
#define MRF_REASONER_H

#include <map>
#include <boost/function.hpp>
#include "graphical_model.h"
#include "hypotheses.h"
#include "reasoning/reasoner.h"

namespace Tracking {
class Traxel;

class SingleTimestepTraxelMrf : public Reasoner {
    public:
    SingleTimestepTraxelMrf(boost::function<double (const Traxel&)> detection,
			    boost::function<double (const Traxel&)> non_detection,
			    boost::function<double (const Traxel&)> appearance,
			    boost::function<double (const Traxel&)> disappearance,
			    boost::function<double (const Traxel&, const Traxel&)> move,
			    boost::function<double (const Traxel&, const Traxel&, const Traxel&)> division,
			    double opportunity_cost = 0,
			    bool with_constraints = true,
			    bool constraints_as_infinite_energy = false
    ) 
    : detection_(detection), 
    non_detection_(non_detection),
    appearance_(appearance),
    disappearance_(disappearance),
    move_(move),
    division_(division),
    opportunity_cost_(opportunity_cost),
    mrf_(NULL),
    optimizer_(NULL),
    with_constraints_(with_constraints),
    constraints_as_infinite_energy_(constraints_as_infinite_energy)
    { };
    ~SingleTimestepTraxelMrf();

    virtual void formulate( const HypothesesGraph& );
    virtual void infer();
    virtual void conclude( HypothesesGraph& );

    /** Return current state of graphical model
     *
     * The returned pointer may be NULL before formulate() is called
     * the first time.
     **/
    const OpengmMrf* get_graphical_model() const;

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
    SingleTimestepTraxelMrf(const SingleTimestepTraxelMrf&) {};
    SingleTimestepTraxelMrf& operator=(const SingleTimestepTraxelMrf&) { return *this;};

    void reset();
    void add_constraints( const HypothesesGraph& );
    void add_detection_nodes( const HypothesesGraph& );
    void add_transition_nodes( const HypothesesGraph& );
    void add_finite_factors( const HypothesesGraph& );

    // helper
    void couple(HypothesesGraph::Node&, HypothesesGraph::Arc&);
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
    
    OpengmMrf* mrf_;
    OpengmMrf::ogmInference* optimizer_;

    std::map<HypothesesGraph::Node, size_t> node_map_;
    std::map<HypothesesGraph::Arc, size_t> arc_map_;

    bool with_constraints_;
    bool constraints_as_infinite_energy_;
};

} /* namespace Tracking */
#endif /* MRF_REASONER_H */
