#ifndef GRAPHICAL_MODEL_H
#define GRAPHICAL_MODEL_H

#include <vector>
#include <utility>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <opengm/adder.hxx>
#include <opengm/discretespace.hxx>
#include <opengm/graphicalmodel.hxx>
#include <opengm/explicitfactor.hxx>
#include <opengm/inference/lpcplex.hxx>
#include <vigra/multi_array.hxx>

//#include "hypotheses.h"
#include "graph.h"

using boost::vecS;
using boost::bidirectionalS;
using boost::shared_ptr;

namespace Tracking {

class OpengmMrf {
    public:
    typedef double Energy;
    typedef opengm::DiscreteSpace ogmSpace;
    typedef opengm::ExplicitFactor<Energy> ogmFactor;
    typedef opengm::GraphicalModel<ogmFactor, opengm::Adder> ogmGraphicalModel;
    typedef opengm::Minimizer ogmAccumulator;
    typedef opengm::Inference<ogmGraphicalModel, ogmAccumulator> ogmInference;
    
    OpengmMrf();
    ~OpengmMrf();

    ogmSpace* Space() {return space_; }
    ogmGraphicalModel* Model() {return model_; }
 
    private:
    OpengmMrf(const OpengmMrf&);
    OpengmMrf& operator=(const OpengmMrf&);

    ogmSpace* space_;
    ogmGraphicalModel* model_;
};




//////////////////////////////////////
/* Legacy graphical model */
//////////////////////////////////////

struct Var {
    Var() : type(InvalidType), state(0), idx(0), 
	    timestep(0), traxel_idx(0),
	    from_timestep(0), to_timestep(0){}
    enum VarType {DetectionType, TransitionType, InvalidType};
    VarType type;
    // will be set to the inferred state after MAP estimation: 0 or 1
    size_t state;
    // discrete space index
    size_t idx;    

    // detection properties
    size_t timestep;
    unsigned int traxel_label; // 'official' traxel idx 
    size_t traxel_idx; // internal idx

    // transition properties
    size_t from_timestep;
    size_t to_timestep;
};

class TertiaryEnergy;
class BinaryEnergy;
class UnaryEnergy;
class Traxel;
class Event;
class GraphicalModel {
	public:
	/**
	* Construct a graphical model for tracking.
	*
	* @param traxels 2-dim. array: first dimension encodes the timestep,
	*        second dimension contains traxels at a timestep
	*/
	static
	shared_ptr<GraphicalModel>
	from(const std::vector< std::vector<Traxel> >& traxels,
	    const TertiaryEnergy& division,
	    const BinaryEnergy& move,
	    const BinaryEnergy& mismove,
	    const UnaryEnergy& detection,
	    const UnaryEnergy& misdetection
	);

	void infer();
	
	// infer() has to be called first
	// outer index is a real timestep, not a volume index
        // i.e. interpret idx 0 as the events, that happened between volume 0 and 1
	shared_ptr<std::vector< std::vector<Event> > > events();

	void write_graphviz(string filename);

	private:
	GraphicalModel() {};

	bool infer_called_;
	size_t n_timesteps_;

	// opengm
	typedef double Energy;
	typedef opengm::DiscreteSpace ogmSpace;
	typedef opengm::ExplicitFactor<Energy> ogmFactor;
	typedef opengm::GraphicalModel<ogmFactor, opengm::Adder> ogmGraphicalModel;
	typedef opengm::LPCplex<ogmGraphicalModel, opengm::Minimizer> ogmOptimizer;
	
	ogmSpace space_;
	ogmGraphicalModel model_;
	boost::shared_ptr<ogmOptimizer> optimizer_;

	// boost graph
	typedef boost::adjacency_list<vecS, vecS, bidirectionalS, Var> boostGraph;
	typedef boost::graph_traits<boostGraph>::vertex_descriptor boostVertex;
	typedef boost::graph_traits<boostGraph>::edge_descriptor boostEdge;

	//timestep filters	
	template <typename VarMap>
	struct detections_at {
	    detections_at() : t_(0) { }
	detections_at(VarMap m, size_t t) : t_(t), m_(m) { }
	    template <typename Vertex>
	    bool operator()(const Vertex& v) const {
		return m_[v].type == Var::DetectionType && m_[v].timestep == t_ ? true : false;
	    }
	    size_t t_;
	    VarMap m_;
	};
	template <typename VarMap>
	struct transitions_between {
	    transitions_between() : from_t_(0), to_t_(0) { }
	    transitions_between(VarMap m, size_t from, size_t to) : m_(m), from_t_(from), to_t_(to) { }
	    template <typename Vertex>
	    bool operator()(const Vertex& v) const {
		return m_[v].type == Var::TransitionType && m_[v].from_timestep == from_t_ && m_[v].to_timestep == to_t_? true : false;
	    }
	    size_t from_t_, to_t_;
	    VarMap m_;
	};
    
	boostGraph graph_;
	std::map<size_t, std::vector<boostVertex> > detections_at_;


	// implementation
	//(timestep, idx in traxel vector)
	typedef std::pair<size_t, size_t> TraxelId;
	boostVertex add_detection_var(TraxelId tr, unsigned int traxel_label);
	boostVertex add_transition_var(size_t from_t, size_t to_t);
	void add_transition_vars(size_t from_t);
	void add_transition_energies(size_t from_t, 
	    const BinaryEnergy& move,
	    const BinaryEnergy& mismove,
	    const TertiaryEnergy& division,
	    const vector< vector<Traxel> >& traxels);
};

} /* namespace Tracking */


#endif /* GRAPHICAL_MODEL_H */
