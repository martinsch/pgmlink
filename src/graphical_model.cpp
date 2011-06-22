#include <assert.h>
#include <ostream>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <utility>
#include <map>
#include <stdexcept>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graphviz.hpp>

#include "energy.h"
#include "graphical_model.h"
#include "event.h"
#include "traxels.h"

using namespace std;
using namespace boost;

namespace Tracking {
////
//// class OpengmMrf
////
    OpengmMrf::OpengmMrf() {
	space_ = new ogmSpace();
	model_ = new ogmGraphicalModel();
    }

    OpengmMrf::~OpengmMrf() {
	delete space_;
	delete model_;
    }

////
//// class GraphicalModel
////
    boost::shared_ptr<GraphicalModel>
    GraphicalModel::from(const vector< vector<Traxel> >& traxels,
    	    const TertiaryEnergy& division,
	    const BinaryEnergy& move,
	    const BinaryEnergy& mismove,
	    const UnaryEnergy& detection,
	    const UnaryEnergy& misdetection) {
	shared_ptr<GraphicalModel> ret(new GraphicalModel);
	ret->infer_called_ = false;
	ret->n_timesteps_ = traxels.size();	

	// iterate through the timesteps and add detection variables and energies
	for(size_t t = 0; t < (ret->n_timesteps_); ++t) {
	    // collect detections at a given timestep for easy reference later
	    vector<boostVertex> detections;
	    for(size_t curr_idx = 0; curr_idx < traxels[t].size(); ++curr_idx) {
	      boostVertex v = ret->add_detection_var(make_pair(t, curr_idx), traxels[t][curr_idx].Id);
		detections.push_back(v);

		size_t vi[] = {ret->graph_[v].idx};
		GraphicalModel::ogmFactor f(ret->space_, vi, vi+1);

		// assign detection energies
		Traxels curr = traxel_map_from_traxel_sequence(traxels[t].begin(), traxels[t].end());
		Traxels prev;
		if(!(t == 0)) {
		    prev = traxel_map_from_traxel_sequence(traxels[t-1].begin(), traxels[t-1].end());
		};
		
		f(0) = misdetection(traxels[t][curr_idx], prev, curr);
		f(1) = detection(traxels[t][curr_idx], prev, curr);
		ret->model_.addFactor(f);
	    }
	    ret->detections_at_[t] = detections;
	}
	
	// add transition vars
	for(size_t t = 0; t < (ret->n_timesteps_ - 1); ++t) {
	    ret->add_transition_vars(t);
	}

	// add transition energies
	for(size_t t = 0; t < (ret->n_timesteps_ - 1); ++t) {
	  ret->add_transition_energies(t, move, mismove, division, traxels);
	}

	// add hard constraints
	ogmOptimizer::Parameter param;
	param.verbose_ = false;
	param.integerConstraint_ = true; // FIXME: investigate meaning of other parameters
	ret->optimizer_ = shared_ptr<ogmOptimizer>(new ogmOptimizer(ret->model_, param));

	//wade through outgoing transitions
	typedef graph_traits<boostGraph>::adjacency_iterator V_it;
	typedef pair<V_it, V_it> Vp; 

	for(size_t t = 0; t < ret->n_timesteps_ - 1; ++t) {
	    for(size_t i = 0; i < ret->detections_at_[t].size(); ++i) {
		boostVertex detection = ret->detections_at_[t][i];
		Vp vp = adjacent_vertices(detection, ret->graph_);
		vector<size_t> cplex_idxs;
		for(; vp.first != vp.second; ++vp.first) {
		    cplex_idxs.push_back(2*ret->graph_[*vp.first].idx+1);
		    // couple detection and transition var
		    vector<size_t> temp;
		    temp.push_back((2*ret->graph_[detection].idx)+1);
		    temp.push_back(cplex_idxs.back());
		    vector<int> temp_coeff;
		    temp_coeff.push_back(1);
		    temp_coeff.push_back(-1);
		    ret->optimizer_->addConstraint(temp.begin(), temp.end(), temp_coeff.begin(), 0, 1);		    
		}
		vector<int> coefficients(cplex_idxs.size(), 1);
		ret->optimizer_->addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coefficients.begin(), 0, 2);
	    }
	}

	//wade through incoming transitions
	typedef boostGraph::inv_adjacency_iterator V_it_inv;
	typedef pair<V_it_inv, V_it_inv> Vp_inv; 


	for(size_t t = 1; t < ret->n_timesteps_; ++t) {
	    for(size_t i = 0; i < ret->detections_at_[t].size(); ++i) {
		boostVertex detection = ret->detections_at_[t][i];
		Vp_inv vp = inv_adjacent_vertices(detection, ret->graph_);
		vector<size_t> cplex_idxs;
		for(; vp.first != vp.second; ++vp.first) {
		    cplex_idxs.push_back(2*ret->graph_[*vp.first].idx + 1);
		    // couple detection and transition var
		    vector<size_t> temp;
		    temp.push_back(2*ret->graph_[detection].idx + 1);
		    temp.push_back(cplex_idxs.back());
		    vector<int> temp_coeff;
		    temp_coeff.push_back(1);
		    temp_coeff.push_back(-1);
		    ret->optimizer_->addConstraint(temp.begin(), temp.end(), temp_coeff.begin(), 0, 1);
		}
		vector<int> coefficients(cplex_idxs.size(), 1);
		ret->optimizer_->addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coefficients.begin(), 0, 1);
	    }
	}

	return ret;
}



void GraphicalModel::infer() {
    opengm::InferenceTermination status = optimizer_->infer();
    if(status != opengm::NORMAL) {
	throw runtime_error("GraphicalModel::infer(): optimizer terminated unnormally");
    }
    infer_called_ = true;
    
    vector<ogmOptimizer::state_type> solution;
    status = optimizer_->arg(solution);
    if(status != opengm::NORMAL) {
	throw runtime_error("GraphicalModel::infer(): solution extraction terminated unnormally");
    }

    // travese vertices and set state
    typedef boost::property_map<boostGraph, vertex_bundle_t>::type VarMap;
    typedef graph_traits<boostGraph>::vertex_iterator v_it;
    typedef pair<v_it, v_it> Vp;

    Vp vp = vertices(graph_);
    for(;vp.first != vp.second; ++vp.first) {
	if(graph_[*vp.first].type == Var::DetectionType) {
	//cout << "Detection "<< graph_[*vp.first].idx << " set to " << solution[graph_[*vp.first].idx] << endl;
	} else {
	//cout << "Transition "<< graph_[*vp.first].idx << " set to " << solution[graph_[*vp.first].idx] << endl;
	}
	graph_[*vp.first].state = solution[graph_[*vp.first].idx];
    }
}



boost::shared_ptr<std::vector< std::vector<Event> > > GraphicalModel::events() {
    if(!infer_called_) {
	throw runtime_error("GraphicalModel::events(): infer() wasn't called yet");
    }
    shared_ptr<std::vector< std::vector<Event> > > ret(new vector< vector<Event> >);

    typedef graph_traits<boostGraph>::adjacency_iterator V_it;
    typedef pair<V_it, V_it> Vp;
    typedef boostGraph::inv_adjacency_iterator V_it_inv;
    typedef pair<V_it_inv, V_it_inv> Vp_inv;

    for(size_t t = 0; t < n_timesteps_ - 1; ++t) {
	ret->push_back(vector<Event>());
	for(size_t i = 0; i < detections_at_[t].size(); ++i) {
	  // forward events
	    boostVertex detection = detections_at_[t][i];
	    Vp vp = adjacent_vertices(detection, graph_);
	    // find active transition vars
	    vector<boostVertex> active;
	    for(; vp.first != vp.second; ++vp.first) {
		if(graph_[*vp.first].state == 1) {
		    active.push_back(*vp.first);
		}
	    }

	    switch(active.size()) {
		// Disappearance
		case 0: {
		    Event e;
		    e.type = Event::Disappearance;
		    e.traxel_ids.push_back(graph_[detection].traxel_label);
		    (*ret)[t].push_back(e);
		    break;
		    }
		// Move
		case 1: {
		    Event e;
		    e.type = Event::Move;
		    e.traxel_ids.push_back(graph_[detection].traxel_label);
		    Vp vp = adjacent_vertices(active[0], graph_);
		    e.traxel_ids.push_back(graph_[*vp.first].traxel_label);
		    assert(distance(vp.first,vp.second) == 1);
		    (*ret)[t].push_back(e);
		    break;
		    }
		// Division
		case 2: {
		    Event e;
		    e.type = Event::Division;
		    e.traxel_ids.push_back(graph_[detection].traxel_label);
		    Vp vp = adjacent_vertices(active[0], graph_);
		    e.traxel_ids.push_back(graph_[*vp.first].traxel_label);
		    assert(distance(vp.first,vp.second) == 1);		    
   		    vp = adjacent_vertices(active[1], graph_);
		    e.traxel_ids.push_back(graph_[*vp.first].traxel_label);
		    assert(distance(vp.first,vp.second) == 1);		   
		    (*ret)[t].push_back(e);
		break;
	        }
		default:
		    throw runtime_error("GraphicalModel::events(): encountered contradicting solution states");
		break;
	    }
	}

	// backward events
	for(size_t i = 0; i < detections_at_[t+1].size(); ++i) {
	    boostVertex detection = detections_at_[t+1][i];
	    Vp_inv vp = inv_adjacent_vertices(detection, graph_);
	    // find active transition vars
	    vector<boostVertex> active;
	    for(; vp.first != vp.second; ++vp.first) {
		if(graph_[*vp.first].state == 1) {
		    active.push_back(*vp.first);
		}
	    }
	    // if all inactive, add appearance
	    if(active.empty()) {
		    Event e;
		    e.type = Event::Appearance;
		    e.traxel_ids.push_back(graph_[detection].traxel_label);
		    (*ret)[t].push_back(e);	      
	    }
	}
    }
    return ret;
}

void GraphicalModel::write_graphviz(string filename) {
    ofstream outfile;
    outfile.open(filename.c_str());
    boost::write_graphviz(outfile, graph_);
    outfile.close();
}



  GraphicalModel::boostVertex GraphicalModel::add_detection_var(TraxelId tr, unsigned int traxel_label) {
    space_.addDimension(2);
    size_t idx = space_.dimension() - 1;

    boostVertex v = add_vertex(graph_);
    graph_[v].type = Var::DetectionType;
    graph_[v].idx = idx;
    graph_[v].timestep = tr.first;
    graph_[v].traxel_idx = tr.second;
    graph_[v].traxel_label = traxel_label;
    return v;
}

GraphicalModel::boostVertex GraphicalModel::add_transition_var(size_t from_t, size_t to_t) {
    space_.addDimension(2);
    size_t idx =  space_.dimension() - 1;

    boostVertex v = add_vertex(graph_);
    graph_[v].type = Var::TransitionType;
    graph_[v].idx = idx;
    graph_[v].from_timestep = from_t;
    graph_[v].to_timestep = to_t;
    return v;
}

void GraphicalModel::add_transition_vars(size_t from_t) {
	typedef boost::property_map<boostGraph, vertex_bundle_t>::type VarMap;
	typedef filtered_graph<boostGraph, keep_all, detections_at<VarMap> > filteredGraph;
        typedef graph_traits<filteredGraph>::vertex_iterator v_it;
	typedef pair<v_it, v_it> vp;

	    filteredGraph curr_detection_vars(graph_, 
					      keep_all(), 
					      detections_at<VarMap>(get(vertex_bundle, graph_), from_t));
	    vp curr_vp = vertices(curr_detection_vars);

	    filteredGraph next_detection_vars(graph_, 
					      keep_all(), 
					      detections_at<VarMap>(get(vertex_bundle, graph_), from_t+1));
	    vp next_vp = vertices(next_detection_vars);

	    for(v_it curr_it = curr_vp.first; curr_it != curr_vp.second; ++curr_it) {
		    for(v_it next_it = next_vp.first; next_it != next_vp.second; ++next_it) {
			boostVertex v = add_transition_var(from_t, from_t+1);
			add_edge(*curr_it, v, graph_);
			add_edge(v, *next_it, graph_);
		    }
	    }   
}

void GraphicalModel::add_transition_energies(size_t from_t, 
	    const BinaryEnergy& move,
	    const BinaryEnergy& mismove,
	    const TertiaryEnergy& division,
	    const vector< vector<Traxel> >& traxels) {
	typedef boost::property_map<boostGraph, vertex_bundle_t>::type VarMap;
	typedef filtered_graph<boostGraph, keep_all, detections_at<VarMap> > filteredGraph;
        typedef graph_traits<filteredGraph>::vertex_iterator v_it;
	typedef pair<v_it, v_it> vp;
        typedef graph_traits<boostGraph>::adjacency_iterator a_it;
	typedef pair<a_it, a_it> ap;
	filteredGraph curr_detection_vars(graph_, 
					      keep_all(), 
					      detections_at<VarMap>(get(vertex_bundle, graph_), from_t));
	    vp curr_vp = vertices(curr_detection_vars);
	    for(v_it curr_it = curr_vp.first; curr_it != curr_vp.second; ++curr_it) {
	      // only one adjacent vertex?
	      if(out_degree(*curr_it, graph_) == 1) {
       	        ap curr_ap = adjacent_vertices(*curr_it, graph_);
		size_t vi[] = {graph_[*curr_ap.first].idx};
		GraphicalModel::ogmFactor f(space_, vi, vi+1);

		// assign energies
		//get traxel_idxs
		size_t ancestor_idx = graph_[*curr_it].traxel_idx;
		ap next_vp = adjacent_vertices(*curr_ap.first, graph_);
		size_t descendant_idx = graph_[*next_vp.first].traxel_idx;
		++next_vp.first;
		assert(next_vp.first == next_vp.second); // tranistion nodes have a single outgoing edge

		//call energy functor
		Traxels prev = traxel_map_from_traxel_sequence(traxels[from_t].begin(), traxels[from_t].end());
		Traxels curr = traxel_map_from_traxel_sequence(traxels[from_t+1].begin(), traxels[from_t+1].end());
		f(1) = move(traxels[from_t][ancestor_idx], traxels[from_t + 1][descendant_idx], prev, curr);
		f(0) = mismove(traxels[from_t][ancestor_idx], traxels[from_t + 1][descendant_idx], prev, curr);;
		//cout << "f(0) " << f(0) << endl;
		//cout << "f(1) " << f(1) << endl;
		model_.addFactor(f);
		// several outgoing vertices (or none)
	      } else {
	    	    ap curr_ap = adjacent_vertices(*curr_it, graph_);
		    for(a_it curr_a = curr_ap.first; curr_a!= curr_ap.second; ++curr_a) {
			a_it next_a = curr_a;
			++next_a;
			    for(; next_a!= curr_ap.second; ++next_a) {
				//get traxel_idxs
				size_t ancestor_idx = graph_[*curr_it].traxel_idx;
				ap next_vp1 = adjacent_vertices(*curr_a, graph_);
				size_t descendant1_idx = graph_[*next_vp1.first].traxel_idx;
				++next_vp1.first;
				assert(next_vp1.first == next_vp1.second); // tranistion nodes have a single outgoing edge
				ap next_vp2 = adjacent_vertices(*next_a, graph_);
				size_t descendant2_idx = graph_[*next_vp2.first].traxel_idx;
				++next_vp2.first;
				assert(next_vp2.first == next_vp2.second); // tranistion nodes have a single outgoing edge

				//call energy functor
				Traxels prev = traxel_map_from_traxel_sequence(traxels[from_t].begin(), traxels[from_t].end());
				Traxels curr = traxel_map_from_traxel_sequence(traxels[from_t+1].begin(), traxels[from_t+1].end());
				size_t vi[] = {graph_[*curr_a].idx, graph_[*next_a].idx};
				GraphicalModel::ogmFactor f(space_, vi, vi+2);
				f(0,0) = mismove(traxels[from_t][ancestor_idx], traxels[from_t + 1][descendant1_idx], prev, curr) +
				         mismove(traxels[from_t][ancestor_idx], traxels[from_t + 1][descendant2_idx], prev, curr);
				f(0,1) = move(traxels[from_t][ancestor_idx], traxels[from_t + 1][descendant2_idx], prev, curr);
				f(1,0) = move(traxels[from_t][ancestor_idx], traxels[from_t + 1][descendant1_idx], prev, curr);
				f(1,1) = division(traxels[from_t][ancestor_idx],
						  traxels[from_t + 1][descendant1_idx],
						  traxels[from_t + 1][descendant2_idx],
						  prev,curr);
				//cout << "f(0,0) " << f(0,0) << endl;
				//cout << "f(0,1) " << f(0,1) << endl;
				//cout << "f(1,0) " << f(1,0) << endl;
				//cout << "f(1,1) " << f(1,1) << endl;
				model_.addFactor(f);
			    }
		    }
	      }   }
}



} /* namespace Tracking */
