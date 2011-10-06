#include <cassert>
#include <set>
#include <string>
#include <iostream>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include "energy.h"
#include "ilp_construction.h"
#include "graphical_model.h"
#include "hypotheses.h"
#include "log.h"
#include "reasoning/mrf_reasoner.h"
#ifdef USE_CPLEX
#include "cplex_solver.h"
#else
#include "lp_solve_solver.h"
#endif
#include "track.h"

using namespace std;
using boost::shared_ptr;
using boost::shared_array;

namespace Tracking {
////
//// class MrfTracking
////
vector<vector<Event> > MrfTracking::operator()(TraxelStore& ts) {
    cout << "-> building energy functions " << endl;
    SquaredDistance move;
    ConstantEnergy appearance(app_);
    ConstantEnergy disappearance(dis_);
    GeometryDivision2 division(mean_div_dist_, min_angle_);

    Traxels empty;
    // random forest?
    boost::function<double (const Traxel&)> detection, misdetection;
    if(use_rf_) {
      LOG(logINFO) << "Loading Random Forest";
      vigra::RandomForest<RF::RF_LABEL_TYPE> rf = RF::getRandomForest(rf_fn_);
      std::vector<std::string> rf_features;
      rf_features.push_back("volume");
      rf_features.push_back("bbox");
      rf_features.push_back("position");
      rf_features.push_back("com");
      rf_features.push_back("pc");
      rf_features.push_back("intensity");
      rf_features.push_back("intminmax");
      rf_features.push_back("pair");
      rf_features.push_back("sgf");
      rf_features.push_back("lcom");
      rf_features.push_back("lpc");
      rf_features.push_back("lintensity");
      rf_features.push_back("lintminmax");
      rf_features.push_back("lpair");
      rf_features.push_back("lsgf");

      LOG(logINFO) << "Predicting cellness";
      RF::predict_traxels(ts, rf, rf_features, 1, "cellness");

      detection = NegLnCellness(det_);
      misdetection = NegLnOneMinusCellness(mis_);
    } else {
      detection = bind<double>(ConstantEnergy(det_), _1, empty, empty);
      misdetection = bind<double>(ConstantEnergy(mis_), _1, empty, empty);;
    }

    cout << "-> building hypotheses" << endl;
    SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts);
    HypothesesGraph* graph = hyp_builder.build();

    cout << "-> init MRF reasoner" << endl;
    SingleTimestepTraxelMrf mrf(detection,
			        misdetection,
			        bind<double>(appearance, _1, empty, empty),
			        bind<double>(disappearance, _1, empty, empty),
				bind<double>(move, _1, _2, empty, empty),
				division,
				opportunity_cost_
			       );

    cout << "-> formulate MRF model" << endl;        
    mrf.formulate(*graph);
    
    cout << "-> infer" << endl;
    mrf.infer();
    
    cout << "-> conclude" << endl;
    mrf.conclude(*graph);

    cout << "-> pruning inactive hypotheses" << endl;
    prune_inactive(*graph);

    cout << "-> constructing events" << endl;
    
    return *events(*graph);
  }



////
//// class Track
////
vector<Event> Track::operator()(const Traxels& prev, const Traxels& curr) {
  ////
  //// Formulate ilp
  ////
  pair<shared_ptr<AdaptiveEnergiesIlp>, vector<Event> > ret = formulation_.formulate_ilp(prev, curr);
  shared_ptr<AdaptiveEnergiesIlp> ilp = ret.first;
  vector<Event> ilp_events = ret.second;

  ////
  //// Solve the ilp
  ////
  // output variables
  double finalCost = 0;
  shared_array<double> finalVars;
  int solverOutput = 0;
  string outputString = "solver not called";

  // traxels too sparsely distributed? 
  if(ilp->nVars > 0) {
#ifdef USE_CPLEX
    CplexSolver solver;
#else
    LpSolveSolver solver;
#endif
        
    int err;
    err = solver.solve(ilp->nVars, ilp->nConstr, ilp->nNonZero,
		       ilp->costs,
		       ilp->rhs,
		       ilp->matbeg,
		       ilp->matcnt,
		       ilp->matind,
		       ilp->matval,
		       string("Adaptive Energies Approach"),
		       finalCost,
		       finalVars,
		       solverOutput,
		       outputString
		       );
    if(err != 0) {
      throw "solver returned error code";
    }
    cout << "Final Cost: " << finalCost << endl;
    cout << "Solver Return Value: " << solverOutput << endl;
    cout << "Solver Output: " << outputString << endl;
  }
  ////
  //// Construct Events from ilp output
  ////
  shared_ptr<vector<Event> > events = events_from(prev, curr, finalVars, ilp->nVars, ilp_events);

  return *events;
}



  shared_ptr<vector<Event> > Track::events_from(const Traxels& prev, const Traxels& curr,
				       shared_array<double> ilp_finalVars,
				       int ilp_nVars,
				       const vector<Event>& ilp_events) const{
  // collect traxel ids
  pair<set<unsigned int>::iterator, bool > id_insert_ret; 
  set<unsigned int> prev_ids;
  for(Traxels::const_iterator it = prev.begin(); it != prev.end(); ++it) {
    id_insert_ret = prev_ids.insert(it->first); 
    if(id_insert_ret.second == false) {
      throw "duplicate traxel ids detected";
    }
  }
  set<unsigned int> curr_ids;
  for(Traxels::const_iterator it = curr.begin(); it != curr.end(); ++it) {
    id_insert_ret = curr_ids.insert(it->first); 
    if(id_insert_ret.second == false) {
      throw "duplicate traxel ids detected";
    }
  }
  

  // remove ids, that belong to Move or Division events and add the
  // Events, that actually happened to the output
  //
  // reset the energy:
  // energy has a different interpretation in the ilp context
  // ilp: energy are the costs corresponding to a binary variable
  // else: energy as a measure of likelihood for the occurence of an event

  Event e; // working var
  shared_ptr<vector<Event> > events(new vector<Event>);
  assert(static_cast<vector<Event>::size_type>(ilp_nVars) == ilp_events.size());
  for(int i = 0; i < ilp_nVars; ++i) {
    // the Move or Divison event actually happend
    if(ilp_finalVars[i] == 1) {
      if(ilp_events[i].type == Event::Move) {
	e = ilp_events[i];
	assert(e.traxel_ids.size() == 2);
	prev_ids.erase(e.traxel_ids[0]);
	curr_ids.erase(e.traxel_ids[1]);

	shared_ptr<const BinaryEnergy> energy_fct = formulation_.Move_energy();
	e.energy = energy_fct->operator()(prev.find(e.traxel_ids[0])->second,
					  curr.find(e.traxel_ids[1])->second, 
					  prev, curr);
	events->push_back(e);

      } else if(ilp_events[i].type == Event::Division) {
	e = ilp_events[i];
	assert(e.traxel_ids.size() == 3);
	prev_ids.erase(e.traxel_ids[0]);
	curr_ids.erase(e.traxel_ids[1]);
	curr_ids.erase(e.traxel_ids[2]);

	shared_ptr<const TertiaryEnergy> energy_fct = formulation_.Division_energy();
	e.energy = energy_fct->operator()(prev.find(e.traxel_ids[0])->second,
					  curr.find(e.traxel_ids[1])->second, 
					  curr.find(e.traxel_ids[2])->second, 
					  prev, curr);
	events->push_back(e);
      } else {
	throw "Division or Move expected";
      }
    }
  }
  
  // every id left has to be a (dis-)appearance
  // construct disappearances
  for(set<unsigned int>::const_iterator it = prev_ids.begin(); it!=prev_ids.end(); ++it) {
    e.type = Event::Disappearance;
    e.traxel_ids.clear();
    e.traxel_ids.push_back(*it);
    shared_ptr<const UnaryEnergy> energy_fct = formulation_.Disappearance_energy();
    e.energy = energy_fct->operator()(prev.find(*it)->second, prev, curr);
    events->push_back(e);
  }
  // construct appearances
  for(set<unsigned int>::const_iterator it = curr_ids.begin(); it!=curr_ids.end(); ++it) {
    e.type = Event::Appearance;
    e.traxel_ids.clear();
    e.traxel_ids.push_back(*it);
    shared_ptr<const UnaryEnergy> energy_fct = formulation_.Appearance_energy();
    e.energy = energy_fct->operator()(curr.find(*it)->second, prev, curr);
    events->push_back(e);
  }
  return events;
}



  Track& Track::Formulation(const AdaptiveEnergiesFormulation& f) {
    formulation_ = f;
    return *this;
  }

  const AdaptiveEnergiesFormulation& Track::Formulation() const {
    return formulation_;
  } 



////
//// class FixedCostTracking
////
FixedCostTracking::FixedCostTracking(double division, double move, double disappearance, double appearance,
    double distance_threshold) {
  shared_ptr<const ConstantEnergy> div(new ConstantEnergy(division));
  shared_ptr<const ConstantEnergy> mov(new ConstantEnergy(move));
  shared_ptr<const ConstantEnergy> app(new ConstantEnergy(appearance));
  shared_ptr<const ConstantEnergy> dis(new ConstantEnergy(disappearance));
  formulation_ = AdaptiveEnergiesFormulation(div, mov, dis, app, distance_threshold);
}



vector<Event> FixedCostTracking::operator()(const Traxels& prev, const Traxels& curr) {
    Track track(formulation_);
    return track(prev, curr);
}



////
//// class ShortestDistanceTracking
////
ShortestDistanceTracking::ShortestDistanceTracking(double division, double disappearance, double appearance,
						   double distance_threshold, unsigned int max_nearest_neighbors) {
  shared_ptr<const KasterDivision> div(new KasterDivision(division));
  shared_ptr<const SquaredDistance> mov(new SquaredDistance());
  shared_ptr<const ConstantEnergy> app(new ConstantEnergy(appearance));
  shared_ptr<const ConstantEnergy> dis(new ConstantEnergy(disappearance));
  formulation_ = AdaptiveEnergiesFormulation(div, mov, dis, app, distance_threshold, "com", max_nearest_neighbors);
}



vector<Event> ShortestDistanceTracking::operator()(const Traxels& prev, const Traxels& curr) {
    Track track(formulation_);
    return track(prev, curr);
}



////
//// class CellnessTracking
////
CellnessTracking::CellnessTracking(const std::string& random_forest_file,
                                   double w_div1, double w_div2,
                                   double w_move, double w_app, double w_disapp,
				   double distance_threshold)
                                       : random_forest_(RF::getRandomForest(random_forest_file))
{
    // set up the content and order of the random forest feature vector   
    rf_features_.push_back("volume");
    rf_features_.push_back("bbox");
    rf_features_.push_back("position");
    rf_features_.push_back("com");
    rf_features_.push_back("pc");
    rf_features_.push_back("intensity");
    rf_features_.push_back("intminmax");
    rf_features_.push_back("pair");
    rf_features_.push_back("sgf");
    rf_features_.push_back("lcom");
    rf_features_.push_back("lpc");
    rf_features_.push_back("lintensity");
    rf_features_.push_back("lintminmax");
    rf_features_.push_back("lpair");
    rf_features_.push_back("lsgf");

    // setup cellness energy functors
    shared_ptr<const CellnessDivision> div(new CellnessDivision(w_div1,w_div2));
    shared_ptr<const CellnessMove> mov(new CellnessMove(w_move));
    shared_ptr<const CellnessAppearance> app(new CellnessAppearance(w_app));
    shared_ptr<const CellnessDisappearance> dis(new CellnessDisappearance(w_disapp));
    formulation_ = AdaptiveEnergiesFormulation(div, mov, dis, app, distance_threshold);
}




vector<Event> CellnessTracking::operator()(Traxels prev, Traxels curr) {
    // predict cellness and add as feature to the traxels
    RF::predictTracklets(prev, random_forest_, rf_features_, 1, "cellness");
    RF::predictTracklets(curr, random_forest_, rf_features_, 1, "cellness");

    Track track(formulation_);
    return track(prev, curr);
}
} // namespace tracking
