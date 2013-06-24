#include <cassert>
#include <memory>
#include <set>
#include <string>
#include <iostream>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include "pgmlink/feature.h"
#include "pgmlink/pgm.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_pgm.h"
#include "pgmlink/tracking.h"

using namespace std;
using boost::shared_ptr;
using boost::shared_array;

namespace pgmlink {
////
//// class ChaingraphTracking
////

void ChaingraphTracking::set_with_divisions(bool state) {
	with_divisions_ = state;
}

void ChaingraphTracking::set_cplex_timeout(double seconds) {
	cplex_timeout_ = seconds;
}

vector<vector<Event> > ChaingraphTracking::operator()(TraxelStore& ts) {
  LOG(logINFO) << "Calling chaingraph tracking with the following parameters:\n"
	       << "\trandom forest filename: " << rf_fn_ << "\n"
	       << "\tappearance: " << app_ << "\n"
               << "\tdisappearance: " << dis_ << "\n"
    	       << "\tdetection: " << det_ << "\n"
    	       << "\tmisdetection: " << mis_  << "\n"
    	       << "\tcellness_by_random_forest: " << use_rf_  << "\n"
    	       << "\topportunity cost: " << opportunity_cost_ << "\n"
    	       << "\tforbidden cost: " << forbidden_cost_ << "\n"
    	       << "\twith constraints: " << with_constraints_ << "\n"
    	       << "\tfixed detections: " << fixed_detections_ << "\n"
    	       << "\tmean division distance: " << mean_div_dist_ << "\n"
    	       << "\tminimal division angle: " << min_angle_  << "\n"
    	       << "\tcplex ep gap: " << ep_gap_ << "\n"
    	       << "\tn neighbors: " <<  n_neighbors_ << "\n"
   	       << "\twith divisions: " << with_divisions_  << "\n"
   	       << "\tcplex timeout: " << cplex_timeout_ << "\n"
   	       << "\talternative builder: " << alternative_builder_;

	cout << "-> building feature functions " << endl;
	SquaredDistance move;
	BorderAwareConstant appearance(app_, earliest_timestep(ts), true, 0);
	BorderAwareConstant disappearance(dis_, latest_timestep(ts), false, 0);
	GeometryDivision2 division(mean_div_dist_, min_angle_);

	Traxels empty;
	// random forest?
	boost::function<double(const Traxel&)> detection, misdetection;
	if (use_rf_) {
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
	  detection = ConstantFeature(det_);
	  misdetection = ConstantFeature(mis_);
	}

	cout << "-> building hypotheses" << endl;
	SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(n_neighbors_, 50);
	SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
	boost::shared_ptr<HypothesesGraph> graph = boost::shared_ptr<HypothesesGraph>(hyp_builder.build());

	cout << "-> init MRF reasoner" << endl;
	std::auto_ptr<Chaingraph> mrf;

	if(alternative_builder_) {
	  pgm::chaingraph::TrainableModelBuilder b(appearance,
						 disappearance,
						 move,
						 opportunity_cost_,
						 forbidden_cost_);
	  
	  if (with_divisions_) {
		  b.with_divisions(division);
	  }

	  b.with_detection_vars(detection, misdetection);
	  mrf = std::auto_ptr<Chaingraph>(new Chaingraph(b, with_constraints_, ep_gap_, fixed_detections_, cplex_timeout_));
	} else {
	  pgm::chaingraph::ECCV12ModelBuilder b(appearance,
					      disappearance,
					      move,
					      opportunity_cost_,
					      forbidden_cost_);
	  
	  if (with_divisions_) {
		  b.with_divisions(division);
	  }

	  b.with_detection_vars(detection, misdetection);
	  mrf = std::auto_ptr<Chaingraph>(new Chaingraph(b, with_constraints_, ep_gap_, fixed_detections_, cplex_timeout_));
	}

	cout << "-> formulate MRF model" << endl;
	mrf->formulate(*graph);

	cout << "-> infer" << endl;
	mrf->infer();

	cout << "-> conclude" << endl;
	mrf->conclude(*graph);

	cout << "-> storing state of detection vars" << endl;
	last_detections_ = state_of_nodes(*graph);

	cout << "-> pruning inactive hypotheses" << endl;
	prune_inactive(*graph);

	cout << "-> constructing events" << endl;

	return *events(*graph);
}

vector<map<unsigned int, bool> > ChaingraphTracking::detections() {
	vector<map<unsigned int, bool> > res;
	if (last_detections_) {
		return *last_detections_;
	} else {
		throw std::runtime_error(
				"ChaingraphTracking::detections(): previous tracking result required");
	}
}



} // namespace tracking
