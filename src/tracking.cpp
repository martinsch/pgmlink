#include <cassert>
#include <set>
#include <string>
#include <iostream>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include "pgmlink/energy.h"
#include "pgmlink/graphical_model.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_opengm.h"
#include "pgmlink/tracking.h"

using namespace std;
using boost::shared_ptr;
using boost::shared_array;

namespace pgmlink {
////
//// class ChaingraphTracking
////
vector<vector<Event> > ChaingraphTracking::operator()(TraxelStore& ts) {
	cout << "-> building energy functions " << endl;
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
		detection = bind<double>(ConstantEnergy(det_), _1, empty, empty);
		misdetection = bind<double>(ConstantEnergy(mis_), _1, empty, empty);
		;
	}

	cout << "-> building hypotheses" << endl;
	SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(6, 50);
	SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
	HypothesesGraph* graph = hyp_builder.build();

	cout << "-> init MRF reasoner" << endl;
	Chaingraph mrf(detection, misdetection, appearance,
			disappearance, bind<double>(move, _1, _2, empty, empty), division,
			opportunity_cost_, forbidden_cost_, with_constraints_,
			fixed_detections_, ep_gap_);

	cout << "-> formulate MRF model" << endl;
	mrf.formulate(*graph);

	cout << "-> infer" << endl;
	mrf.infer();

	cout << "-> conclude" << endl;
	mrf.conclude(*graph);

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
