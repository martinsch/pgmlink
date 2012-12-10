#include <cassert>
#include <memory>
#include <set>
#include <string>
#include <iostream>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include "pgmlink/feature.h"
#include "pgmlink/graphical_model.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_opengm.h"
#include "pgmlink/tracking.h"
#include "pgmlink/reasoner_nearestneighbor.h"
#include "pgmlink/reasoner_nntracklets.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/track.h"

using namespace std;
using boost::shared_ptr;
using boost::shared_array;

namespace pgmlink {
////
//// class ChaingraphTracking
////
vector<vector<Event> > ChaingraphTracking::operator()(TraxelStore& ts) {
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
	SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(6, 50);
	SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
	shared_ptr<HypothesesGraph> graph = shared_ptr<HypothesesGraph>(hyp_builder.build());

	cout << "-> init MRF reasoner" << endl;
	std::auto_ptr<Chaingraph> mrf;

	if(alternative_builder_) {
	  pgm::TrainableChaingraphModelBuilder b(appearance,
						 disappearance,
						 move,
						 opportunity_cost_,
						 forbidden_cost_);
	  
	  b.with_divisions(division)
	   .with_detection_vars(detection, misdetection);
	  mrf = std::auto_ptr<Chaingraph>(new Chaingraph(b, with_constraints_, ep_gap_, fixed_detections_));
	} else {
	  pgm::ChaingraphModelBuilderECCV12 b(appearance,
					      disappearance,
					      move,
					      opportunity_cost_,
					      forbidden_cost_);
	  
	  b.with_divisions(division)
	   .with_detection_vars(detection, misdetection);
	  mrf = std::auto_ptr<Chaingraph>(new Chaingraph(b, with_constraints_, ep_gap_, fixed_detections_));
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


////
//// class NNTracking
////
vector<vector<Event> > NNTracking::operator()(TraxelStore& ts) {
	cout << "-> building hypotheses" << endl;
	SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(2, max(divDist_,movDist_));
	SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
	HypothesesGraph* graph = hyp_builder.build();
	HypothesesGraph& g = *graph;

	LOG(logDEBUG1) << "NNTracking: adding offered property to nodes";
	// adding 'offered' property and set it true for each node
	g.add(node_offered());
	property_map<node_offered, HypothesesGraph::base_graph>::type& offered_nodes = g.get(node_offered());
	// adding 'split_into' property
	g.add(split_from());
	property_map<split_from, HypothesesGraph::base_graph>::type& split_from_map = g.get(split_from());
	for(HypothesesGraph::NodeIt n(g); n!=lemon::INVALID; ++n) {
		offered_nodes.set(n,true);
		split_from_map.set(n,-1);
	}
	LOG(logDEBUG1) << "NNTracking: adding distance property to edges";
	g.add(arc_distance());
	property_map<arc_distance, HypothesesGraph::base_graph>::type& arc_distances = g.get(arc_distance());
	property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
	for(HypothesesGraph::ArcIt a(g); a!=lemon::INVALID; ++a) {
		HypothesesGraph::Node from = g.source(a);
		HypothesesGraph::Node to = g.target(a);
		Traxel from_tr = traxel_map[from];
		Traxel to_tr = traxel_map[to];

		double dist = 0;
		// if we want to add another dimension to the norm, we remove the sqrt, add the dimensions and sqrt again

		for(std::vector<std::string>::const_iterator it = distanceFeatures_.begin(); it!=distanceFeatures_.end(); ++it) {
			if (*it == "com") {
				// com is already considered in Traxel::distance_to()
				// Traxel::distance_to computes: sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
				double d = from_tr.distance_to(to_tr);
				LOG(logDEBUG2) << "NNTracking:: com distance from " << from_tr.Id << " to " << to_tr.Id << " = " << d;
				dist += (d*d);
			} else {
				std::vector<float> from_feat = from_tr.features.find(*it)->second;
				std::vector<float> to_feat = to_tr.features.find(*it)->second;
				for (size_t i = 0; i<from_feat.size(); ++i) {
					// TODO: do we have to consider x/y/z scale for some features?
					double d = (from_feat[i] - to_feat[i]);
					dist += (d*d);
				}
			}
		}
		dist = sqrt(dist);
		LOG(logDEBUG2) << "NNTracking:: combined distance from " << from_tr.Id << " to " << to_tr.Id << " = " << dist;

		arc_distances.set(a, dist);
	}


	cout << "-> init NN reasoner" << endl;
	SingleTimestepTraxelNN nn_reasoner(divDist_,movDist_,divisionThreshold_,splitterHandling_, mergerHandling_, maxTraxelIdAt_);

	cout << "-> formulate NN model" << endl;
	nn_reasoner.formulate(*graph);

	cout << "-> infer" << endl;
	nn_reasoner.infer();

	cout << "-> conclude" << endl;
	nn_reasoner.conclude(*graph);

	cout << "-> storing state of detection vars" << endl;
	last_detections_ = state_of_nodes(*graph);

	cout << "-> pruning inactive hypotheses" << endl;
	prune_inactive(*graph);

	cout << "-> constructing events" << endl;

	return *events(*graph);
}

vector<map<unsigned int, bool> > NNTracking::detections() {
	vector<map<unsigned int, bool> > res;
	if (last_detections_) {
		return *last_detections_;
	} else {
		throw std::runtime_error(
				"NNTracking::detections(): previous tracking result required");
	}
}



////
//// class NNTrackletsTracking
////
vector<vector<Event> > NNTrackletsTracking::operator()(TraxelStore& ts) {
	cout << "-> building hypotheses" << endl;
	SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(1, // max_nearest_neighbors
			maxDist_,
			true, // forward_backward
			true, // consider_divisions
			divisionThreshold_
			);
	SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
	HypothesesGraph* graph = hyp_builder.build();
	HypothesesGraph& g = *graph;

	LOG(logDEBUG1) << "NNTrackletsTracking: adding offered property to nodes";
	// adding 'offered' property and set it true for each node
	g.add(node_offered());
	property_map<node_offered, HypothesesGraph::base_graph>::type& offered_nodes = g.get(node_offered());
	// adding 'split_into' property
	g.add(split_from());
	property_map<split_from, HypothesesGraph::base_graph>::type& split_from_map = g.get(split_from());
	for(HypothesesGraph::NodeIt n(g); n!=lemon::INVALID; ++n) {
		offered_nodes.set(n,true);
		split_from_map.set(n,-1);
	}
	LOG(logDEBUG1) << "NNTrackletsTracking: adding distance property to edges";
	g.add(arc_distance());
	property_map<arc_distance, HypothesesGraph::base_graph>::type& arc_distances = g.get(arc_distance());
	property_map<arc_distance, HypothesesGraph::base_graph>::type& arc_vol_ratios = g.get(arc_vol_ratio());
	property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
	for(HypothesesGraph::ArcIt a(g); a!=lemon::INVALID; ++a) {
		HypothesesGraph::Node from = g.source(a);
		HypothesesGraph::Node to = g.target(a);
		Traxel from_tr = traxel_map[from];
		Traxel to_tr = traxel_map[to];

//		double dist = 0;
//		// if we want to add another dimension to the norm, we remove the sqrt, add the dimensions and sqrt again
//
//		for(std::vector<std::string>::const_iterator it = distanceFeatures_.begin(); it!=distanceFeatures_.end(); ++it) {
//			if (*it == "com") {
//				// com is already considered in Traxel::distance_to()
//				// Traxel::distance_to computes: sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
//				double d = from_tr.distance_to(to_tr);
//				LOG(logDEBUG3) << "NNTrackletsTracking: com distance from " << from_tr.Id << " to " << to_tr.Id << " = " << d;
//				dist += (d*d);
//			} else {
//				std::vector<float> from_feat = from_tr.features.find(*it)->second;
//				std::vector<float> to_feat = to_tr.features.find(*it)->second;
//				for (size_t i = 0; i<from_feat.size(); ++i) {
//					// TODO: do we have to consider x/y/z scale for some features?
//					double d = (from_feat[i] - to_feat[i]);
//					dist += (d*d);
//				}
//			}
//		}
//		dist = sqrt(dist);
		arc_distances.set(a, from_tr.distance_to(to_tr));
		LOG(logDEBUG2) << "NNTrackletsTracking: combined distance from " << from_tr.Id << " to " << to_tr.Id << " = " << from_tr.distance_to(to_tr);

		std::vector<float> from_vol = from_tr.features.find("count")->second;
		std::vector<float> to_vol = to_tr.features.find("count")->second;
		double ratio = from_vol[0] / double(to_vol[0]);
		arc_vol_ratios.set(a, ratio);
		LOG(logDEBUG2) << "NNTrackletsTracking: volume ratio for arc from " << from_tr.Id << " to " << to_tr.Id << " = " << ratio;

	}


	cout << "-> init NN reasoner" << endl;
	SingleTimestepTraxelNNTracklets nn_reasoner(maxDist_,divisionThreshold_,splitterHandling_, mergerHandling_, maxTraxelIdAt_);

	cout << "-> formulate NN model" << endl;
	nn_reasoner.formulate(*graph);

	cout << "-> infer" << endl;
	nn_reasoner.infer();

	cout << "-> conclude" << endl;
	nn_reasoner.conclude(*graph);

	cout << "-> storing state of detection vars" << endl;
	last_detections_ = state_of_nodes(*graph);

	cout << "-> pruning inactive hypotheses" << endl;
	prune_inactive(*graph);

	cout << "-> constructing events" << endl;

	return *events(*graph);
}

vector<map<unsigned int, bool> > NNTrackletsTracking::detections() {
	vector<map<unsigned int, bool> > res;
	if (last_detections_) {
		return *last_detections_;
	} else {
		throw std::runtime_error(
				"NNTracking::detections(): previous tracking result required");
	}
}


namespace {
std::vector<double> computeDetProb(double vol, double avg_vol, vector<double> s2) {
	std::vector<double> result;

	double sum = 0;
	size_t k = 0;
	for (vector<double>::const_iterator it = s2.begin(); it != s2.end(); ++it) {
		double val = vol - k*avg_vol;
		val = exp(-(val*val)/(2*(*it)));
		if (k==0) {
			val *= 1.5;
		}
		result.push_back(val);
		sum += val;
		++k;
	}

	// normalize
	for(std::vector<double>::iterator it = result.begin(); it!=result.end(); ++it) {
		(*it) /= sum;
	}

	return result;
}
}


////
//// class ConsTracking
////
vector<vector<Event> > ConsTracking::operator()(TraxelStore& ts) {
	cout << "-> building energy functions " << endl;

	Traxels empty;
	boost::function<double(const Traxel&, const size_t)> detection, division;
	boost::function<double(const double)> transition;

	if (use_size_dependent_detection_) {
//		LOG(logINFO) << "Loading Random Forest";
//		vigra::RandomForest<RF::RF_LABEL_TYPE> rf = RF::getRandomForest(detection_rf_fn_);
//		std::vector<std::string> rf_features;
//		rf_features.push_back("volume");
//		rf_features.push_back("bbox");
//		rf_features.push_back("position");
//		rf_features.push_back("com");
//		rf_features.push_back("pc");
//		rf_features.push_back("intensity");
//		rf_features.push_back("intminmax");
//		rf_features.push_back("pair");
//		rf_features.push_back("sgf");
//		rf_features.push_back("lcom");
//		rf_features.push_back("lpc");
//		rf_features.push_back("lintensity");
//		rf_features.push_back("lintminmax");
//		rf_features.push_back("lpair");
//		rf_features.push_back("lsgf");

//		LOG(logINFO) << "Predicting cellness";
//		RF::predict_traxels(ts, rf, rf_features, 1, "cellness");

		double s2 = (avg_obj_size_*avg_obj_size_)/4.0;
		LOG(logDEBUG) << "sigmas are all set to " << s2;
		vector<double> sigma2(max_number_objects_+1,s2);
		for(size_t i = 0; i < sigma2.size(); ++i) {
			cout << sigma2[i] << endl;
		}

		for(TraxelStore::iterator tr = ts.begin(); tr != ts.end(); ++tr) {
			Traxel trax = *tr;
			FeatureMap::const_iterator it = trax.features.find("count");
			if(it == trax.features.end()) {
				throw runtime_error("get_detection_prob(): cellness feature not in traxel");
			}
			double vol = it->second[0];
			vector<double> detProb;
			detProb = computeDetProb(vol,avg_obj_size_,sigma2);
			feature_array detProbFeat(feature_array::difference_type(max_number_objects_+1));
			for(int i = 0; i<=max_number_objects_; ++i) {
				LOG(logDEBUG2) << "detection probability for " << trax.Id << "[" << i << "] = " << detProb[i];
				detProbFeat[i] = detProb[i];
			}
			trax.features["detProb"] = detProbFeat;
			ts.replace(tr, trax);
		}
		detection = NegLnDetection(1); // weight 1
	} else {
		// assume a quasi geometric distribution
		vector<double> prob_vector;
		double p = 0.7; // e.g. for max_number_objects_=3, p=0.7: P(X=(0,1,2,3)) = (0.027, 0.7, 0.21, 0.063)
		double sum = 0;
		for(double state = 0; state < max_number_objects_; ++state) {
			double prob = p*pow(1-p,state);
			prob_vector.push_back(prob);
			sum += prob;
		}
		prob_vector.insert(prob_vector.begin(), 1-sum);

		detection = bind<double>(NegLnConstant(1,prob_vector), _2);
	}

	division = NegLnDivision(1); // weight 1
	transition = NegLnTransition(1); // weight 1

	cout << "-> building hypotheses" << endl;
	SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(1, // max_nearest_neighbors
				max_dist_,
				true, // forward_backward
				true, // consider_divisions
				division_threshold_
				);
	SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
	HypothesesGraph* graph = hyp_builder.build();


	LOG(logDEBUG1) << "ConsTracking(): adding distance property to edges";
	HypothesesGraph& g = *graph;
	g.add(arc_distance());
	property_map<arc_distance, HypothesesGraph::base_graph>::type& arc_distances = g.get(arc_distance());
	property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
	for(HypothesesGraph::ArcIt a(g); a!=lemon::INVALID; ++a) {
		HypothesesGraph::Node from = g.source(a);
		HypothesesGraph::Node to = g.target(a);
		Traxel from_tr = traxel_map[from];
		Traxel to_tr = traxel_map[to];

		double dist = from_tr.distance_to(to_tr);
		arc_distances.set(a, dist);
	}

	cout << "-> init ConservationTracking reasoner" << endl;
	SingleTimestepTraxelConservation pgm(
			max_number_objects_,
			detection,
			division,
			transition,
			forbidden_cost_,
			with_constraints_,
			fixed_detections_,
			ep_gap_
			);

	cout << "-> formulate ConservationTracking model" << endl;
	pgm.formulate(*graph);

	cout << "-> infer" << endl;
	pgm.infer();

	cout << "-> conclude" << endl;
	pgm.conclude(*graph);

	cout << "-> storing state of detection vars" << endl;
	last_detections_ = state_of_nodes(*graph);

	cout << "-> pruning inactive hypotheses" << endl;
	prune_inactive(*graph);

	cout << "-> constructing events" << endl;

	return *events(*graph);
}

vector<map<unsigned int, bool> > ConsTracking::detections() {
	vector<map<unsigned int, bool> > res;
	if (last_detections_) {
		return *last_detections_;
	} else {
		throw std::runtime_error(
				"MrfTracking::detections(): previous tracking result required");
	}
}


} // namespace tracking
