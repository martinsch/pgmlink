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
#include "pgmlink/reasoner_nearestneighbor.h"
#include "pgmlink/reasoner_nntracklets.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/merger_resolving.h"

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
	} else if (ts.begin()->features.find("detProb") != ts.begin()->features.end()) {
          for (TraxelStore::iterator it = ts.begin(); it != ts.end(); ++it) {
            Traxel trax = *it;
            trax.features["cellness"] = trax.features["detProb"];
            assert(trax.features["detProb"].size() == 2);
            ts.replace(it, trax);
          }
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
	NnTracking nn_reasoner(divDist_,movDist_,divisionThreshold_,splitterHandling_, mergerHandling_, maxTraxelIdAt_);

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
	NnTrackletTracking nn_reasoner(maxDist_,divisionThreshold_,splitterHandling_, mergerHandling_, maxTraxelIdAt_);

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
std::vector<double> computeDetProb(double vol, vector<double> means, vector<double> s2) {
	std::vector<double> result;

	double sum = 0;
	for (size_t k = 0; k < means.size(); ++k) {
		double val = vol - means[k];
		val = exp(-(val*val)/s2[k]);
		result.push_back(val);
		sum += val;
	}

	// normalize
	for(std::vector<double>::iterator it = result.begin(); it!=result.end(); ++it) {
		(*it) /= sum;
	}

	return result;
}

double dot(double x1,double y1,double z1, double x2,double y2,double z2) {
      return x1*x2 + y1*y2 + z1*z2;
}

double norm(double x,double y,double z) {
      return sqrt(dot(x,y,z, x,y,z));
}

double getCorrectedDistance(Traxel from, Traxel to) {
	FeatureMap::const_iterator it = from.features.find("com_corrected");
	if (it == from.features.end()) {
		throw runtime_error("getCorrectedDistance(): com_corrected feature not found in traxel");
	}
	return norm(it->second[0]-to.X(),it->second[1]-to.Y(),it->second[2]-to.Z());
}
}

////
//// helper function needed for boost::algorithm::all_of
//// 
template <typename T>
bool has_data(const std::vector<T>& vector) {
  return vector.size() > 0;
}


////
//// helper function equivalent to all_of (c++11)
////
template<class InputIterator, class UnaryPredicate>
bool all_true (InputIterator first, InputIterator last, UnaryPredicate pred) {
  while (first!=last) {
    if (!pred(*first)) return false;
    ++first;
  }
  return true;
}


////
//// class ConsTracking
////
vector<vector<Event> > ConsTracking::operator()(TraxelStore& ts) {
	cout << "-> building energy functions " << endl;

	double detection_weight = 10;
	Traxels empty;
	boost::function<double(const Traxel&, const size_t)> detection, division;
	boost::function<double(const double)> transition;

	bool use_classifier_prior = false;
	Traxel trax = *(ts.begin());
	FeatureMap::const_iterator it = trax.features.find("detProb");
	if(it != trax.features.end()) {
		use_classifier_prior = true;
	}
	if (use_classifier_prior) {
		LOG(logINFO) << "Using classifier prior";
		detection = NegLnDetection(detection_weight);
	} else if (use_size_dependent_detection_) {
		LOG(logINFO) << "Using size dependent prior";
		vector<double> means;
		if (means_.size() == 0 ) {
			for(int i = 0; i<max_number_objects_+1; ++i) {
				means.push_back(i*avg_obj_size_);
				LOG(logINFO) << "mean[" << i << "] = " << means[i];
			}
		} else {
			assert(sigmas_.size() != 0);
			for(int i = 0; i<max_number_objects_+1; ++i) {
				means.push_back(means_[i]);
				LOG(logINFO) << "mean[" << i << "] = " << means[i];
			}
		}

		vector<double> sigma2;
		if (sigmas_.size() == 0) {
			double s2 = (avg_obj_size_*avg_obj_size_)/4.0;
			if (s2 < 0.0001) {
				s2 = 0.0001;
			}
			for(int i = 0; i<max_number_objects_+1; ++i) {
				sigma2.push_back(s2);
				LOG(logINFO) << "sigma2[" << i << "] = "  << sigma2[i];
			}
		} else {
			for (int i = 0; i<max_number_objects_+1; ++i) {
				sigma2.push_back(sigmas_[i]);
				LOG(logINFO) << "sigma2[" << i << "] = "  << sigma2[i];
			}
		}

		for(TraxelStore::iterator tr = ts.begin(); tr != ts.end(); ++tr) {
			Traxel trax = *tr;
			FeatureMap::const_iterator it = trax.features.find("count");
			if(it == trax.features.end()) {
				throw runtime_error("get_detection_prob(): cellness feature not in traxel");
			}
			double vol = it->second[0];
			vector<double> detProb;
			detProb = computeDetProb(vol,means,sigma2);
			feature_array detProbFeat(feature_array::difference_type(max_number_objects_+1));
			for(int i = 0; i<=max_number_objects_; ++i) {
				double d = detProb[i];
				if (d < 0.01) {
					d = 0.01;
				} else if (d > 0.99) {
					d = 0.99;
				}
				LOG(logDEBUG2) << "detection probability for " << trax.Id << "[" << i << "] = " << d;
				detProbFeat[i] = d;
			}
			trax.features["detProb"] = detProbFeat;
			ts.replace(tr, trax);
		}
		detection = NegLnDetection(detection_weight); // weight 1
	} else {
		LOG(logINFO) << "Using hard prior";
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

		detection = bind<double>(NegLnConstant(detection_weight,prob_vector), _2);
	}

	LOG(logDEBUG1) << "division_weight_ = " << division_weight_;
	LOG(logDEBUG1) << "transition_weight_ = " << transition_weight_;
	division = NegLnDivision(division_weight_);
	transition = NegLnTransition(transition_weight_);

	cout << "-> building hypotheses" << endl;
	SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(1, // max_nearest_neighbors
				max_dist_,
				true, // forward_backward
				with_divisions_, // consider_divisions
				division_threshold_
				);
	SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
	HypothesesGraph* graph = hyp_builder.build();


	LOG(logDEBUG1) << "ConsTracking(): adding distance property to edges";
	HypothesesGraph& g = *graph;
	g.add(arc_distance()).add(tracklet_intern_dist()).add(node_tracklet()).add(tracklet_intern_arc_ids()).add(traxel_arc_id());
	property_map<arc_distance, HypothesesGraph::base_graph>::type& arc_distances = g.get(arc_distance());
	property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
	bool with_optical_correction = false;
	Traxel some_traxel = (*traxel_map.beginValue());
	if (some_traxel.features.find("com_corrected") != some_traxel.features.end()) {
		LOG(logINFO) << "optical correction enabled";
		with_optical_correction = true;
	}

	for(HypothesesGraph::ArcIt a(g); a!=lemon::INVALID; ++a) {
		HypothesesGraph::Node from = g.source(a);
		HypothesesGraph::Node to = g.target(a);
		Traxel from_tr = traxel_map[from];
		Traxel to_tr = traxel_map[to];

		if (with_optical_correction) {
			arc_distances.set(a, getCorrectedDistance(from_tr,to_tr));
		} else {
			arc_distances.set(a, from_tr.distance_to(to_tr));
		}
	}
	//border_width_ is given in normalized scale, 1 corresponds to a maximal distance of dim_range/2
	boost::function<double(const Traxel&)> appearance_cost_fn, disappearance_cost_fn;
	LOG(logINFO) << "  using border-aware appearance and disappearance costs, with margin: " << border_width_;
	appearance_cost_fn = SpatialBorderAwareWeight(appearance_cost_, border_width_, fov_);
	disappearance_cost_fn = SpatialBorderAwareWeight(disappearance_cost_, border_width_, fov_);

	cout << "-> init ConservationTracking reasoner" << endl;
	ConservationTracking pgm(
			max_number_objects_,
			detection,
			division,
			transition,
			forbidden_cost_,
			ep_gap_,
			with_tracklets_,
			with_divisions_,
			disappearance_cost_fn,
			appearance_cost_fn,
			transition_parameter_
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

        cout << "-> constructing unresolved events" << endl;
        boost::shared_ptr<std::vector< std::vector<Event> > > ev = events(*graph);
        

        if (with_merger_resolution_ && all_true(ev->begin(), ev->end(), has_data<Event>)) {
          cout << "-> resolving mergers" << endl;
          MergerResolver m(graph);
          // FeatureExtractorMCOMsFromKMeans extractor;
          FeatureExtractorMCOMsFromGMM extractor(number_of_dimensions_);
          DistanceFromCOMs distance;
          FeatureHandlerFromTraxels handler(extractor, distance);
          m.resolve_mergers(handler);
        
          HypothesesGraph g_res;
          resolve_graph(*graph, g_res, transition, ep_gap_, with_tracklets_);
          prune_inactive(*graph);

          cout << "-> constructing resolved events" << endl;
          boost::shared_ptr<std::vector< std::vector<Event> > > multi_frame_moves = multi_frame_move_events(*graph);

          cout << "-> merging unresolved and resolved events" << endl;
          return *merge_event_vectors(*ev, *multi_frame_moves);
        }

        else {
          return *ev;
        }

        

	

	
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
