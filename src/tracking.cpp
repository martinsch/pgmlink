#include <cassert>
#include <memory>
#include <set>
#include <string>
#include <iostream>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "pgmlink/feature.h"
#include "pgmlink/pgm.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_pgm.h"
#include "pgmlink/tracking.h"
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
  vector<vector<Event> > ConsTracking::operator()(TraxelStore& ts,
						  double forbidden_cost,
						  double ep_gap,
						  bool with_tracklets,
						  double division_weight,
						  double transition_weight,
						  double disappearance_cost,
						  double appearance_cost,
						  bool with_merger_resolution,
						  int n_dim,
						  double transition_parameter,
						  double border_width,
						  bool with_constraints,
						  double cplex_timeout,
						  TimestepIdCoordinateMapPtr coordinates) {
    build_hypo_graph(ts);
    return track(forbidden_cost,ep_gap,with_tracklets,division_weight,transition_weight,disappearance_cost,appearance_cost,with_merger_resolution,n_dim,transition_parameter,border_width,with_constraints,cplex_timeout,coordinates);
}

shared_ptr<HypothesesGraph> ConsTracking::build_hypo_graph(TraxelStore& ts) {

  
  LOG(logDEBUG3) << "enering build_hypo_graph"<< endl;;
  
	traxel_store_ = &ts;

	use_classifier_prior_ = false;
	Traxel trax = *(traxel_store_->begin());
	FeatureMap::const_iterator it = trax.features.find("detProb");
	if(it != trax.features.end()) {
	        use_classifier_prior_ = true;
	        LOG(logDEBUG3) << "could not find detProb, falling back to classifier prior";
	} else {
	        LOG(logDEBUG3) << "COULD find detProb!!";
	}

	if(not use_classifier_prior_ and use_size_dependent_detection_){
	        LOG(logDEBUG3) << "creating detProb feature in traxel store!!";
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

		for(TraxelStore::iterator tr = traxel_store_->begin(); tr != traxel_store_->end(); ++tr) {
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
			traxel_store_->replace(tr, trax);
		}
	}

	LOG(logDEBUG1) << "-> building hypotheses" << endl;
	SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(1, // max_nearest_neighbors
				max_dist_,
				true, // forward_backward
				with_divisions_, // consider_divisions
				division_threshold_
				);
	SingleTimestepTraxel_HypothesesBuilder hyp_builder(traxel_store_, builder_opts);
	hypotheses_graph_ = boost::shared_ptr<HypothesesGraph>(hyp_builder.build());

	hypotheses_graph_->add(arc_distance()).add(tracklet_intern_dist()).add(node_tracklet()).add(tracklet_intern_arc_ids()).add(traxel_arc_id());
 	
	property_map<arc_distance, HypothesesGraph::base_graph>::type& arc_distances = (hypotheses_graph_)->get(arc_distance());
	property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = (hypotheses_graph_)->get(node_traxel());
	bool with_optical_correction = false;
	Traxel some_traxel = (*traxel_map.beginValue());
	if (some_traxel.features.find("com_corrected") != some_traxel.features.end()) {
		LOG(logINFO) << "optical correction enabled";
		with_optical_correction = true;
	}

	for(HypothesesGraph::ArcIt a(*hypotheses_graph_); a!=lemon::INVALID; ++a) {
		HypothesesGraph::Node from = (hypotheses_graph_)->source(a);
		HypothesesGraph::Node to = (hypotheses_graph_)->target(a);
		Traxel from_tr = traxel_map[from];
		Traxel to_tr = traxel_map[to];

		if (with_optical_correction) {
			arc_distances.set(a, from_tr.distance_to_corr(to_tr));
		} else {
			arc_distances.set(a, from_tr.distance_to(to_tr));
		}
	}

        if(event_vector_dump_filename_ != "none")
	  {
	    // store the traxel store and the resulting event vector
	    std::ofstream ofs(event_vector_dump_filename_.c_str());
	    boost::archive::text_oarchive out_archive(ofs);
	    out_archive << ts;
	  }
	return hypotheses_graph_;
    
  }

  std::vector<std::vector<Event> >ConsTracking::track(double forbidden_cost,
						      double ep_gap,
						      bool with_tracklets,
						      double division_weight,
						      double transition_weight,
						      double disappearance_cost,
						      double appearance_cost,
						      bool with_merger_resolution,
						      int n_dim,
						      double transition_parameter,
						      double border_width,
						      bool with_constraints,
						      double cplex_timeout,
						      TimestepIdCoordinateMapPtr coordinates){
    
    LOG(logDEBUG1) <<"max_number_objects  \t"<< max_number_objects_  ; 
    LOG(logDEBUG1) <<"size_dependent_detection_prob\t"<<  use_size_dependent_detection_ ; 
    LOG(logDEBUG1) <<"forbidden_cost\t"<<      forbidden_cost; 
    LOG(logDEBUG1) <<"ep_gap\t"<<      ep_gap; 
    LOG(logDEBUG1) <<"avg_obj_size\t"<<      avg_obj_size_; 
    LOG(logDEBUG1) <<"with_tracklets\t"<<      with_tracklets; 
    LOG(logDEBUG1) <<"division_weight\t"<<      division_weight; 
    LOG(logDEBUG1) <<"transition_weight\t"<<      transition_weight; 
    LOG(logDEBUG1) <<"with_divisions\t"<<      with_divisions_; 
    LOG(logDEBUG1) <<"disappearance_cost\t"<<      disappearance_cost; 
    LOG(logDEBUG1) <<"appearance_cost\t"<<      appearance_cost; 
    LOG(logDEBUG1) <<"with_merger_resolution\t"<<      with_merger_resolution; 
    LOG(logDEBUG1) <<"n_dim\t"<<      n_dim; 
    LOG(logDEBUG1) <<"transition_parameter\t"<<      transition_parameter; 
    LOG(logDEBUG1) <<"border_width\t"<<      border_width; 
    LOG(logDEBUG1) <<"with_constraints\t"<<      with_constraints; 
    LOG(logDEBUG1) <<"cplex_timeout\t"<<      cplex_timeout;
    
    

	double detection_weight = 10;
	Traxels empty;
	boost::function<double(const Traxel&, const size_t)> detection, division;
	boost::function<double(const double)> transition;


	if (use_classifier_prior_) {
		LOG(logINFO) << "Using classifier prior";
		detection = NegLnDetection(detection_weight);
	} else if (use_size_dependent_detection_) {
		LOG(logINFO) << "Using size dependent prior";
		detection = NegLnDetection(detection_weight); // weight 
	} else {
		LOG(logINFO) << "Using hard prior";
		// assume a quasi geometric distribution
		vector<double> prob_vector;
		double p = 0.7; // e.g. for max_number_objects=3, p=0.7: P(X=(0,1,2,3)) = (0.027, 0.7, 0.21, 0.063)
		double sum = 0;
		for(double state = 0; state < max_number_objects_; ++state) {
			double prob = p*pow(1-p,state);
			prob_vector.push_back(prob);
			sum += prob;
		}
		prob_vector.insert(prob_vector.begin(), 1-sum);

		detection = boost::bind<double>(NegLnConstant(detection_weight,prob_vector), _2);
	}

	

	LOG(logDEBUG1) << "division_weight = " << division_weight;
	LOG(logDEBUG1) << "transition_weight = " << transition_weight;
	division = NegLnDivision(division_weight);
	transition = NegLnTransition(transition_weight);

	//border_width_ is given in normalized scale, 1 corresponds to a maximal distance of dim_range/2
	boost::function<double(const Traxel&)> appearance_cost_fn, disappearance_cost_fn;
	LOG(logINFO) << "using border-aware appearance and disappearance costs, with absolute margin: " << border_width;
	appearance_cost_fn = SpatialBorderAwareWeight(appearance_cost,
												border_width,
												false, // true if relative margin to border
												fov_);
	disappearance_cost_fn = SpatialBorderAwareWeight(disappearance_cost,
												border_width,
												false, // true if relative margin to border
												fov_);

	cout << "-> init ConservationTracking reasoner" << endl;
	
	ConservationTracking pgm(
			max_number_objects_,
			detection,
			division,
			transition,
			forbidden_cost,
			ep_gap,
			with_tracklets,
			with_divisions_,
			disappearance_cost_fn,
			appearance_cost_fn,
			true, // with_misdetections_allowed
			true, // with_appearance
			true, // with_disappearance
			transition_parameter,
            with_constraints,
            cplex_timeout
			);

	cout << "-> formulate ConservationTracking model" << endl;
	pgm.formulate(*hypotheses_graph_);

	cout << "-> infer" << endl;
	pgm.infer();

	cout << "-> conclude" << endl;
	pgm.conclude(*hypotheses_graph_);

	cout << "-> storing state of detection vars" << endl;
	last_detections_ = state_of_nodes(*hypotheses_graph_);

	cout << "-> pruning inactive hypotheses" << endl;
	prune_inactive(*hypotheses_graph_);

	cout << "-> constructing unresolved events" << endl;
	boost::shared_ptr<std::vector< std::vector<Event> > > ev = events(*hypotheses_graph_);


	if (max_number_objects_ > 1 && with_merger_resolution && all_true(ev->begin()+1, ev->end(), has_data<Event>)) {
	  cout << "-> resolving mergers" << endl;
	  MergerResolver m(hypotheses_graph_.get());
	  FeatureExtractorBase* extractor;
	  DistanceFromCOMs distance;
	  if (coordinates) {
	    extractor = new FeatureExtractorArmadillo(coordinates);
	  } else {
	    calculate_gmm_beforehand(*hypotheses_graph_, 1, n_dim);
	    extractor = new FeatureExtractorMCOMsFromMCOMs;
	  }
	  FeatureHandlerFromTraxels handler(*extractor, distance);
      
	  m.resolve_mergers(handler);
	  
	  HypothesesGraph g_res;
	  resolve_graph(*hypotheses_graph_, g_res, transition, ep_gap, with_tracklets, transition_parameter, with_constraints);
	  prune_inactive(*hypotheses_graph_);

	  cout << "-> constructing resolved events" << endl;
	  boost::shared_ptr<std::vector< std::vector<Event> > > multi_frame_moves = multi_frame_move_events(*hypotheses_graph_);

	  cout << "-> merging unresolved and resolved events" << endl;
	  // delete extractor; // TO DELETE FIRST CREATE VIRTUAL DTORS
	  ev = merge_event_vectors(*ev, *multi_frame_moves);
	}

	if(event_vector_dump_filename_ != "none")
	  {
	    // store the traxel store and the resulting event vector
	    std::ofstream ofs(event_vector_dump_filename_.c_str());
	    boost::archive::text_oarchive out_archive(ofs);
	    out_archive << *ev;
	  }

	return *ev;

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
