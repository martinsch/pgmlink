// stl
#include <string>
#include <map>
#include <vector>
#include <stdexcept>

// boost
#include <boost/shared_ptr.hpp>

// vigra
#include <vigra/random_forest.hxx>
#include <vigra/random_forest_hdf5_impex.hxx>

// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/multi_hypotheses_graph.h"
#include "pgmlink/reasoner_multi_hypotheses.h"
#include "pgmlink/pgm_multi_hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/multi_hypotheses_tracking.h"
#include "pgmlink/classifier_auxiliary.h"

// serialization

#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace pgmlink {

////
//// class MultiHypothesesTracking
////
MultiHypothesesTracking::Options::Options() :
    weights(),
    classifiers(),
    paths(),
    feature_lists(),
    with_divisions(false),
    with_constraints(false),
    with_detection_vars(false),
    with_classifiers(false),
    with_constant_classifiers(false),
    with_maximal_conflict_cliques(false),
    forward_backward(false),
    constant_classifier_fallback(false),
    hierarchical_count_factor(false),
    counting_incoming_factor(false),
    classifier_count_precomputed(false),
    with_maximum_arcs(false),
    restrict_timestep_range(false),
  with_one_active_constraint(false),
  with_conflict_factors(false),
  transition_parameter(0)
  {

}

double MultiHypothesesTracking::Options::get_weight(const std::string& name) const {
  std::map<std::string, double>::const_iterator ret = weights.find(name);
  if (ret == weights.end()) {
    throw std::runtime_error("Weight " + name + " not found in MultiHypothesesTracking options!");
  }
  return ret->second;
}

void MultiHypothesesTracking::print_info() {
  LOG(logINFO) << "Calling multihypotheses tracking with the following parameters:\n"
               << "\tappearance: " << options_.weights["app"] << '\n'
               << "\tdisapparance: " << options_.weights["dis"] << '\n'
               << "\tdetection: " << options_.weights["det"] << '\n'
               << "\tmisdetection: " << options_.weights["mis"] << '\n'
               << "\tdvision: " << options_.weights["div"] << '\n'
               << "\tforbidden cost: " << options_.weights["forbidden"] << '\n'
               << "\twith constraints: " << options_.with_constraints << '\n'
               << "\twith divisions: " << options_.with_divisions << '\n'
               << "\twith detection variables: " << options_.with_detection_vars << '\n'
               << "\tnumber of neighbors: " << options_.weights["neighbors"] << '\n'
               << "\tmaximum neighbor distance: " << options_.weights["distance"] << '\n'
               << "\tcplex timeout: " << options_.weights["timeout"] << '\n'
               << "\tep gap: " << options_.weights["gap"] << '\n'
               << "\tmaximum division level: " << options_.weights["max_div"] << '\n'
               << "\twith constant classifiers: " << options_.with_constant_classifiers << '\n'
               << "\twith hierarchical count factor: " << options_.hierarchical_count_factor << '\n'
               << "\twith counting incoming factor: " << options_.counting_incoming_factor << '\n'
               << "\twith maximum arcs: " << options_.with_maximum_arcs << '\n'
               << "\tmaximum arcs: " << options_.weights["arc"] << '\n'
               << "\twith one active region per component constraint: " << options_.with_one_active_constraint  << '\n'    
               << "\twith conflict factors: " << options_.with_conflict_factors << '\n'
               << "\twith transition parameter: " << options_.transition_parameter << '\n'
               << "\topporunity cost: " << options_.get_weight("opportunity")
      ;
}

class my_iarchive : public boost::archive::detail::common_iarchive<my_iarchive>,
                    public boost::archive::detail::shared_ptr_helper {

};

boost::shared_ptr<std::vector<std::vector<Event> > >
MultiHypothesesTracking::operator()(MultiHypothesesTraxelStore& ts, std::string serialize_to_fn) {

  print_info();

  boost::shared_ptr<std::vector<std::vector<Event> > >
      events (new std::vector<std::vector<Event> >);

  SingleTimestepTraxel_HypothesesBuilder::Options builder_options(options_.get_weight("neighbors"), // maximum number of nearest neighbors
                                                                  options_.get_weight("distance"), //distance threshold for nearest neighbors
                                                                  options_.forward_backward, // forward backward,
                                                                  false //options_.with_divisions // dividing objects
                                                                  );
  SingleTimestepTraxel_MultiHypothesesBuilder mult_builder(&ts.ts, builder_options);

  std::cout << " -> workflow: building graph from traxelstore" << std::endl;
  // ts.make_sane();
  MultiHypothesesGraphPtr graph = mult_builder.build_multi_hypotheses_graph();
  graph->add_conflicts(ts.conflicts);
  graph->add_cardinalities();

  if (options_.with_constant_classifiers) {
    LOG(logINFO) << "MultiHypothesesTracking: using constant classifiers";
    ClassifierConstant mov(options_.weights["const_prob"], "move");
    ClassifierConstant div(options_.weights["const_prob"], "division");
    ClassifierConstant det(options_.weights["const_prob"], "detProb");
    ClassifierConstant cnt(options_.weights["const_prob"], "count_prediction");
    graph->add_classifier_features(&mov, &div, &cnt, &det);
  } else if (options_.with_classifiers) {
    LOG(logINFO) << "MultiHypothesesTracking: using classifiers";
    ClassifierStrategyBuilder classifier_builder;
    ClassifierStrategyBuilder::Options classifier_options;
    classifier_options.constant_probability = options_.weights["const_prob"];
    classifier_options.rf_filename = options_.paths["classifier_file"];

    LOG(logINFO) << "MultiHypothesesTracking: creating move classifier";
    classifier_options.name = "move";
    classifier_options.rf_internal_path = options_.paths["classifier_move"];
    classifier_options.type = ClassifierStrategyBuilder::RF_MOVE;
    classifier_options.feature_list.insert(classifier_options.feature_list.end(),
                                           options_.feature_lists["move"].begin(),
                                           options_.feature_lists["move"].end());
    boost::shared_ptr<ClassifierStrategy> mov = classifier_builder.build(classifier_options);
    classifier_options.feature_list.clear();

    LOG(logINFO) << "MultiHypothesesTracking: creating division classifier";
    classifier_options.name = "division";
    classifier_options.rf_internal_path = options_.paths["classifier_division"];
    classifier_options.type = ClassifierStrategyBuilder::RF_DIVISION;
    classifier_options.feature_list.insert(classifier_options.feature_list.end(),
                                           options_.feature_lists["division"].begin(),
                                           options_.feature_lists["division"].end());
    boost::shared_ptr<ClassifierStrategy> div = classifier_builder.build(classifier_options);
    classifier_options.feature_list.clear();

    boost::shared_ptr<ClassifierStrategy> cnt;
    if (options_.classifier_count_precomputed) {
      LOG(logINFO) << "MultiHypothesesTracking: count probability already calculated - using \"lazy\" classifier ";
      cnt = boost::shared_ptr<ClassifierStrategy>(new ClassifierLazy);
    } else {
    LOG(logINFO) << "MultiHypothesesTracking: creating count classifier";
    classifier_options.name = "count_prediction";
    classifier_options.rf_internal_path = options_.paths["classifier_count"];
    classifier_options.type = ClassifierStrategyBuilder::RF_COUNT;
    classifier_options.feature_list.insert(classifier_options.feature_list.end(),
                                           options_.feature_lists["count"].begin(),
                                           options_.feature_lists["count"].end());
    cnt = classifier_builder.build(classifier_options);
    classifier_options.feature_list.clear();
    }
    
    LOG(logINFO) << "MultiHypothesesTracking: creating detection classifier";
    classifier_options.name = "detProb";
    classifier_options.rf_internal_path = options_.paths["classifier_detection"];
    classifier_options.type = ClassifierStrategyBuilder::RF_DETECTION;
    classifier_options.feature_list.insert(classifier_options.feature_list.end(),
                                           options_.feature_lists["detection"].begin(),
                                           options_.feature_lists["detection"].end());
    boost::shared_ptr<ClassifierStrategy> det = classifier_builder.build(classifier_options);
    classifier_options.feature_list.clear();
    
    LOG(logINFO) << "MultiHypothesesTracking: adding classifier features";
    graph->add_classifier_features(mov.get(), div.get(), cnt.get(), det.get());
  }

  graph->remove_traxel_features();

  if (serialize_to_fn.length() != 0) {
	  // get rid of all the traxel features which are not longer needed

	  // create and open a character archive for output
	  std::ofstream ofs(serialize_to_fn.c_str(), std::fstream::out | std::fstream::binary);

	  // save data to archive
	  {
		  boost::archive::text_oarchive oa(ofs);
		  // write class instance to archive
		  oa << *graph;
		// archive and stream closed when destructors are called
	  }
  }

  track(*graph, events);

  return events;
}


boost::shared_ptr<std::vector<std::vector<Event> > >
MultiHypothesesTracking::operator()(std::string deserialize_from_fn) {

  print_info();

  boost::shared_ptr<std::vector<std::vector<Event> > >
      events (new std::vector<std::vector<Event> >);

  std::cout << " -> workflow: deserialize graph" << std::endl;

  MultiHypothesesGraph graph;
  	{
	    // create and open an archive for input
	  	std::ifstream ifs(deserialize_from_fn.c_str(), std::fstream::binary | std::fstream::in);
	  	// read class state from archive
  		boost::archive::text_iarchive ia(ifs);
  		ia >> graph;
  		// archive and stream closed when destructors are called
  	}

  track(graph, events);
  return events;
}

void MultiHypothesesTracking::track(MultiHypothesesGraph& g, boost::shared_ptr<std::vector<std::vector<Event> > >& event_vector) {
  std::cout << " -> workflow: initializing builder" << std::endl;
  
  MultiHypothesesGraph* graph = &g;

  // IMPORTANT FIXME!
  // FIXME: Build proper cost functions for the cases where no classifier priors are available
  // If classifiers are present, the appropriate cost functions will be included in lines 273-295

  ConstantFeature app(options_.get_weight("app"));  // appearance
  ConstantFeature dis(options_.get_weight("dis"));  // disappearance
  ConstantFeature det(options_.get_weight("det"));  // detection
  ConstantFeature mis(options_.get_weight("mis"));  // misdetection
  ConstantFeature div(options_.get_weight("div"));  // division
  ConstantFeature count(options_.get_weight("count")); // count
  SquaredDistance mov; // move

  pgm::multihypotheses::CVPR2014ModelBuilder builder( app, // appearance
                                                      dis, // disappearance,
                                                      mov, // move
                                                      count, // count
                                                      options_.get_weight("forbidden"), // forbidden cost
                                                      options_.get_weight("opportunity"), // opportunity cost
                                                      options_.get_weight("max_div"), // maximum division level
                                                      10000, // max_count -> not neccessary -> fix!
                                                      options_.fov,
                                                      options_.get_weight("border_margin")
  );

  builder.with_maximal_conflict_cliques(options_.with_maximal_conflict_cliques);
  builder.with_hierarchical_counting_factor(options_.hierarchical_count_factor);
  builder.with_counting_incoming_factor(options_.counting_incoming_factor);

  if (options_.with_maximum_arcs) {
    builder.with_maximum_arcs(options_.weights["arc_limit"]);
  }

  if (options_.with_detection_vars) {
    builder.with_detection_vars(det, mis);
  }
  
  if (options_.with_divisions) {
    builder.with_divisions(div);
  }

  if (options_.with_one_active_constraint) {
    builder.with_one_active_per_component_constraint(true);
  }

  if (options_.restrict_timestep_range) {
    builder.with_timestep_range(options_.get_weight("first_timestep"), options_.get_weight("last_timestep"));
    LOG(logINFO) << "MultiHypothesesTracking: restricting timestep range to ["
                 << builder.first_timestep() << ',' << builder.last_timestep() << "] (inclusive)";
  }
  
  if (options_.with_constant_classifiers || options_.with_classifiers) {
    builder
        .with_classifier_priors(NegLnTransition(options_.get_weight("mov")),
                                NegLn(options_.get_weight("count"))
                                );
    if (options_.with_detection_vars) {
      if (options_.with_conflict_factors) {
        LOG(logINFO) << "MultiHypothesesTracking: using conflict factors";
        builder
            .with_conflict_factors(NegLnCardinalityDetection(options_.get_weight("det")));
      } else {
        builder
            .with_detection_vars(NegLnDetection(options_.get_weight("det")),
                                 NegLnDetection(options_.get_weight("det")));
      }
    }
    if (options_.with_divisions) {
      builder.with_divisions(NegLnDivision(options_.get_weight("div"), options_.get_weight("max_div")));
    }
  }
  if (options_.transition_parameter != 0) {
     builder.with_transition_parameter(options_.transition_parameter);
  }

  std::cout << " -> workflow: initializing reasoner" << std::endl;
  MultiHypotheses reasoner(builder,
                           options_.with_constraints,
                           options_.get_weight("gap"), // ep_gap
                           false, // fixed_detections
                           options_.weights["timeout"] // cplex timeout
                           );

  std::cout << " -> workflow: formulating model" << std::endl;
  reasoner.formulate( *graph );

  std::cout << " -> workflow: infer" << std::endl;
  double objective = reasoner.infer();
  LOG(logDEBUG) << "MultiHypothesesTracking(): optimal solution: " << objective;

  std::cout << " -> workflow: conclude" << std::endl;
  reasoner.conclude( *graph );

  prune_inactive( *graph );

  boost::shared_ptr<std::vector<std::vector<Event> > > event_vector_tmp = events( *graph );
  event_vector->push_back(std::vector<Event>()); // need to create empty event at first time frame
  event_vector->insert(event_vector->end(), event_vector_tmp->begin(), event_vector_tmp->end());

 


  std::cout << " -> workflow: return events" << std::endl;

}


}






