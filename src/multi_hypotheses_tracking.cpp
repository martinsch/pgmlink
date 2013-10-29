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


namespace pgmlink {

////
//// class MultiHypothesesTracking
////
double MultiHypothesesTracking::Options::get_weight(const std::string& name) const {
  std::map<std::string, double>::const_iterator ret = weights.find(name);
  if (ret == weights.end()) {
    throw std::runtime_error("Weight " + name + " not found in MultiHypothesesTracking options!");
  }
  return ret->second;
}


boost::shared_ptr<std::vector<std::vector<Event> > >
MultiHypothesesTracking::operator()(MultiHypothesesTraxelStore& ts) {
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
               << "\twith hierarchical count factor: " << options_.hierarchical_count_factor;

  boost::shared_ptr<std::vector<std::vector<Event> > >
      events (new std::vector<std::vector<Event> >);

  MultiHypothesesGraphBuilder mult_builder(MultiHypothesesGraphBuilder::Options(options_.get_weight("neighbors"),
                                                                                options_.get_weight("distance"),
                                                                                options_.forward_backward // forward backward
                                                                                )
                                           );
  std::cout << " -> workflow: building graph from traxelstore" << std::endl;
  ts.make_sane();
  MultiHypothesesGraphPtr graph = mult_builder.build(ts);

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

  std::cout << " -> workflow: initializing builder" << std::endl;
  
  // IMPORTANT FIXME!
  // FIXME: Build proper cost functions!

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
                                                      10000 // max_count -> not neccessary -> fix!
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

  
  if (options_.with_constant_classifiers || options_.with_classifiers) {
    builder
        .with_classifier_priors(NegLnTransition(options_.get_weight("mov")),
                                NegLn(options_.get_weight("count"))
                                );
    if (options_.with_detection_vars) {
      builder
          .with_detection_vars(NegLnDetection(options_.get_weight("det")),
                               NegLnDetection(options_.get_weight("det")));
    }
    if (options_.with_divisions) {
      builder.with_divisions(NegLnDivision(options_.get_weight("div")));
    }
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
  
  std::cout << " -> workflow: creating events" << std::endl;

  // create all neccessary vectors:
  /* for (int i = graph->earliest_timestep(); i < graph->latest_timestep(); ++i) {
    events->push_back(std::vector<Event>());
  } */


  
  MultiHypothesesGraph::ContainedRegionsMap& regions = graph->get(node_regions_in_component());
  MultiHypothesesGraph::node_timestep_map& timesteps = graph->get(node_timestep());

  // CHECK INDICES AT PUSH_BACK!
  // first timestep should be empty!
  events->push_back(std::vector<Event>());
  for (MultiHypothesesGraph::node_timestep_map::ValueIt timestep = timesteps.beginValue();
       timestep != timesteps.endValue();
       ++timestep) {
    // iterating over timesteps. current timestep t == *timestep
    // do nothing in last timestep. there must not be any events there anyway
    MultiHypothesesGraph::node_timestep_map::ValueIt timestep_check = timestep;
    if (++timestep_check == timesteps.endValue()) {
      break;
    }
    events->push_back(std::vector<Event>());
    for (MultiHypothesesGraph::node_timestep_map::ItemIt n(timesteps, *timestep);
         n!= lemon::INVALID;
         ++n) {
      std::vector<Traxel>& traxels = regions.get_value(n);
      LOG(logDEBUG4) << "Region " << traxels[0].Id << " at time " << traxels[0].Timestep << '\n';
      for (std::vector<Traxel>::iterator t = traxels.begin(); t != traxels.end(); ++t) {
        assert(t->Timestep == *timestep);
        std::vector<Event>& events_at = *(events->rbegin());
        if (t->features["active"][0] > 0.) {        
          if (t->features["outgoing"].size() == 2) {
            // division: parent cell in timestep t, children cells in t+1
            Event e;
            e.type = Event::Division;
            e.traxel_ids.push_back(t->Id);
            e.traxel_ids.push_back(t->features["outgoing"][0]);
            e.traxel_ids.push_back(t->features["outgoing"][1]);
            events_at.push_back(e);
          } else if (t->features["outgoing"].size() == 1) {
            Event e;
            e.type = Event::Move;
            e.traxel_ids.push_back(t->Id);
            e.traxel_ids.push_back(t->features["outgoing"][0]);
            events_at.push_back(e);
          } else if (t->features["outgoing"].size() == 0 && t->Timestep < graph->latest_timestep()) {
            Event e;
            e.type = Event::Disappearance;
            e.traxel_ids.push_back(t->Id);
            events_at.push_back(e);
          }
        }
      }
    }
    // appears in the next timestep
    for (MultiHypothesesGraph::node_timestep_map::ItemIt n(timesteps, *timestep + 1);
         n != lemon::INVALID;
         ++n) {
      std::vector<Traxel>& traxels = regions.get_value(n);
      for (std::vector<Traxel>::iterator t = traxels.begin(); t != traxels.end(); ++t) {
        std::vector<Event>& events_at = *(events->rbegin());
        if (t->features["active"][0] > 0.) {
          if (t->features["parent"].size() == 0 && t->Timestep > graph->earliest_timestep()) {
            Event e;
            e.type = Event::Appearance;
            e.traxel_ids.push_back(t->Id);
            events_at.push_back(e);
          }
        }
      }
    }
  }


  std::cout << " -> workflow: return events" << std::endl;
  return events;
}


}






