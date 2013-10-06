// stl
#include <string>
#include <map>
#include <vector>
#include <stdexcept>

// boost
#include <boost/shared_ptr.hpp>

// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/multi_hypotheses_graph.h"
#include "pgmlink/reasoner_multi_hypotheses.h"
#include "pgmlink/pgm_multi_hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/multi_hypotheses_tracking.h"


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
               << "\tmaximum division level " << options_.weights["max_div"];

  boost::shared_ptr<std::vector<std::vector<Event> > >
      events (new std::vector<std::vector<Event> >);

  MultiHypothesesGraphBuilder mult_builder(MultiHypothesesGraphBuilder::Options(options_.get_weight("neighbors"),
                                                                                options_.get_weight("distance"),
                                                                                false // forward backward
                                                                                )
                                           );
  std::cout << " -> workflow: building graph from traxelstore" << std::endl;
  MultiHypothesesGraphPtr graph = mult_builder.build(ts);

  std::cout << " -> workflow: initializing builder" << std::endl;
  
  // IMPORTANT FIXME!
  // FIXME: Build proper cost functions!n

  ConstantFeature app(options_.get_weight("app"));  // appearance
  ConstantFeature dis(options_.get_weight("dis"));  // disappearance
  ConstantFeature det(options_.get_weight("det"));  // detection
  ConstantFeature mis(options_.get_weight("mis"));  // misdetection
  ConstantFeature div(options_.get_weight("div"));  // division
  SquaredDistance mov; // move

  pgm::multihypotheses::CVPR2014ModelBuilder builder( app, // appearance
                                                      dis, // disappearance,
                                                      mov, // move
                                                      options_.get_weight("forbidden"), // forbidden cost
                                                      options_.get_weight("max_div"), // maximum division level
                                                      10000 // max_count -> not neccessary -> fix!
  );
  if (options_.with_detection_vars) {
    builder.with_detection_vars(det, mis);
  }
  if (options_.with_divisions) {
    builder.with_divisions(div);
  }

  std::cout << " -> workflow: initializing reasoner" << std::endl;
  MultiHypotheses reasoner(builder,
                           options_.with_constraints,
                           options_.get_weight("gap"), // ep_gap
                           false, // fixed_detections
                           options_.weights["timeout"] // cplex timeout
                           );

  std::cout << " -> workflow: formulationg model" << std::endl;
  reasoner.formulate( *graph );

  std::cout << " -> workflow: infer" << std::endl;
  double objective = reasoner.infer();
  LOG(logDEBUG) << "MultiHypothesesTracking(): optimal solution: " << objective;

  std::cout << " -> workflow: conclude" << std::endl;
  reasoner.conclude( *graph );
  
  std::cout << " -> workflow: creating events" << std::endl;

  // create all neccessary vectors:
  for (int i = graph->earliest_timestep(); i < graph->latest_timestep(); ++i) {
    events->push_back(std::vector<Event>());
  }
  
  MultiHypothesesGraph::ContainedRegionsMap& regions = graph->get(node_regions_in_component());

  // CHECK INDICES AT PUSH_BACK!
  for (MultiHypothesesGraph::NodeIt n(*graph); n != lemon::INVALID; ++n) {
    std::vector<Traxel>& traxels = regions.get_value(n);
    LOG(logDEBUG4) << "Region " << traxels[0].Id << " at time " << traxels[0].Timestep << '\n';
    for (std::vector<Traxel>::iterator t = traxels.begin(); t != traxels.end(); ++t) {
      break; // no event creation for now!
      if (t->features["active"][0] > 0.) {        
        if (t->features["outgoing"].size() == 2) {
          Event e;
          e.type = Event::Division;
          e.traxel_ids.push_back(t->Id);
          e.traxel_ids.push_back(t->features["outgoing"][0]);
          e.traxel_ids.push_back(t->features["outgoing"][1]);
          (*events)[t->Timestep-graph->earliest_timestep()].push_back(e);
        } else if (t->features["outgoing"].size() == 1) {
          Event e;
          e.type = Event::Move;
          e.traxel_ids.push_back(t->Id);
          e.traxel_ids.push_back(t->features["outgoing"][0]);
          (*events)[t->Timestep-graph->earliest_timestep()].push_back(e);
        } else if (t->features["outgoing"].size() == 0 && t->Timestep < graph->latest_timestep()) {
          Event e;
          e.type = Event::Disappearance;
          e.traxel_ids.push_back(t->Id);
          (*events)[t->Timestep-graph->earliest_timestep()].push_back(e);
        }
        if (t->features["incoming"].size() == 0 && t->Timestep > graph->earliest_timestep()) {
          Event e;
          e.type = Event::Appearance;
          e.traxel_ids.push_back(t->Id);
          (*events)[t->Timestep-graph->earliest_timestep()].push_back(e);
        }
      }
    }
  }

  std::cout << " -> workflow: return events" << std::endl;
  return events;
}


}






