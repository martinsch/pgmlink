/**
@file
@ingroup multi_hypotheses_tracking
@brief multi hypotheses tracking API
 */

#ifndef MULTI_HYPOTHESES_TRACKING_H
#define MULTI_HYPOTHESES_TRACKING_H

// stl
#include <map>
#include <string>
#include <vector>

// boost
#include <boost/shared_ptr.hpp>

// vigra
#include <vigra/random_forest.hxx>

// pgmlink
#include "pgmlink/event.h"
#include "pgmlink/traxels.h"
#include "pgmlink/pgmlink_export.h"
#include "pgmlink/multi_hypotheses_graph.h"


namespace pgmlink {

////
//// class MultiHypothesesTracking
////
class PGMLINK_EXPORT MultiHypothesesTracking {
public:
  struct Options {
    Options();
    double get_weight(const std::string& name) const;
    std::map<std::string, double> weights;
    std::map<std::string, vigra::RandomForest<> > classifiers;
    std::map<std::string, std::string> paths;
    std::map<std::string, std::vector<std::pair<std::string, std::string> > > feature_lists;
    bool with_divisions;
    bool with_constraints;
    bool with_detection_vars;
    bool with_classifiers;
    bool with_constant_classifiers;
    bool with_maximal_conflict_cliques;
    bool forward_backward;
    bool constant_classifier_fallback;
    bool hierarchical_count_factor;
    bool counting_incoming_factor;
    bool classifier_count_precomputed;
    bool with_maximum_arcs;
    bool restrict_timestep_range;
    bool with_one_active_constraint;
    bool with_conflict_factors;
  };
  MultiHypothesesTracking(const Options& options) : options_(options) {}
  boost::shared_ptr<std::vector<std::vector<Event> > > operator()(MultiHypothesesTraxelStore& ts, std::string serialize_to_fn="");
  boost::shared_ptr<std::vector<std::vector<Event> > > operator()(std::string deserialize_from_fn);
private:
  void print_info();
  void track(MultiHypothesesGraph& graph, boost::shared_ptr<std::vector<std::vector<Event> > >& events);
  Options options_;
};
}
  

#endif /* MULTI_HYPOTHESES_TRACKING_H */
