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
    double get_weight(const std::string& name) const;
    std::map<std::string, double> weights;
    std::map<std::string, vigra::RandomForest<> > classifiers;
    bool with_divisions;
    bool with_constraints;
    bool with_detection_vars;
    bool with_classifiers;
    bool with_constant_classifiers;
    bool with_maximal_conflict_cliques;
  };
  MultiHypothesesTracking(const Options& options) : options_(options) {}
  boost::shared_ptr<std::vector<std::vector<Event> > > operator()(MultiHypothesesTraxelStore& ts);
private:
  Options options_;
};
}
  

#endif /* MULTI_HYPOTHESES_TRACKING_H */
