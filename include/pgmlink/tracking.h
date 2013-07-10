/**
@file
@ingroup tracking
@brief tracking API
 */

#ifndef TRACKING_H
#define TRACKING_H

#include "pgmlink/randomforest.h"
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "pgmlink/event.h"
#include "pgmlink/pgmlink_export.h"
#include "pgmlink/traxels.h"
#include "pgmlink/field_of_view.h"

namespace pgmlink {
  class PGMLINK_EXPORT ChaingraphTracking {
  public:
    ChaingraphTracking(const std::string& random_forest_filename = "none",
	      double appearance = 500, 
	      double disappearance = 500,
	      double detection = 10,
	      double misdetection = 500,
	      bool cellness_by_random_forest = false,
	      double opportunity_cost = 0,
	      double forbidden_cost = 0,
	      bool with_constraints = true,
	      bool fixed_detections = false,
	      double mean_div_dist=25,
	      double min_angle=0,
	      double ep_gap=0.01,
	      int n_neighbors=6,
	      bool with_divisions = true,
	      double cplex_timeout = 1e+75,
	      bool alternative_builder=false 
	      )
      : app_(appearance), dis_(disappearance), det_(detection), mis_(misdetection), 
      rf_fn_(random_forest_filename), use_rf_(cellness_by_random_forest), 
      opportunity_cost_(opportunity_cost), forbidden_cost_(forbidden_cost), with_constraints_(with_constraints),
      fixed_detections_(fixed_detections), mean_div_dist_(mean_div_dist), min_angle_(min_angle),
      ep_gap_(ep_gap), n_neighbors_(n_neighbors), with_divisions_(with_divisions),
      cplex_timeout_(cplex_timeout), alternative_builder_(alternative_builder) {}
    std::vector< std::vector<Event> > operator()(TraxelStore&);

    /**
     * Get state of detection variables after call to operator().
     */
    std::vector< std::map<unsigned int, bool> > detections();
    
    /**
     * Setter functions
     */
    void set_with_divisions(bool);
    void set_cplex_timeout(double);

  private:
    double app_, dis_, det_, mis_;
    const std::string rf_fn_;
    bool use_rf_;
    double opportunity_cost_;
    double forbidden_cost_;
    bool with_constraints_;
    bool fixed_detections_;
    double mean_div_dist_, min_angle_;
    double ep_gap_;
    int n_neighbors_;
    bool with_divisions_;
    double cplex_timeout_;
    bool alternative_builder_;
    shared_ptr<std::vector< std::map<unsigned int, bool> > > last_detections_;
  };

}

#endif /* TRACKING_H */
