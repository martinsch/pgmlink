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
#include "pgmlink/traxels.h"
#include "pgmlink/field_of_view.h"

namespace pgmlink {
  class ChaingraphTracking {
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

  class NNTracking {
    public:
      NNTracking(double divDist = 30,
    		  double movDist = 10,
    		  std::vector<std::string> features = 0,
    		  double divisionThreshold = 0.5,
    		  bool splitterHandling = true,
    		  bool mergerHandling = true,
    		  std::vector<int> maxTraxelIdAt = 0)
        : divDist_(divDist), movDist_(movDist), distanceFeatures_(features),
          divisionThreshold_(divisionThreshold), splitterHandling_(splitterHandling),
          mergerHandling_(mergerHandling), maxTraxelIdAt_(maxTraxelIdAt){}
      std::vector< std::vector<Event> > operator()(TraxelStore&);

      /**
       * Get state of detection variables after call to operator().
       */
      std::vector< std::map<unsigned int, bool> > detections();

    private:
      double divDist_, movDist_;
      std::vector<std::string> distanceFeatures_;
      double divisionThreshold_;
      bool splitterHandling_, mergerHandling_;
      shared_ptr<std::vector< std::map<unsigned int, bool> > > last_detections_;
      std::vector<int> maxTraxelIdAt_;
    };


  class NNTrackletsTracking {
      public:
	  NNTrackletsTracking(double maxDist = 30,
      		  std::vector<std::string> features = 0,
      		  double divisionThreshold = 0.5,
      		  bool splitterHandling = true,
      		  bool mergerHandling = true,
      		  std::vector<int> maxTraxelIdAt = 0)
          : maxDist_(maxDist), distanceFeatures_(features),
            divisionThreshold_(divisionThreshold), splitterHandling_(splitterHandling),
            mergerHandling_(mergerHandling), maxTraxelIdAt_(maxTraxelIdAt){}
        std::vector< std::vector<Event> > operator()(TraxelStore&);

        /**
         * Get state of detection variables after call to operator().
         */
        std::vector< std::map<unsigned int, bool> > detections();

      private:
        double maxDist_;
        std::vector<std::string> distanceFeatures_;
        double divisionThreshold_;
        bool splitterHandling_, mergerHandling_;
        shared_ptr<std::vector< std::map<unsigned int, bool> > > last_detections_;
        std::vector<int> maxTraxelIdAt_;
  };

  class ConsTracking {
    public:
	  ConsTracking(
		  int max_number_objects=3,
	      double max_neighbor_distance = 20,
		  double division_threshold = 0.3,
		  const std::string& random_forest_filename = "none",
  	      bool size_dependent_detection_prob = false,
  	      double forbidden_cost = 0,
  	      double ep_gap=0.01,
  	      double avg_obj_size=30.0,
  	      bool with_tracklets=true,
  	      double division_weight=10.0,
  	      double transition_weight=10.0,
  	      bool with_divisions=true,
  	      double disappearance_cost = 0,
  	      double appearance_cost = 0,
		  bool with_merger_resolution = true,
		  int n_dim = 3,
		  double transition_parameter = 5.
  	      )
        : max_number_objects_(max_number_objects),
        	max_dist_(max_neighbor_distance), division_threshold_(division_threshold),
        detection_rf_fn_(random_forest_filename), use_size_dependent_detection_(size_dependent_detection_prob),
        forbidden_cost_(forbidden_cost),
        ep_gap_(ep_gap), avg_obj_size_(avg_obj_size),
        with_tracklets_(with_tracklets),
        division_weight_(division_weight),
        transition_weight_(transition_weight),
        with_divisions_(with_divisions),
        disappearance_cost_(disappearance_cost),
        appearance_cost_(appearance_cost),
		means_(std::vector<double>()),
		sigmas_(std::vector<double>()),
		with_merger_resolution_(with_merger_resolution),
		number_of_dimensions_(n_dim),
		transition_parameter_(transition_parameter){}
      std::vector< std::vector<Event> > operator()(TraxelStore&);

      /**
       * Get state of detection variables after call to operator().
       */
      std::vector< std::map<unsigned int, bool> > detections();

    private:
      int max_number_objects_;
      double max_dist_;
      double division_threshold_;
      const std::string detection_rf_fn_;
      bool use_size_dependent_detection_;
      double forbidden_cost_;
      double ep_gap_;
      double avg_obj_size_;
      bool with_tracklets_;
      double division_weight_;
      double transition_weight_;
      bool with_divisions_;
      double disappearance_cost_, appearance_cost_;
      std::vector<double> means_, sigmas_;
      bool with_merger_resolution_;
      shared_ptr<std::vector< std::map<unsigned int, bool> > > last_detections_;
      int number_of_dimensions_;
      double transition_parameter_;
    };
}

#endif /* TRACKING_H */
