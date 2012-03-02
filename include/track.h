#ifndef TRACK_H
#define TRACK_H

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "ilp_construction.h"
#include "randomforest.h"
#include "traxels.h"
#include "field_of_view.h"

namespace Tracking {
  class BotTracking {
  public:
    BotTracking(double detection = 10,
	      double misdetection = 500,
	      double opportunity_cost = 0) : det_(detection), mis_(misdetection), opportunity_cost_(opportunity_cost) {}
    std::vector< std::vector<Event> > operator()(TraxelStore&);
    
  private:
    double det_, mis_;
    double opportunity_cost_;
  };

  class MrfTracking {
  public:
    MrfTracking(const std::string& random_forest_filename = "none",
	      double appearance = 500, 
	      double disappearance = 500,
	      double detection = 10,
	      double misdetection = 500,
	      bool cellness_by_random_forest = false,
	      double opportunity_cost = 0,
	      double forbidden_cost = 0,
	      bool with_constraints = true,
	      double mean_div_dist=25,
	      double min_angle=0)
      : app_(appearance), dis_(disappearance), det_(detection), mis_(misdetection), 
      rf_fn_(random_forest_filename), use_rf_(cellness_by_random_forest), 
      opportunity_cost_(opportunity_cost), forbidden_cost_(forbidden_cost), with_constraints_(with_constraints), mean_div_dist_(mean_div_dist), min_angle_(min_angle) {}
    std::vector< std::vector<Event> > operator()(TraxelStore&);

    /**
     * Get state of detection variables after call to operator().
     */
    std::vector< std::map<unsigned int, bool> > detections();
    
  private:
    double app_, dis_, det_, mis_;
    const std::string rf_fn_;
    bool use_rf_;
    double opportunity_cost_;
    double forbidden_cost_;
    bool with_constraints_;
    double mean_div_dist_, min_angle_;
    shared_ptr<std::vector< std::map<unsigned int, bool> > > last_detections_; 
  };



  class KanadeTracking {
  public:
  KanadeTracking(FieldOfView fov = FieldOfView(0,0,0,0, 100,4000,4000,4000),
		 double misdetection_rate = 0.1265,
		 double temporal_lambda = 5,
		 double spatial_lambda = 30,
		 double link_lambda = 25,
		 double temporal_cutoff = 15,
		 double spatial_cutoff = 40) 
    : fov_(fov),
      misdetection_rate_(misdetection_rate),
      temporal_lambda_(temporal_lambda),
      spatial_lambda_(spatial_lambda),
      link_lambda_(link_lambda),
      temporal_cutoff_(temporal_cutoff),
      spatial_cutoff_(spatial_cutoff) {}

    std::vector< std::vector<Event> > operator()(TraxelStore&);

    /**
     * Get state of detection variables after call to operator().
     */
    std::vector< std::map<unsigned int, bool> > detections();

  private:
    FieldOfView fov_;
    double misdetection_rate_;
    double temporal_lambda_;
    double spatial_lambda_;
    double link_lambda_;
    double temporal_cutoff_;
    double spatial_cutoff_;

    shared_ptr<std::vector< std::map<unsigned int, bool> > > last_detections_; 
  };



     /**
     * Track two timesteps with an adaptive energies-type ilp.
     */
    class Track {
    public:
      Track() {};
      Track(const AdaptiveEnergiesFormulation& f) : formulation_(f) {};

      std::vector<Event> operator()(const Traxels& prev, const Traxels& curr);

      // getter/setter
      Track& Formulation(const AdaptiveEnergiesFormulation& f);
      const AdaptiveEnergiesFormulation& Formulation() const;

    protected:
      boost::shared_ptr<std::vector<Event> > events_from(const Traxels& prev, const Traxels& curr,
							 boost::shared_array<double> ilp_finalVars,
							 int ilp_nVars,
							 const std::vector<Event>& ilp_events) const;

    private:
      AdaptiveEnergiesFormulation formulation_;
    };

    

    class FixedCostTracking {
    public:
	FixedCostTracking(double divison, double move, double disappearance, double appearance,
	    double distance_threshold = 50);
	std::vector<Event> operator()(const Traxels& prev, const Traxels& curr);

    private:
	AdaptiveEnergiesFormulation formulation_;
    };

    class ShortestDistanceTracking {
    public:
	ShortestDistanceTracking(double divison, double disappearance, double appearance,
				 double distance_threshold = 50,
				 unsigned int max_nearest_neighbors = 6);
	std::vector<Event> operator()(const Traxels& prev, const Traxels& curr);

    private:
	AdaptiveEnergiesFormulation formulation_;
    };

    

    class CellnessTracking {
    public:
        // Weights description:
        //
        // Division:
        //   w_div1: weight for the difference of the children's cellness
        //   w_div2: weight for the parent cellness
        //
        // Move:
        //   w_move: weight for the cellness difference
        //
        // Appearance:
        //   w_app: scale of the appearance cost
        //
        // Disappearance:
        //   w_disapp: scale of the disappearance cost
        CellnessTracking(const std::string& random_forest_file,
                         double w_div1 = 1461,
                         double w_div2 = 190,
                         double w_move = 1140,
                         double w_app = 1000,
                         double w_disapp = 1000,
                         double distance_threshold = 50
                         );
        std::vector<Event> operator()(Traxels prev, Traxels curr);
    private:
        AdaptiveEnergiesFormulation formulation_;
        vigra::RandomForest<RF::RF_LABEL_TYPE> random_forest_;
        std::vector<std::string> rf_features_;
    };

}

#endif /* TRACK_H */
