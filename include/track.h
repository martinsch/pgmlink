#ifndef TRACK_H
#define TRACK_H

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "ilp_construction.h"
#include "randomforest.h"
#include "traxels.h"

namespace Tracking {
  /**
   * Track several timesteps at once.
   */
  class MultiTrack {
  public:
    MultiTrack();
    ~MultiTrack();
    std::vector< std::vector<Event> > operator()(const TraxelStore&);
  private:
    TertiaryEnergy* division_;
    BinaryEnergy* move_;
    BinaryEnergy* mismove_;
    UnaryEnergy* detection_;
    UnaryEnergy* misdetection_;
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
	      double mean_div_dist=25,
	      double min_angle=0)
    : app_(appearance), dis_(disappearance), det_(detection), mis_(misdetection), rf_fn_(random_forest_filename), use_rf_(cellness_by_random_forest), opportunity_cost_(opportunity_cost), mean_div_dist_(mean_div_dist), min_angle_(min_angle) {}
    std::vector< std::vector<Event> > operator()(TraxelStore&);
    
  private:
    double app_, dis_, det_, mis_;
    const std::string rf_fn_;
    bool use_rf_;
    double opportunity_cost_;
    double mean_div_dist_, min_angle_;
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
