/**
   @file
   @ingroup tracking
   @brief feature functions
*/

#ifndef FEATURE_H
#define FEATURE_H

#include <cmath>
#include <stdexcept>

#include "pgmlink/traxels.h"
#include "pgmlink/pgmlink_export.h"


namespace pgmlink {

  class PGMLINK_EXPORT GeometryDivision2 {
    /**
     * Division feature based on geometric properties of a ancestor-child-child configuration.
     *
     * The feature is the sum of the squared distances from the ancestor object to the two children shifted by
     * a mean distance [distance dependence].
     * Furthermore, there may be an additional hard constraint on the division angle. If the angle is to small,
     * the feature will be set to a very high value. [angle constraint]
     */
  public:
    /**
     * @param mean_div_dist expected moving distance during divisions
     * @param min_angle minimal angle to accept configuration as a division
     * @param distance_dependence if turned off, feature will not depend on distance between ancestor and children
     * @param angle_constraint if turned on, a minimal division angle is demanded
     */
  GeometryDivision2(double mean_div_dist, double min_angle, bool distance_dependence=true, bool angle_constraint=true) 
    : mean_div_dist_(mean_div_dist), min_angle_(min_angle), 
      distance_dependence_(distance_dependence), angle_constraint_(angle_constraint) {};
    double operator()(const Traxel& ancestor,
		      const Traxel& child1,
		      const Traxel& child2) const;
  private:
    double mean_div_dist_, min_angle_;
    bool distance_dependence_, angle_constraint_; 
  };

  class PGMLINK_EXPORT KasterDivision2 {
  public:
  KasterDivision2(double weight, double div_cost) : w_(weight), div_cost_(div_cost) {};
    double operator()(const Traxel& ancestor,
		      const Traxel& child1,
		      const Traxel& child2) const;
  private:
    double w_;
    double div_cost_;
  };

  class PGMLINK_EXPORT NegLnCellness {
  public:
  NegLnCellness(double weight) : w_(weight) {}
    double operator()( const Traxel& ) const;
  private:
    double w_;
  };

  class PGMLINK_EXPORT NegLnOneMinusCellness {
  public:
  NegLnOneMinusCellness(double weight) : w_(weight) {}
    double operator()( const Traxel& ) const;
  private:
    double w_;
  };
 
  /**
     @brief Zero near temporal border, else 1.
  */
  class PGMLINK_EXPORT BorderAwareConstant {
  public:
  BorderAwareConstant( double weight,
		       int at,
		       bool early,
		       int margin_t=1) 
    : w_(weight), at_(at), early_(early), margin_t_(margin_t) {}
   
    double operator()( const Traxel& tr ) const {
      int t = tr.Timestep;
      if(early_) {
	if( at_ <= t && t < at_ + margin_t_) { 
	  return 0;
	} else return w_;
      } else {
	if (at_ - margin_t_ < t && t <= at_) {
	  return 0.;
	} else return w_;
      }
    }

  private:
    double w_;
    int at_;
    bool early_;
    int margin_t_;
  };

  /**
     @brief Primitive fixed value feature.
  */
  class PGMLINK_EXPORT ConstantFeature {
  public:
  ConstantFeature( double value = 0. ) : value(value) {}
   
    double operator()( const Traxel&,
		       const Traxel&,
		       const Traxel& ) const { return value; }
    double operator()( const Traxel&,
		       const Traxel&) const { return value; }
    double operator()( const Traxel& ) const { return value; }
    double operator()() const { return value; };
   
    double value;
  };
    
  class PGMLINK_EXPORT SquaredDistance {
  public:
    double operator()(const Traxel& from, const Traxel& to) const;
  };

  class PGMLINK_EXPORT KasterDivision {
  public:
  KasterDivision(double division_cost) : div_cost_(division_cost) {};

    double operator()(const Traxel& ancestor,
		      const Traxel& child1,
		      const Traxel& child2) const;
  private:
    double div_cost_;
  };

  class PGMLINK_EXPORT GeometryDivision {
  public:
    double operator()(const Traxel& ancestor,
		      const Traxel& child1,
		      const Traxel& child2) const;
  };

  ////
  //// Cellness based mlinder-type energies
  ////
  class PGMLINK_EXPORT CellnessDivision {
  public:
    CellnessDivision( double diffCellness = 1461, double absCellness = 190 );
    
    double operator()(const Traxel& ancestor,
		      const Traxel& child1,
		      const Traxel& child2) const;
  private:
    // weight for the difference of daughter cells' cellness
    double param_diff_c;
    // weight for the absolute cellness of the parent cell
    double param_abs_c;
  };
  
  class PGMLINK_EXPORT CellnessMove {
  public:
    CellnessMove( double diffCellness = 1140 );
    
    double operator()(const Traxel& from,
		      const Traxel& to) const;

  private:
    // weight for the difference of the cells' cellness
    double param_diff_c;
  };

  class PGMLINK_EXPORT CellnessDisappearance {
  public:
    CellnessDisappearance( double disappWeight = 1000 );
    
    double operator()(const Traxel& from) const;

  private:
    // weight for the absolute cellness
    double param_abs_c;
  };

  class PGMLINK_EXPORT CellnessAppearance {
  public:
    CellnessAppearance( double appWeight = 1000 );

    double operator()(const Traxel& to) const;

  private:
    // weight for the absolute cellness
    double param_abs_c;
  };
  
} /* Namespace pgmlink */

#endif /* FEATURE_H */
