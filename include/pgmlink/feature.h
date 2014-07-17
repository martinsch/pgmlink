/**
   @file
   @ingroup tracking
   @brief feature functions
*/

#ifndef FEATURE_H
#define FEATURE_H

#include <cmath>
#include <stdexcept>
#include "pgmlink/log.h"
#include "pgmlink/traxels.h"
#include "pgmlink/field_of_view.h"
#include "pgmlink/pgmlink_export.h"

namespace pgmlink {

class GeometryDivision2 
{
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
    PGMLINK_EXPORT
    GeometryDivision2(double mean_div_dist, double min_angle, bool distance_dependence=true, bool angle_constraint=true) 
    : mean_div_dist_(mean_div_dist), min_angle_(min_angle), 
      distance_dependence_(distance_dependence), angle_constraint_(angle_constraint)
    {}

    PGMLINK_EXPORT double operator()(const Traxel& ancestor,
                                     const Traxel& child1,
                                     const Traxel& child2) const;
  private:
    double mean_div_dist_, min_angle_;
    bool distance_dependence_, angle_constraint_; 
};

class KasterDivision2
{
  public:
    PGMLINK_EXPORT KasterDivision2(double weight, double div_cost)
    : w_(weight), div_cost_(div_cost)
    {}
    
    PGMLINK_EXPORT double operator()(const Traxel& ancestor,
                                     const Traxel& child1,
                                     const Traxel& child2) const;
  private:
    double w_;
    double div_cost_;
};

class NegLnCellness 
{
  public:
    PGMLINK_EXPORT NegLnCellness(double weight) 
    : w_(weight)
    {}
    
    PGMLINK_EXPORT double operator()( const Traxel& ) const;
  private:
    double w_;
};

class NegLnOneMinusCellness
{
  public:
    PGMLINK_EXPORT NegLnOneMinusCellness(double weight) : w_(weight) {}
    PGMLINK_EXPORT double operator()( const Traxel& ) const;
  private:
    double w_;
};
 
class NegLnDetection 
{
public:
    PGMLINK_EXPORT NegLnDetection(double weight)
    : w_(weight)
    {}
    
    PGMLINK_EXPORT double operator()( const Traxel&, const size_t state ) const;
private:
    double w_;
};

class NegLnConstant 
{
 public:
    PGMLINK_EXPORT NegLnConstant(double weight, std::vector<double> prob_vector)
    : w_(weight), prob_vector_(prob_vector)
    {}

    PGMLINK_EXPORT double operator()(const size_t state ) const;
private:
    double w_;
    std::vector<double> prob_vector_;
};

class NegLnDivision 
{
 public:
    PGMLINK_EXPORT NegLnDivision(double weight) 
    : w_(weight) 
    {}
    
    PGMLINK_EXPORT double operator()( const Traxel&, const size_t state ) const;
private:
    double w_;
};

class NegLnTransition 
{
 public:
    PGMLINK_EXPORT NegLnTransition(double weight)
    : w_(weight)
    {}

    PGMLINK_EXPORT double operator()( const double ) const;
private:
    double w_;
};

/**
   @brief Zero near temporal border, else 1.
*/
class BorderAwareConstant 
{
 public:
  PGMLINK_EXPORT BorderAwareConstant( double weight,
                                      int at,
                                      bool early,
                                      int margin_t=1) 
  : w_(weight), at_(at), early_(early), margin_t_(margin_t)
  {}
 
  PGMLINK_EXPORT double operator()( const Traxel& tr ) const 
  {
    int t = tr.Timestep;
    if(early_) 
    {
      if( at_ <= t && t < at_ + margin_t_) 
      { 
        return 0;
      } 
      else 
        return w_;
    } else {
      if (at_ - margin_t_ < t && t <= at_) 
      {
        return 0.;
      } 
      else 
        return w_;
    }
  }

private:
  double w_;
  int at_;
  bool early_;
  int margin_t_;
};

class SpatialBorderAwareWeight
{
  public:
    PGMLINK_EXPORT SpatialBorderAwareWeight( double cost, double margin, bool relative, FieldOfView& fov) 
    : cost_(cost), margin_(margin), relative_(relative), fov_(fov)
    {
        if (relative && margin > 0.5) {
            throw std::runtime_error("The relative margin may not exceed 0.5.");
        }
    }

    PGMLINK_EXPORT double operator()( const Traxel& tr ) const;

  private:
    double cost_;
    double margin_;
    bool relative_;
    FieldOfView fov_;
  };


/**
   @brief Primitive fixed value feature.
*/
class ConstantFeature 
{
 public:
  PGMLINK_EXPORT ConstantFeature( double value = 0. ) 
  : value(value) 
  {}
 
  PGMLINK_EXPORT double operator()( const Traxel&,
                                    const Traxel&,
                                    const Traxel& ) const { return value; }

  PGMLINK_EXPORT double operator()( const Traxel&,
                                    const Traxel&) const { return value; }

  PGMLINK_EXPORT double operator()( const Traxel& ) const { return value; }

  PGMLINK_EXPORT double operator()() const { return value; };
 
  double value;
};
  
class SquaredDistance 
{
public:
  PGMLINK_EXPORT double operator()(const Traxel& from, const Traxel& to) const;
};

class KasterDivision 
{
 public:
  PGMLINK_EXPORT KasterDivision(double division_cost)
  : div_cost_(division_cost)
  {}

  PGMLINK_EXPORT double operator()(const Traxel& ancestor,
                                   const Traxel& child1,
                                   const Traxel& child2) const;
private:
  double div_cost_;
};

class GeometryDivision 
{
 public:
  PGMLINK_EXPORT double operator()(const Traxel& ancestor,
                                   const Traxel& child1,
                                   const Traxel& child2) const;
};

////
//// Cellness based mlinder-type energies
////
class CellnessDivision 
{
public:
  PGMLINK_EXPORT CellnessDivision( double diffCellness = 1461, double absCellness = 190 );
  
  PGMLINK_EXPORT double operator()(const Traxel& ancestor,
                                   const Traxel& child1,
                                   const Traxel& child2) const;
private:
  // weight for the difference of daughter cells' cellness
  double param_diff_c;
  // weight for the absolute cellness of the parent cell
  double param_abs_c;
};

class CellnessMove 
{
 public:
  PGMLINK_EXPORT CellnessMove( double diffCellness = 1140 );
  
  PGMLINK_EXPORT double operator()(const Traxel& from,
                                   const Traxel& to) const;

private:
  // weight for the difference of the cells' cellness
  double param_diff_c;
};

class CellnessDisappearance 
{
 public:
  PGMLINK_EXPORT CellnessDisappearance( double disappWeight = 1000 );
  
  PGMLINK_EXPORT double operator()(const Traxel& from) const;

private:
  // weight for the absolute cellness
  double param_abs_c;
};

class CellnessAppearance 
{
 public:
  PGMLINK_EXPORT CellnessAppearance( double appWeight = 1000 );

  PGMLINK_EXPORT double operator()(const Traxel& to) const;

private:
  // weight for the absolute cellness
  double param_abs_c;
};
  
} /* Namespace pgmlink */

#endif /* FEATURE_H */
