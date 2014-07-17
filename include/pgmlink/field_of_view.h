/**
   @file
   @ingroup tracking
   @brief field-of-view filtering of a traxelstore 
*/

#ifndef FIELD_OF_VIEW_H
#define FIELD_OF_VIEW_H

#include <vector>
#include "pgmlink/pgmlink_export.h"



namespace pgmlink {

  /** Field of view in 3d+t space as a rectangular cuboid. */
  class FieldOfView 
  {
   public:
    PGMLINK_EXPORT FieldOfView() 
    : lb_(4, 0), ub_(4, 0)
    {}

    PGMLINK_EXPORT FieldOfView(double lt,
		                       double lx,
		                       double ly,
		                       double lz,
		                       double ut,
		                       double ux,
		                       double uy,
		                       double uz ) 
    : lb_(4, 0), ub_(4, 0) 
    {
      set_boundingbox(lt, lx, ly, lz, ut, ux, uy, uz);
    }

    /**
     * Set lower and upper bound of cuboid.
     * One has for all coordinates l. <= u.
     */
    PGMLINK_EXPORT FieldOfView& set_boundingbox(double lt,
		                         double lx,
		                         double ly,
		                         double lz,
		                         double ut,
		                         double ux,
		                         double uy,
		                         double uz );

    PGMLINK_EXPORT bool contains( double t, double x, double y, double z ) const;
     
    PGMLINK_EXPORT const std::vector<double>& lower_bound() const;
    PGMLINK_EXPORT const std::vector<double>& upper_bound() const;

    /**
     * Distance to nearest face of the cuboid.
     * The point can be inside or outside the field of view.
     */
    PGMLINK_EXPORT double spatial_margin( double t, double x, double y, double z ) const;
    PGMLINK_EXPORT double spatial_distance_to_border( double t, double x, double y, double z, bool relative ) const;
    
    /** Shortest distance to the temporal boundary of the field of view. */
    PGMLINK_EXPORT double temporal_margin( double t, double x, double y, double z ) const;

  private:
    std::vector<double> lb_; // lower bound
    std::vector<double> ub_; // upper bound
  };

} /* namespace pgmlink */

#endif /* FIELD_OF_VIEW_H */
