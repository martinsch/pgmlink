#ifndef FIELD_OF_VIEW_H
#define FIELD_OF_VIEW_H

#include <vector>


namespace Tracking {

  class FieldOfView {
  public:
  FieldOfView() : lb_(4, 0), ub_(4, 0) {}

  FieldOfView& set_boundingbox(double lt,
		  double lx,
		  double ly,
		  double lz,
		  double ut,
		  double ux,
		  double uy,
		  double uz );

  bool contains( double t, double x, double y, double z ) const;

  const std::vector<double>& lower_bound() const;
  const std::vector<double>& upper_bound() const;

  double spatial_margin( double t, double x, double y, double z ) const;
  double temporal_margin( double t, double x, double y, double z ) const;

  private:
    std::vector<double> lb_; // lower bound
    std::vector<double> ub_; // upper bound
  };

} /* namespace Tracking */

#endif /* FIELD_OF_VIEW_H */
