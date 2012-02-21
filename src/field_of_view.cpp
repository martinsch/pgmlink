#include "field_of_view.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>

using namespace std;

namespace Tracking {

  FieldOfView& FieldOfView::set_boundingbox(double lt,
					    double lx,
					    double ly,
					    double lz,
					    double ut,
					    double ux,
					    double uy,
					    double uz ){
    if(lt>ut || lx>ux || ly>uy || lz>uz) {
      throw runtime_error("FieldOfView:set_boundingbox: illegal lower/upper bounds");}
    lb_.clear();
    lb_.push_back(lt);
    lb_.push_back(lx);
    lb_.push_back(ly);
    lb_.push_back(lz);

    ub_.clear();
    ub_.push_back(ut);
    ub_.push_back(ux);
    ub_.push_back(uy);
    ub_.push_back(uz);

    return *this;
  }

  bool FieldOfView::contains( double t, double x, double y, double z ) const {
    if(    lb_[0] <= t && t <= ub_[0]
	&& lb_[1] <= x && x <= ub_[1]
	&& lb_[2] <= y && y <= ub_[2]
	&& lb_[3] <= z && z <= ub_[3]
	   ) {
      return true;
    } else {
      return false;
    }
  }

  const std::vector<double>& FieldOfView::lower_bound() const {
    return lb_;
  }

  const std::vector<double>& FieldOfView::upper_bound() const {
    return ub_;
  }

  namespace {
    void print_vec(const double v[]) {
      cout << "<" << v[0] << " " << v[1] << " " << v[2] << ">" << "\n";
    }

    // b-a
    void diff(const double a[], const double b[], double res[]) {
      res[0] = b[0] - a[0];
      res[1] = b[1] - a[1];
      res[2] = b[2] - a[2];
    }

    double dot(const double v1[], const double v2[]) {
      return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];      
    }

    double norm(const double v[]) {
      return sqrt(dot(v, v));
    }

    void cross(const double a[], const double b[], double res[]) {
      res[0] = a[1]*b[2] - a[2]*b[1];
      res[1] = a[2]*b[0] - a[0]*b[2];
      res[2] = a[0]*b[1] - a[1]*b[0];
    }

    void hesse_normal(const double v1[], const double v2[], double res[]) {
      double temp[3];
      cross(v1, v2, temp);
      double n = norm(temp);
      assert(n != 0);
      res[0] = temp[0] / n;
      res[1] = temp[1] / n;
      res[2] = temp[2] / n;
    }
    
    double abs_distance(const double p1[], const double p2[], const double p3[], const double q[]) {
      double u[3], v[3];
      diff(p1, p2, u);
      diff(p1, p3, v);

      double normal[3];
      hesse_normal(u, v, normal);

      double w[3];
      diff(p1, q, w);
      
      return abs(dot(w,normal));
    }
  }

  double FieldOfView::spatial_margin( double /*t*/, double x, double y, double z ) const {
    // the eight corners of the fov cube
    double c1[3], c2[3], c3[3], c4[3], c5[3], c6[3], /*c7[3],*/ c8[3];
    c1[0] = lb_[0]; c1[1] = lb_[1]; c1[2] = lb_[2];
    c2[0] = ub_[0]; c2[1] = lb_[1]; c2[2] = lb_[2];
    c3[0] = ub_[0]; c3[1] = ub_[1]; c3[2] = lb_[2];
    c4[0] = lb_[0]; c4[1] = ub_[1]; c4[2] = lb_[2];
    
    c5[0] = lb_[0]; c5[1] = lb_[1]; c5[2] = ub_[2];
    c6[0] = ub_[0]; c6[1] = lb_[1]; c6[2] = ub_[2];
    //c7[0] = ub_[0]; c7[1] = ub_[1]; c7[2] = ub_[2]; // unused
    c8[0] = lb_[0]; c8[1] = ub_[1]; c8[2] = ub_[2];

    // distances to the six faces of the cube
    double ds[6];
    double q[3] = {x, y, z};
    ds[0] = abs_distance(c1, c2, c5, q);
    ds[1] = abs_distance(c2, c3, c6, q);
    ds[2] = abs_distance(c4, c3, c8, q);
    ds[3] = abs_distance(c1, c4, c5, q);
    ds[4] = abs_distance(c1, c2, c4, q);
    ds[5] = abs_distance(c5, c6, c8, q);
    cout << ds[0] << " " << ds[1] << "\n";
    return *min_element(ds, ds+6);
  }
  
  double FieldOfView::temporal_margin( double t, double /*x*/, double /*y*/, double /*z*/ ) const {
    return min(abs(t - lb_[0]), abs(t - ub_[0]));
  }
  
}
