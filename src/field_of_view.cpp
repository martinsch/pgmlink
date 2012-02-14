#include "field_of_view.h"
#include <stdexcept>

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
    lb.push_back(lt);
    lb.push_back(lx);
    lb.push_back(ly);
    lb.push_back(lz);

    ub_.clear();
    ub.push_back(ut);
    ub.push_back(ux);
    ub.push_back(uy);
    ub.push_back(uz);

    return *this;
  }

  boolean FieldOfView::contains( double t, double x, double y, double z ) {
    
    return false;
  }
}
