#include <boost/python.hpp>
#include <boost/utility.hpp>

#include "../include/pgmlink/field_of_view.h"

using namespace pgmlink;
using namespace boost::python;

void export_field_of_view() {
    class_< FieldOfView >("FieldOfView")
      .def(init<double, double, double, double, double, double, double, double>(
	   args("lt","lx","ly","lz","ut","ux","uy","uz")))
      .def("set_boundingbox", &FieldOfView::set_boundingbox, return_self<>())
      //.def("contains", &FieldOfView::contains)
      //.def("lower_bound", &FieldOfView::lower_bound, return_value_policy<copy_const_reference>())
      //.def("upper_bound", &FieldOfView::upper_bound, return_value_policy<copy_const_reference>())
      //.def("spatial_margin", &FieldOfView::spatial_margin)
      //.def("temporal_margin", &FieldOfView::temporal_margin)
    ;
}
