#include <vector>

#include "../include/pgmlink/tracking.h"
#include <boost/utility.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python.hpp>
#include <Python.h>


using namespace std;
using namespace pgmlink;
using namespace boost::python;

vector<vector<Event> > pythonChaingraphTracking(ChaingraphTracking& tr, TraxelStore& ts) {
	vector<vector<Event> > result = std::vector<std::vector<Event> >(0);
	// release the GIL
	Py_BEGIN_ALLOW_THREADS
	try {
		result = tr(ts);
	} catch (std::exception& e) {
		Py_BLOCK_THREADS
		throw;
	}
	Py_END_ALLOW_THREADS
	return result;
}

void export_track() {
    class_<vector<Event> >("EventVector")
	.def(vector_indexing_suite<vector<Event> >())
    ;

    class_<vector<vector<Event> > >("NestedEventVector")
	.def(vector_indexing_suite<vector<vector<Event> > >())
    ;

    class_<map<unsigned int, bool> >("DetectionMap")
      .def(map_indexing_suite<map<unsigned int, bool> >())
    ;

    class_<vector<map<unsigned int, bool> > >("DetectionMapsVector")
      .def(vector_indexing_suite<vector<map<unsigned int, bool> > >())
    ;

    class_<ChaingraphTracking>("ChaingraphTracking", 
			       init<string,double,double,double,double,
			       	   bool,double,double,bool,
			       	   bool,double,double,double,double>(
							  args("random_forest_filename", "appearance", "disappearance", "detection", "misdetection",
									  "use_random_forest", "opportunity_cost", "forbidden_cost", "with_constraints",
									  "fixed_detections", "mean_div_dist", "min_angle", "ep_gap", "n_neighbors"
									  )))
      .def("__call__", &pythonChaingraphTracking)
      .def("detections", &ChaingraphTracking::detections)
      .def("set_with_divisions", &ChaingraphTracking::set_with_divisions)
      .def("set_cplex_timeout", &ChaingraphTracking::set_cplex_timeout)
    ;

    enum_<Event::EventType>("EventType")
	.value("Move", Event::Move)
	.value("Division", Event::Division)
	.value("Appearance", Event::Appearance)
	.value("Disappearance", Event::Disappearance)
	.value("Void", Event::Void)
    ;

    class_<vector<unsigned int> >("IdVector")
	.def(vector_indexing_suite<vector<unsigned int> >())
    ;
    class_<Event>("Event")
	.def_readonly("type", &Event::type)
	.def_readonly("traxel_ids", &Event::traxel_ids)
	.add_property("energy", &Event::energy)
    ;
}
