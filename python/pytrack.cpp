#define BOOST_PYTHON_MAX_ARITY 20

#include <vector>

#include <boost/python.hpp>
#include <boost/utility.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "../include/pgmlink/tracking.h"

using namespace std;
using namespace pgmlink;
using namespace boost::python;

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
			       init<string,double,double,double,double,bool,double,double,bool,bool,double,double,double,int>(
						        		    args("random_forest_filename", "appearance", "disappearance", "detection", "misdetection", "use_random_forest", "opportunity_cost", "forbidden_cost", "with_constraints", "fixed_detections", "mean_div_dist", "min_angle", "ep_gap", "nneighbors")))
      .def("__call__", &ChaingraphTracking::operator())
      .def("detections", &ChaingraphTracking::detections) 
    ;

    class_<ConsTracking>("ConsTracking",
			init<int,double,double,string,bool,double,double,double,bool,double,double>(
						args("max_number_objects", "max_neighbor_distance", "division_threshold",
							"detection_rf_filename", "size_dependent_detection_prob", "forbidden_cost",
							"ep_gap", "avg_obj_size",
							"with_tracklets",
							"division_weight", "transition_weight")))
	  .def("__call__", &ConsTracking::operator())
	  .def("detections", &ConsTracking::detections)
	;

    class_<NNTracking>("NNTracking",
			init<double,double,std::vector<std::string>, double,bool,bool,std::vector<int> >(
				 args("divDist", "movDist", "features", "divisionThreshold", "splitterHandling", "mergerHandling", "maxTraxelIdAt")))
	  .def("__call__", &NNTracking::operator())
	  .def("detections", &NNTracking::detections)
        ;

    enum_<Event::EventType>("EventType")
	.value("Move", Event::Move)
	.value("Division", Event::Division)
	.value("Appearance", Event::Appearance)
	.value("Disappearance", Event::Disappearance)
	.value("Merger", Event::Merger)
	.value("Void", Event::Void)
    ;

    class_<vector<unsigned int> >("IdVector")
	.def(vector_indexing_suite<vector<unsigned int> >())
    ;
    class_<Event>("Event")
	.def_readonly("type", &Event::type)
	.def_readonly("traxel_ids", &Event::traxel_ids)
	.def_readonly("energy", &Event::energy)
    ;
}
