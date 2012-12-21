#include <boost/python.hpp>

//forward declarations
void export_field_of_view();
void export_hypotheses();
void export_track();
void export_traxels();

BOOST_PYTHON_MODULE( pgmlink )
{
    export_field_of_view();
    export_hypotheses();
    export_track();
    export_traxels();
}
