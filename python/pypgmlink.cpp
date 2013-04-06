#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray

#include <boost/python.hpp>
#include <vigra/numpy_array.hxx>
#include <iostream>
#include <vector>
#include <complex>


//forward declarations
void export_field_of_view();
void export_hypotheses();
void export_track();
void export_traxels();
void export_cross_correlation();
void export_gmm();


BOOST_PYTHON_MODULE( pgmlink )
{
	vigra::import_vigranumpy();
    export_field_of_view();
    export_hypotheses();
    export_track();
    export_traxels();
    export_cross_correlation();
    export_gmm();
}
