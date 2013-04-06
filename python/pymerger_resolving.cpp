#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

#include <boost/python.hpp>
#include "pgmlink/merger_resolving.h"


using namespace std;
using namespace pgmlink;
using namespace boost::python;


void export_gmm() {
  def("gmm_priors_and_centers", gmm_priors_and_centers);
}
