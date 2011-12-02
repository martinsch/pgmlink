#include <assert.h>
#include <ostream>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <utility>
#include <map>
#include <stdexcept>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graphviz.hpp>

#include "energy.h"
#include "graphical_model.h"
#include "event.h"
#include "traxels.h"

using namespace std;
using namespace boost;

namespace Tracking {

////
//// class OpengmMrf
////
    OpengmMrf::OpengmMrf() {
	model_ = new ogmGraphicalModel();
    }

    OpengmMrf::~OpengmMrf() {
	delete model_;
    }

} /* namespace Tracking */
