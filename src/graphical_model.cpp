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

#include "pgmlink/energy.h"
#include "pgmlink/graphical_model.h"
#include "pgmlink/event.h"
#include "pgmlink/traxels.h"

using namespace std;
using namespace boost;

namespace pgmlink {

  ////
  //// class OpengmModel
  ////
  OpengmModel::OpengmModel() {
    model_ = new ogmGraphicalModel();
  }
  
  OpengmModel::~OpengmModel() {
    delete model_;
  }

} /* namespace pgmlink */
