#define BOOST_TEST_MODULE learning_systemtest

#include <vector>
#include <iostream>
#include <fstream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/archive/text_oarchive.hpp> // has to be include even though we dont use oarchives here
#include <boost/archive/text_iarchive.hpp>

#include "pgmlink/traxels.h"
#include "pgmlink/hypotheses.h"

using namespace std;
using namespace boost;

BOOST_AUTO_TEST_CASE( learning_from_autolabels ) {
  ifstream ifs("@PROJECT_SOURCE_DIR@/tests/system/traxelstore_dros_0-5.boost");
  pgmlink::TraxelStore ts;
  boost::archive::text_iarchive ia(ifs);
  ia >> ts;

  SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(6, 50);
  SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
  
}


