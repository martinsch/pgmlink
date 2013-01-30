#define BOOST_TEST_MODULE learning_systemtest

#include <vector>
#include <iostream>
#include <fstream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/archive/text_oarchive.hpp> // has to be include even though we dont use oarchives here
#include <boost/archive/text_iarchive.hpp>
#include <boost/shared_ptr.hpp>
#include <lemon/maps.h>

#include "pgmlink/reasoner_pgm.h"
#include "pgmlink/pgm_chaingraph.h"
#include "pgmlink/traxels.h"
#include "pgmlink/hypotheses.h"

using namespace std;
using namespace boost;
using namespace pgmlink;

BOOST_AUTO_TEST_CASE( learning_from_autolabels ) {
  ifstream ifs("@PROJECT_SOURCE_DIR@/tests/system/traxelstore_dros_0-5.boost");
  pgmlink::TraxelStore ts;
  boost::archive::text_iarchive ia(ifs);
  ia >> ts;

  SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(6, 50);
  SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
  shared_ptr<HypothesesGraph> graph = shared_ptr<HypothesesGraph>(hyp_builder.build()); 
  pgm::chaingraph::ECCV12ModelBuilder b;
  b.with_detection_vars().with_divisions();
  Chaingraph c(b, true);
  c.formulate(*graph);
  c.infer();
  c.conclude(*graph);

  property_map<node_active, HypothesesGraph::base_graph>::type node_labels(*graph);
  lemon::mapCopy(*graph, graph->get(node_active()), node_labels);
  property_map<arc_active, HypothesesGraph::base_graph>::type arc_labels(*graph);
  lemon::mapCopy(*graph, graph->get(arc_active()), arc_labels);
  vector<pgm::OpengmModel::ValueType> weights;

  pgm::chaingraph::ModelTrainer trainer;
  weights = trainer.train(graph.get(), graph.get()+1, &node_labels, &arc_labels);
}



// Local Variables:
// mode: c++
// End:
