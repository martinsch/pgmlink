#define BOOST_TEST_MODULE learning_systemtest

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/archive/text_oarchive.hpp> // has to be include even though we dont use oarchives here
#include <boost/archive/text_iarchive.hpp>
#include <boost/shared_ptr.hpp>
#include <lemon/maps.h>

#include "pgmlink/event.h"
#include "pgmlink/feature.h"
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

  { // print some stats
  SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(6, 50);
  SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
  shared_ptr<HypothesesGraph> graph = shared_ptr<HypothesesGraph>(hyp_builder.build()); 
  pgm::chaingraph::TrainableModelBuilder b;
  b.without_detection_vars().without_divisions()
    .appearance(ConstantFeature(500))
    .disappearance(ConstantFeature(500))
    .move(SquaredDistance())
    ;
  //b.with_detection_vars().with_divisions();
  Chaingraph c(b, true);
  c.formulate(*graph);
  c.infer();
  c.conclude(*graph);
  prune_inactive(*graph);

  boost::shared_ptr<std::vector< std::vector<Event> > > es = events(*graph);
  EventsStatistics total;
  cout << "Tracking results:\n";
  for(size_t t=0; t < es->size(); ++t) {
    EventsStatistics stats = collect_events_statistics(es->operator[](t).begin(), es->operator[](t).end());
    total += stats;
    cout << "At t="<<t<<": " << stats  << "\n";
  }
  cout << "Total: " << total << "\n";
  }

  SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(6, 50);
  SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
  shared_ptr<HypothesesGraph> graph = shared_ptr<HypothesesGraph>(hyp_builder.build()); 
  pgm::chaingraph::TrainableModelBuilder b;
  b.without_detection_vars().without_divisions()
    .appearance(ConstantFeature(500))
    .disappearance(ConstantFeature(500))
    .move(SquaredDistance())
    ;

  //b.with_detection_vars().with_divisions();
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
