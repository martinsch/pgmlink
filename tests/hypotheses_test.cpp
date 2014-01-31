#define BOOST_TEST_MODULE hypotheses_test

#include <vector>
#include <string>
#include <iostream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <lemon/core.h>
#include <lemon/concepts/digraph.h>
#include <lemon/list_graph.h>
#include <lemon/maps.h>

#include "pgmlink/hypotheses.h"
#include "pgmlink/traxels.h"

using namespace pgmlink;
using namespace std;
using namespace boost;

BOOST_AUTO_TEST_CASE( HypothesesGraph_add_node ) {
    HypothesesGraph graph;
    graph.add_node(13);
}

BOOST_AUTO_TEST_CASE( HypothesesGraph_serialize ) {
  HypothesesGraph g;
  HypothesesGraph::Node n00 = g.add_node(0);
  HypothesesGraph::Node n01 = g.add_node(1);
  HypothesesGraph::Node n11 = g.add_node(1);
  /*HypothesesGraph::Arc a0 =*/ g.addArc(n00, n01);
  /*HypothesesGraph::Arc a1 =*/ g.addArc(n00, n11);

  Traxel tr00, tr01, tr11;
  feature_array com00(3);
  feature_array com01(3);
  feature_array com11(3);

  com00[0] = 0;
  com00[1] = 0;
  com00[2] = 0;
  tr00.features["com"] = com00;
  tr00.Id = 5;
  tr00.Timestep = 0;

  com01[0] = 0;
  com01[1] = 1;
  com01[2] = 0;
  tr01.features["com"] = com01;
  tr01.Id = 7;
  tr01.Timestep = 1;

  com11[0] = 1;
  com11[1] = 0;
  com11[2] = 0;
  tr11.features["com"] = com11;
  tr11.Id = 9;
  tr11.Timestep = 1;

  g.add(node_traxel());
  g.get(node_traxel()).set(n00, tr00);
  g.get(node_traxel()).set(n01, tr01);
  g.get(node_traxel()).set(n11, tr11);

  //////////

  // save to string
  string s;
  {
    stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa & g;
    s = ss.str();
  }
  
  // load from string and compare
  HypothesesGraph loaded;
  {
    stringstream ss(s);
    boost::archive::text_iarchive ia(ss);
    ia & loaded;
  }

  loaded.get(node_traxel()); // shouldn't throw an exception

  BOOST_CHECK_EQUAL_COLLECTIONS(g.timesteps().begin(),
                                g.timesteps().end(),
                                loaded.timesteps().begin(),
                                loaded.timesteps().end());
  BOOST_CHECK_EQUAL(countNodes(loaded), countNodes(g));
  BOOST_CHECK_EQUAL(countArcs(loaded), countArcs(g));
  
}

BOOST_AUTO_TEST_CASE( lgf_serialization ) {
  HypothesesGraph g;
  HypothesesGraph::Node n00 = g.add_node(0);
  HypothesesGraph::Node n01 = g.add_node(1);
  HypothesesGraph::Node n11 = g.add_node(1);
  /*HypothesesGraph::Arc a0 =*/ g.addArc(n00, n01);
  /*HypothesesGraph::Arc a1 =*/ g.addArc(n00, n11);

  cout << "without traxels\n";
  write_lgf(g);
  cout << "\n";

  // now with some optional node/arc maps
  Traxel tr00, tr01, tr11;
  feature_array com00(3);
  feature_array com01(3);
  feature_array com11(3);

  com00[0] = 0;
  com00[1] = 0;
  com00[2] = 0;
  tr00.features["com"] = com00;
  tr00.Id = 5;
  tr00.Timestep = 0;

  com01[0] = 0;
  com01[1] = 1;
  com01[2] = 0;
  tr01.features["com"] = com01;
  tr01.Id = 7;
  tr01.Timestep = 1;

  com11[0] = 1;
  com11[1] = 0;
  com11[2] = 0;
  tr11.features["com"] = com11;
  tr11.Id = 9;
  tr11.Timestep = 1;

  g.add(node_traxel());
  g.get(node_traxel()).set(n00, tr00);
  g.get(node_traxel()).set(n01, tr01);
  g.get(node_traxel()).set(n11, tr11);

  stringstream ss;
  write_lgf(g, ss, true);
  cout << "with traxels\n";
  cout << ss.str();
  cout << "\n";

  HypothesesGraph g_new;
  read_lgf(g_new, ss, true);
  cout << "restored from lgf\n";
  stringstream ss2;
  write_lgf(g_new, ss2, true);
  cout << ss2.str();
  cout << "\n";

  BOOST_CHECK_EQUAL(ss.str(), ss2.str());
}

BOOST_AUTO_TEST_CASE( SingleTimestepTraxel_HypothesesBuilder_build ) {
    Traxel tr11, tr12, tr21, tr22;
    feature_array com11(3);
    feature_array com12(3);
    feature_array com21(3);
    feature_array com22(3);

    com11[0] = 0;
    com11[1] = 0;
    com11[2] = 0;
    tr11.features["com"] = com11;
    tr11.Id = 5;
    tr11.Timestep = 0;

    com12[0] = 0;
    com12[1] = 1;
    com12[2] = 0;
    tr12.features["com"] = com12;
    tr12.Id = 7;
    tr12.Timestep = 0;

    com21[0] = 1;
    com21[1] = 0;
    com21[2] = 0;
    tr21.features["com"] = com21;
    tr21.Id = 9;
    tr21.Timestep = 1;

    com22[0] = 1;
    com22[1] = 1;
    com22[2] = 0;
    tr22.features["com"] = com22;
    tr22.Id = 11;
    tr22.Timestep = 1;

    TraxelStore ts;
    add(ts, tr11);
    add(ts, tr12);
    add(ts, tr21);
    add(ts, tr22);

    SingleTimestepTraxel_HypothesesBuilder builder(&ts);

    //builder.build();
}

// EOF
