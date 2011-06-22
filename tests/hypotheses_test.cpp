#define BOOST_TEST_MODULE hypotheses_test

#include <vector>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <lemon/core.h>
#include <lemon/concepts/digraph.h>
#include <lemon/list_graph.h>
#include <lemon/maps.h>

#include "hypotheses.h"
#include "traxels.h"

using namespace Tracking;
using namespace std;
using namespace boost;

BOOST_AUTO_TEST_CASE( HypothesesGraph_add_node ) {
    HypothesesGraph graph;
    graph.add_node(13);
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
