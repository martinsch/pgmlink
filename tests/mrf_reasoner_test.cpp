#define BOOST_TEST_MODULE mrf_reasoner_test

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/bind.hpp>

#include <lemon/color.h>
#include <lemon/graph_to_eps.h>

#include "graph.h"
#include "hypotheses.h"
#include "energy.h"
#include "reasoning/mrf_reasoner.h"
#include "traxels.h"

using namespace Tracking;
using namespace std;
using namespace boost;
using namespace lemon;

BOOST_AUTO_TEST_CASE( SingleTimestepTraxelMrf_construction )
{
    Traxels empty;
    ConstantEnergy e(25);
    
    SingleTimestepTraxelMrf mrf(bind<double>(e, _1, empty, empty),
			        bind<double>(e, _1, empty, empty),
			        bind<double>(e, _1, empty, empty),
			        bind<double>(e, _1, empty, empty),
				bind<double>(e, _1, _2, empty, empty),
				bind<double>(e, _1, _2, _3, empty, empty)
			       );
}



BOOST_AUTO_TEST_CASE( SingleTimestepTraxelMrf_workflow )
{
  cout << "-> workflow: scaffolding" << endl; 
  // scaffolding
    Traxel tr11, tr12, tr21, tr22;
    feature_array com11(feature_array::difference_type(3));
    feature_array com12(feature_array::difference_type(3));
    feature_array com21(feature_array::difference_type(3));
    feature_array com22(feature_array::difference_type(3));

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
    add(ts,tr11);
    add(ts,tr12);
    add(ts,tr21);
    add(ts,tr22);

    SingleTimestepTraxel_HypothesesBuilder builder(&ts);
    HypothesesGraph* graph = builder.build();

    Traxels empty;
    ConstantEnergy e(25);
    
    cout << "-> workflow: construction MRF reasoner" << endl; 
    SingleTimestepTraxelMrf mrf(bind<double>(e, _1, empty, empty),
			        bind<double>(e, _1, empty, empty),
			        bind<double>(e, _1, empty, empty),
			        bind<double>(e, _1, empty, empty),
				bind<double>(e, _1, _2, empty, empty),
				bind<double>(e, _1, _2, _3, empty, empty)
			       );

    // test
    cout << "-> workflow: formulating model" << endl; 
    mrf.formulate(*graph);
    cout << "-> workflow: infer" << endl; 
    mrf.infer();
    cout << "-> workflow: conclude" << endl; 
    mrf.conclude(*graph);
    prune_inactive(*graph);
}

BOOST_AUTO_TEST_CASE( diamond_pattern ) {
    /*
    //      --- tr21 ---
    //    /              \
    // tr11 --- tr22 --- tr31
    //    \              /
    //      --- tr23 ---
    */
    // trXY: traxel Y at timestep X
  Traxel tr11, tr21, tr22, tr23, tr31;
    feature_array com11(3);
    feature_array com21(3);
    feature_array com22(3);
    feature_array com23(3);
    feature_array com31(3);

    com11[0] = 0;
    com11[1] = 0;
    com11[2] = 0;
    tr11.features["com"] = com11;
    tr11.Id = 11;
    tr11.Timestep = 1;

    com21[0] = 1;
    com21[1] = 0;
    com21[2] = 0;
    tr21.features["com"] = com21;
    tr21.Id = 21;
    tr21.Timestep = 2;

    com22[0] = 1;
    com22[1] = 1;
    com22[2] = 0;
    tr22.features["com"] = com22;
    tr22.Id = 22;
    tr22.Timestep = 2;

    com23[0] = 1;
    com23[1] = 2;
    com23[2] = 0;
    tr23.features["com"] = com23;
    tr23.Id = 23;
    tr23.Timestep = 2;


    com31[0] = 2;
    com31[1] = 0;
    com31[2] = 0;
    tr31.features["com"] = com31;
    tr31.Id = 31;
    tr31.Timestep = 3;
 
    TraxelStore ts;
    add(ts,tr11);
    add(ts,tr21);
    add(ts,tr22);
    add(ts,tr23);
    add(ts,tr31);

    SingleTimestepTraxel_HypothesesBuilder builder(&ts);
    HypothesesGraph* graph = builder.build();

    typedef Tracking::property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map_t;
    node_traxel_map_t& node_traxel_map = graph->get(node_traxel());
    for(HypothesesGraph::NodeIt it(*graph); it!=lemon::INVALID; ++it) {
	cout << "Node id: " << graph->id(it) << endl;
	cout << "Traxel id: " << node_traxel_map[it].Id << endl;

	for(HypothesesGraph::OutArcIt ait(*graph, it); ait!=lemon::INVALID; ++ait) {
	    cout << "Arc id: " << graph->id(ait) << endl;
	    cout << "Target node id: " << graph->id(graph->target(ait)) << endl;
	}
	cout << endl;
    }


    Traxels empty;
    ConstantEnergy move(1);
    ConstantEnergy detection(1);
    ConstantEnergy non_detection(25);
    ConstantEnergy appearance(1);
    ConstantEnergy disappearance(1);
    ConstantEnergy division(25);
    
    cout << "-> diamond_pattern: construction MRF reasoner" << endl; 
    SingleTimestepTraxelMrf mrf(bind<double>(detection, _1, empty, empty),
			        bind<double>(non_detection, _1, empty, empty),
			        bind<double>(appearance, _1, empty, empty),
			        bind<double>(disappearance, _1, empty, empty),
				bind<double>(move, _1, _2, empty, empty),
				bind<double>(division, _1, _2, _3, empty, empty)
			       );

    // test
    cout << "-> diamond_pattern: formulating model" << endl; 
    mrf.formulate(*graph);
    cout << "-> diamond_pattern: infer" << endl; 
    mrf.infer();
    cout << "-> diamond_pattern: conclude" << endl; 
    mrf.conclude(*graph);

      //Tracking::property_map<node_active, HypothesesGraph::base_graph>::type& active_nodes = (*graph).get(node_active());
      Tracking::property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = (*graph).get(arc_active());
    
      // prune inactive arcs
      typedef Tracking::property_map<arc_active, HypothesesGraph::base_graph>::type::FalseIt inactive_arc_it;
      for(inactive_arc_it it(active_arcs); it!=lemon::INVALID; ++it) {
	cout << "Inactive arc: " << (*graph).id(it) << endl;
      } 

      typedef Tracking::property_map<arc_active, HypothesesGraph::base_graph>::type::TrueIt active_arc_it;
      for(active_arc_it it(active_arcs); it!=lemon::INVALID; ++it) {
	cout << "Active arc: " << (*graph).id(it) << endl;
      } 


    prune_inactive(*graph);
    events(*graph);
}

// EOF
