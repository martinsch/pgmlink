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
  std::map<std::string, bool> config; // default value for entries == false!
  write_lgf(g, std::cout, config);
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
  config["node_traxel"] = true;
  write_lgf(g, ss, config);
  cout << "with traxels\n";
  cout << ss.str();
  cout << "\n";

  HypothesesGraph g_new;
  read_lgf(g_new, ss, config);
  cout << "restored from lgf\n";
  stringstream ss2;
  write_lgf(g_new, ss2, config);
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

    builder.build();
}


BOOST_AUTO_TEST_CASE( SingleTimestepTraxel_HypothesesBuilder_build_forward_backward ) {
    Traxel tr11, tr12, tr21, tr22, tr23;
    feature_array com11(3);
    feature_array com12(3);
    feature_array com21(3);
    feature_array com22(3);
    feature_array com23(3);

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

    com21[0] = 0;
    com21[1] = 0;
    com21[2] = 0;
    tr21.features["com"] = com21;
    tr21.Id = 9;
    tr21.Timestep = 1;

    com22[0] = 0;
    com22[1] = 1;
    com22[2] = 0;
    tr22.features["com"] = com22;
    tr22.Id = 11;
    tr22.Timestep = 1;

    com22[0] = 0;
	com23[1] = 2;
	com23[2] = 0;
	tr23.features["com"] = com23;
	tr23.Id = 8;
	tr23.Timestep = 1;

    TraxelStore ts;
    add(ts, tr11);
    add(ts, tr12);
    add(ts, tr21);
    add(ts, tr22);
    add(ts, tr23);

    SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(1, // max_nn
    		10, // max_distance
    		true, // forward_backward
    		false, // consider_divisions
    		 0.5 //division_threshold
    		);
    SingleTimestepTraxel_HypothesesBuilder builder(&ts, builder_opts);
    HypothesesGraph& g = *(builder.build());

    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    size_t number_of_arcs = 0;
    for(HypothesesGraph::ArcIt a(g); a!=lemon::INVALID; ++a) {
    		HypothesesGraph::Node from = g.source(a);
    		HypothesesGraph::Node to = g.target(a);
    		Traxel from_tr = traxel_map[from];
    		Traxel to_tr = traxel_map[to];
    		switch(from_tr.Id) {
    		case 5:
    			BOOST_CHECK_EQUAL(to_tr.Id,9);
    			break;
    		case 7:
    			BOOST_CHECK(to_tr.Id == 8 || to_tr.Id == 11);
    			break;
    		default:
				BOOST_CHECK(false);
				break;
    		}
    		++number_of_arcs;
    }
    BOOST_CHECK_EQUAL(number_of_arcs,3);

}


BOOST_AUTO_TEST_CASE( SingleTimestepTraxel_HypothesesBuilder_build_divisions ) {
    Traxel tr11, tr12, tr21, tr22, tr23;
    feature_array com11(3);
    feature_array com12(3);
    feature_array com21(3);
    feature_array com22(3);
    feature_array com23(3);
    feature_array divProbNo(1);
    feature_array divProbYes(1);

    divProbNo[0] = 0.1;
    divProbYes[0] = 0.6;


    com11[0] = 0;
    com11[1] = 0;
    com11[2] = 0;
    tr11.features["com"] = com11;
    tr11.features["divProb"] = divProbNo;
    tr11.Id = 5;
    tr11.Timestep = 0;

    com12[0] = 0;
    com12[1] = 1;
    com12[2] = 0;
    tr12.features["com"] = com12;
    tr12.features["divProb"] = divProbYes;
    tr12.Id = 7;
    tr12.Timestep = 0;

    com21[0] = 0;
    com21[1] = 0;
    com21[2] = 0;
    tr21.features["com"] = com21;
    tr21.features["divProb"] = divProbNo;
    tr21.Id = 9;
    tr21.Timestep = 1;

    com22[0] = 0;
    com22[1] = 1;
    com22[2] = 0;
    tr22.features["com"] = com22;
    tr22.features["divProb"] = divProbNo;
    tr22.Id = 11;
    tr22.Timestep = 1;

    com22[0] = 0;
	com23[1] = 2;
	com23[2] = 0;
	tr23.features["com"] = com23;
	tr23.features["divProb"] = divProbNo;
	tr23.Id = 8;
	tr23.Timestep = 1;

    TraxelStore ts;
    add(ts, tr11);
    add(ts, tr12);
    add(ts, tr21);
    add(ts, tr22);
    add(ts, tr23);

    SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(1, // max_nn
    		10, // max_distance
    		false, // forward_backward
    		true, // consider_divisions
    		 0.5 //division_threshold
    		);
    SingleTimestepTraxel_HypothesesBuilder builder(&ts, builder_opts);
    HypothesesGraph& g = *(builder.build());

    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    size_t number_of_arcs = 0;
    for(HypothesesGraph::ArcIt a(g); a!=lemon::INVALID; ++a) {
    		HypothesesGraph::Node from = g.source(a);
    		HypothesesGraph::Node to = g.target(a);
    		Traxel from_tr = traxel_map[from];
    		Traxel to_tr = traxel_map[to];
    		switch(from_tr.Id) {
    		case 5:
    			BOOST_CHECK_EQUAL(to_tr.Id, 9);
    			break;
    		case 7:
    			BOOST_CHECK(to_tr.Id == 8 || to_tr.Id == 11);
    			break;
    		default:
				BOOST_CHECK(false);
				break;
    		}
    		++number_of_arcs;
    }
    BOOST_CHECK_EQUAL(number_of_arcs,3);
}



BOOST_AUTO_TEST_CASE( SingleTimestepTraxel_HypothesesGraph_generateTraxelGraph ) {
	HypothesesGraph traxel_graph;
	HypothesesGraph tracklet_graph;

	typedef HypothesesGraph::Node Node;
	typedef HypothesesGraph::Arc Arc;

	std::cout << "Adding nodes and arcs to traxel graph" << std::endl;
	Node n11 = traxel_graph.add_node(1);
	Node n12 = traxel_graph.add_node(1);
	Node n13 = traxel_graph.add_node(1);
	traxel_graph.add_node(1); // n14
	Node n21 = traxel_graph.add_node(2);
	Node n22 = traxel_graph.add_node(2);
	Node n23 = traxel_graph.add_node(2);
	Node n31 = traxel_graph.add_node(3);
	Node n32 = traxel_graph.add_node(3);
	Node n33 = traxel_graph.add_node(3);
	Node n41 = traxel_graph.add_node(4);
	Node n42 = traxel_graph.add_node(4);
	Node n43 = traxel_graph.add_node(4);

	traxel_graph.addArc(n11,n21);
	traxel_graph.addArc(n12,n21);
	traxel_graph.addArc(n12,n22);
	Arc a1_33 = traxel_graph.addArc(n13,n23);
	Arc a2_11 = traxel_graph.addArc(n21,n31);
	Arc a2_22 = traxel_graph.addArc(n22,n32);
	traxel_graph.addArc(n23,n33);
	Arc a3_11 = traxel_graph.addArc(n31,n41);
	traxel_graph.addArc(n32,n42);
	traxel_graph.addArc(n32,n43);
	traxel_graph.addArc(n33,n43);

//	Arc a1_11 = traxel_graph.addArc(n11,n21);
//	Arc a1_21 = traxel_graph.addArc(n12,n21);
//	Arc a1_22 = traxel_graph.addArc(n12,n22);
//	Arc a1_33 = traxel_graph.addArc(n13,n23);
//	Arc a2_11 = traxel_graph.addArc(n21,n31);
//	Arc a2_22 = traxel_graph.addArc(n22,n32);
//	Arc a2_33 = traxel_graph.addArc(n23,n33);
//	Arc a3_11 = traxel_graph.addArc(n31,n41);
//	Arc a3_22 = traxel_graph.addArc(n32,n42);
//	Arc a3_23 = traxel_graph.addArc(n32,n43);
//	Arc a3_33 = traxel_graph.addArc(n33,n43);

	std::cout << "Setting some arcs of traxel graph active" << std::endl;
	// set some arcs active
	traxel_graph.add(arc_active());
	property_map<arc_active, HypothesesGraph::base_graph > ::type& active_arcs = traxel_graph.get(arc_active());
	for(HypothesesGraph::ArcIt a(traxel_graph); a != lemon::INVALID; ++a) {
		active_arcs.set(a, false);
	}
	active_arcs.set(a1_33, true);
	active_arcs.set(a2_11, true);
	active_arcs.set(a2_22, true);
	active_arcs.set(a3_11, true);

	std::cout << "Adding traxels to the nodes of the traxel graph" << std::endl;
	// add traxels to the graph nodes
	traxel_graph.add(node_traxel());
	property_map<node_timestep, HypothesesGraph::base_graph>::type& node_timestep_map = traxel_graph.get(node_timestep());
	property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = traxel_graph.get(node_traxel());
	size_t id = 0;
	for(HypothesesGraph::NodeIt n(traxel_graph); n != lemon::INVALID; ++n) {
		Traxel traxel;
		traxel.Id = ++id;
		traxel.Timestep = node_timestep_map[n];
		traxel_map.set(n, traxel);
	}

	std::cout << "Generating tracklet graph" << std::endl;
	generateTrackletGraph(traxel_graph, tracklet_graph);

	std::cout << "Checking result" << std::endl;
//	tracklet_graph.add(node_tracklet());
	property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = tracklet_graph.get(node_tracklet());
	size_t num_of_nodes = 0;
	size_t num_of_arcs = 0;
	size_t num_of_tracklets_size[] = {0,0,0};

	for(HypothesesGraph::NodeIt n(tracklet_graph); n!=lemon::INVALID; ++n) {
		++num_of_tracklets_size[tracklet_map[n].size()-1];
		++num_of_nodes;
	}
	for(HypothesesGraph::ArcIt a(tracklet_graph); a!=lemon::INVALID; ++a) {
		++num_of_arcs;
	}
	BOOST_CHECK_EQUAL(num_of_nodes,9);
	BOOST_CHECK_EQUAL(num_of_tracklets_size[0],6); // tracklets of length 1
	BOOST_CHECK_EQUAL(num_of_tracklets_size[1],2); // tracklets of length 2
	BOOST_CHECK_EQUAL(num_of_tracklets_size[2],1); // tracklets of length 3
	BOOST_CHECK_EQUAL(num_of_arcs,7);


	std::cout << "Generating tracklet graph from tracklet graph" << std::endl;
	// make a tracklet graph from the tracklet graph -> should be the same!
	HypothesesGraph tracklet_tracklet_graph;
	tracklet_graph.add(arc_active());
	generateTrackletGraph(tracklet_graph, tracklet_tracklet_graph);
//	tracklet_tracklet_graph.add(node_track)
	property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map2 = tracklet_tracklet_graph.get(node_tracklet());
	num_of_nodes = 0;
	num_of_arcs = 0;
	num_of_tracklets_size[0] = 0; num_of_tracklets_size[1] = 0; num_of_tracklets_size[2] = 0;

	for(HypothesesGraph::NodeIt n(tracklet_tracklet_graph); n!=lemon::INVALID; ++n) {
		++num_of_tracklets_size[tracklet_map2[n].size()-1];
		++num_of_nodes;
	}
	for(HypothesesGraph::ArcIt a(tracklet_tracklet_graph); a!=lemon::INVALID; ++a) {
		++num_of_arcs;
	}
	BOOST_CHECK_EQUAL(num_of_nodes,9);
	BOOST_CHECK_EQUAL(num_of_tracklets_size[0],6); // tracklets of length 1
	BOOST_CHECK_EQUAL(num_of_tracklets_size[1],2); // tracklets of length 2
	BOOST_CHECK_EQUAL(num_of_tracklets_size[2],1); // tracklets of length 3
	BOOST_CHECK_EQUAL(num_of_arcs,7);
}

BOOST_AUTO_TEST_CASE( HypothesesGraph_Map_Serialization )
{
    // build example and run conservation tracking

    //  t=1      2      3       4
    //  o                       o
    //    |                    |
    //      ---- o ---- o ----
    //    |                    |
    //  o                       o
    TraxelStore ts;
    Traxel n11, n12, n21, n31, n41, n42;
    feature_array com(feature_array::difference_type(3));
    feature_array divProb(feature_array::difference_type(1));
    n11.Id = 1; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
    n11.features["com"] = com; n11.features["divProb"] = divProb;
    add(ts,n11);
    n12.Id = 3; n12.Timestep = 1; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1;
    n12.features["com"] = com; n12.features["divProb"] = divProb;
    add(ts,n12);
    n21.Id = 10; n21.Timestep = 2; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.1;
    n21.features["com"] = com; n21.features["divProb"] = divProb;
    add(ts,n21);
    n31.Id = 11; n31.Timestep = 3; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.1;
    n31.features["com"] = com; n31.features["divProb"] = divProb;
    add(ts,n31);
    n41.Id = 12; n41.Timestep = 4; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
    n41.features["com"] = com; n41.features["divProb"] = divProb;
    add(ts,n41);
    n42.Id = 13; n42.Timestep = 4; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
    n42.features["com"] = com; n42.features["divProb"] = divProb;
    add(ts,n42);

    SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(1, // max_nearest_neighbors
                20,
                true, // forward_backward
                true, // consider_divisions
                0.3
                );
    SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
    HypothesesGraph* graph = hyp_builder.build();

    {
        // set some property map values
        HypothesesGraph::NodeIt node(*graph); ++node;
        HypothesesGraph::ArcIt arc(*graph); ++arc;

        graph->add(arc_active_count());
        property_map<arc_active_count, HypothesesGraph::base_graph>::type& map_arc_active_cout = graph->get(arc_active_count());
        std::vector<bool> arc_active_vec; arc_active_vec.push_back(false);arc_active_vec.push_back(true);
        map_arc_active_cout.set(arc, arc_active_vec);

        graph->add(arc_active());
        property_map<arc_active, HypothesesGraph::base_graph>::type& map_arc_active = graph->get(arc_active());
        map_arc_active.set(arc, false);

        graph->add(node_active());
        property_map<node_active, HypothesesGraph::base_graph>::type& map_node_active = graph->get(node_active());
        map_node_active.set(node, true);

        graph->add(tracklet_intern_dist());
        property_map<tracklet_intern_dist, HypothesesGraph::base_graph>::type& map_node_tracklet_intern_dist = graph->get(tracklet_intern_dist());
        std::vector<double> tracklet_dist_vec; tracklet_dist_vec.push_back(7.32);
        map_node_tracklet_intern_dist.set(node, tracklet_dist_vec);
    }

    // serialize and deserialize
    HypothesesGraph graph_loaded;
    std::stringstream buf;
    {
        // now store the results
        boost::archive::text_oarchive out_archive(buf);
        out_archive << *graph;
    }

    {
        // read again
        boost::archive::text_iarchive in_archive(buf);
        in_archive >> graph_loaded;
    }

    // check which property maps have been read
    BOOST_CHECK(graph_loaded.has_property(arc_active_count()));
    BOOST_CHECK(graph_loaded.has_property(arc_active()));
    BOOST_CHECK(graph_loaded.has_property(node_active()));
    BOOST_CHECK(graph_loaded.has_property(tracklet_intern_dist()));
    BOOST_CHECK(graph_loaded.has_property(division_active()) == false);

    // get property maps
    property_map<arc_active_count, HypothesesGraph::base_graph>::type& map_arc_active_cout = graph_loaded.get(arc_active_count());
    property_map<node_active, HypothesesGraph::base_graph>::type& map_node_active = graph_loaded.get(node_active());
    property_map<tracklet_intern_dist, HypothesesGraph::base_graph>::type& map_node_tracklet_intern_dist = graph_loaded.get(tracklet_intern_dist());
    property_map<arc_active, HypothesesGraph::base_graph>::type& map_arc_active = graph_loaded.get(arc_active());

    // check whether the properties have been read again
    size_t arc_id = 0;
    for(HypothesesGraph::ArcIt a(graph_loaded); a!=lemon::INVALID; ++a)
    {
        // check transition map stuff
        if(arc_id == 1)
        {
            BOOST_CHECK(map_arc_active_cout[a].size() == 2);
            BOOST_CHECK(map_arc_active_cout[a][0] == false);
            BOOST_CHECK(map_arc_active_cout[a][1] == true);

            BOOST_CHECK(map_arc_active[a] == false);
        }

        arc_id++;
    }

    size_t node_id = 0;
    for(HypothesesGraph::NodeIt n(graph_loaded); n!=lemon::INVALID; ++n)
    {
        // check transition map stuff
        if(node_id == 1)
        {
            BOOST_CHECK(map_node_active[n] == true);

            BOOST_CHECK(map_node_tracklet_intern_dist[n].size() == 1);
            BOOST_CHECK(fabs(map_node_tracklet_intern_dist[n][0] - 7.32) < 0.000001);
        }

        node_id++;
    }
}
