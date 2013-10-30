#define BOOST_TEST_MODULE multi_hypotheses_test

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

#include "pgmlink/multi_hypotheses_graph.h"
#include "pgmlink/traxels.h"

using namespace pgmlink;
using namespace std;
using namespace boost;

typedef MultiHypothesesGraph::Node Node;
typedef MultiHypothesesGraph::Arc Arc;

BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_serialize ) {
	MultiHypothesesGraph g;
	MultiHypothesesGraph::ContainedRegionsMap& regions = g.get(node_regions_in_component());

	std::cout << "Adding nodes and arcs" << std::endl;
	std::cout << std::endl;

	Node n1 = g.add_node(0);
	Node n2 = g.add_node(1);
	Node n3 = g.add_node(1);

	g.addArc(n1, n2);
	g.addArc(n1, n3);

	std::cout << "Adding contained regions to nodes" << std::endl;
	std::cout << std::endl;
	std::vector<Traxel>& t1 = regions.get_value(n1);
	std::vector<Traxel>& t2 = regions.get_value(n2);
	std::vector<Traxel>& t3 = regions.get_value(n3);

	feature_array com(3,1.);
	feature_array count(1, 0.1);


	t1.push_back(Traxel(1, 0));
	float conflict_arr11[] = {2., 3., 4., 5.};
	(t1.end()-1)->features["conflicts"] = feature_array(conflict_arr11, conflict_arr11 + 4);
	(t1.end()-1)->features["level"].push_back(0.);
	(t1.end()-1)->features["com"] = com;
	(t1.end()-1)->features["count_prediction"] = count;

	t1.push_back(Traxel(2, 0));
	float conflict_arr12[] = {1., 4., 5.};
	std::fill(com.begin(), com.end(), 2.);
	(t1.end()-1)->features["conflicts"] = feature_array(conflict_arr12, conflict_arr12 + 3);
	(t1.end()-1)->features["level"].push_back(1.);
	(t1.end()-1)->features["com"] = com;

	t1.push_back(Traxel(3, 0));
	float conflict_arr13[] = {1.};
	std::fill(com.begin(), com.end(), 0.);
	(t1.end()-1)->features["conflicts"] = feature_array(conflict_arr13, conflict_arr13 + 1);
	(t1.end()-1)->features["level"].push_back(1.);
	(t1.end()-1)->features["com"] = com;

	t1.push_back(Traxel(4, 0));
	float conflict_arr14[] = {1., 2.};
	std::fill(com.begin(), com.end(), 2.); com[0] = 1;
	(t1.end()-1)->features["conflicts"] = feature_array(conflict_arr14, conflict_arr14 + 2);
	(t1.end()-1)->features["level"].push_back(2.);
	(t1.end()-1)->features["com"] = com;

	t1.push_back(Traxel(5, 0));
	float conflict_arr15[] = {1., 2.};
	com[0] = 3;
	(t1.end()-1)->features["conflicts"] = feature_array(conflict_arr15, conflict_arr15 + 2);
	(t1.end()-1)->features["level"].push_back(2.);
	(t1.end()-1)->features["com"] = com;


	t2.push_back(Traxel(1, 1));
	float conflict_arr21[] = {3., 4.};
	std::fill(com.begin(), com.end(), 0.); com[0] = 1;
	(t2.end()-1)->features["conflicts"] = feature_array(conflict_arr21, conflict_arr21 + 2);
	(t2.end()-1)->features["level"].push_back(0.);
	(t2.end()-1)->features["com"] = com;
	(t2.end()-1)->features["count_prediction"] = count;

	t2.push_back(Traxel(3, 1));
	float conflict_arr23[] = {1.};
	com[0] = 0;
	(t2.end()-1)->features["conflicts"] = feature_array(conflict_arr23, conflict_arr23 + 1);
	(t2.end()-1)->features["level"].push_back(1.);
	(t2.end()-1)->features["com"] = com;

	t2.push_back(Traxel(4, 1));
	com[0] = 2;
	float conflict_arr24[] = {1.};
	(t2.end()-1)->features["conflicts"] = feature_array(conflict_arr24, conflict_arr24 + 1);
	(t2.end()-1)->features["level"].push_back(1.);
	(t2.end()-1)->features["com"] = com;


	t3.push_back(Traxel(2, 1));
	std::fill(com.begin(), com.end(), 2.);
	float conflict_arr22[] = {5., 6.};
	(t3.end()-1)->features["conflicts"] = feature_array(conflict_arr22, conflict_arr22 + 2);
	(t3.end()-1)->features["level"].push_back(0);
	(t3.end()-1)->features["com"] = com;
	(t3.end()-1)->features["count_prediction"] = count;

	t3.push_back(Traxel(5, 1));
	com[0] = 1;
	float conflict_arr25[] = {2.};
	(t3.end()-1)->features["conflicts"] = feature_array(conflict_arr25, conflict_arr25 + 1);
	(t3.end()-1)->features["level"].push_back(1);
	(t3.end()-1)->features["com"] = com;

	t3.push_back(Traxel(6, 1));
	com[0] = 3;
	float conflict_arr26[] = {2.};
	(t3.end()-1)->features["conflicts"] = feature_array(conflict_arr26, conflict_arr26 + 1);
	(t3.end()-1)->features["level"].push_back(1);
	(t3.end()-1)->features["com"] = com;

	//////////

	// save to string
	string s;
	{
		stringstream ss;
		boost::archive::text_oarchive oa(ss);
		oa & g;
		s = ss.str();
	}

	cout << s << endl;

	// load from string and compare
	MultiHypothesesGraph loaded;
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

	Traxel tr = loaded.get(node_traxel())[n1];
	cout << "deserialized traxel: " << tr << endl;
	FeatureMap f_loaded = tr.features;
	FeatureMap f_orig = t1[0].features;
	std::string name = "conflicts";
	BOOST_CHECK_EQUAL_COLLECTIONS(f_loaded[name].begin(), f_loaded[name].end(), f_orig[name].begin(), f_orig[name].end());
	name = "level";
	BOOST_CHECK_EQUAL_COLLECTIONS(f_loaded[name].begin(), f_loaded[name].end(), f_orig[name].begin(), f_orig[name].end());
	name = "com";
	BOOST_CHECK_EQUAL_COLLECTIONS(f_loaded[name].begin(), f_loaded[name].end(), f_orig[name].begin(), f_orig[name].end());
	name = "count_prediction";
	BOOST_CHECK_EQUAL_COLLECTIONS(f_loaded[name].begin(), f_loaded[name].end(), f_orig[name].begin(), f_orig[name].end());
  
}

//BOOST_AUTO_TEST_CASE( lgf_serialization ) {
//  HypothesesGraph g;
//  HypothesesGraph::Node n00 = g.add_node(0);
//  HypothesesGraph::Node n01 = g.add_node(1);
//  HypothesesGraph::Node n11 = g.add_node(1);
//  /*HypothesesGraph::Arc a0 =*/ g.addArc(n00, n01);
//  /*HypothesesGraph::Arc a1 =*/ g.addArc(n00, n11);
//
//  cout << "without traxels\n";
//  write_lgf(g);
//  cout << "\n";
//
//  // now with some optional node/arc maps
//  Traxel tr00, tr01, tr11;
//  feature_array com00(3);
//  feature_array com01(3);
//  feature_array com11(3);
//
//  com00[0] = 0;
//  com00[1] = 0;
//  com00[2] = 0;
//  tr00.features["com"] = com00;
//  tr00.Id = 5;
//  tr00.Timestep = 0;
//
//  com01[0] = 0;
//  com01[1] = 1;
//  com01[2] = 0;
//  tr01.features["com"] = com01;
//  tr01.Id = 7;
//  tr01.Timestep = 1;
//
//  com11[0] = 1;
//  com11[1] = 0;
//  com11[2] = 0;
//  tr11.features["com"] = com11;
//  tr11.Id = 9;
//  tr11.Timestep = 1;
//
//  g.add(node_traxel());
//  g.get(node_traxel()).set(n00, tr00);
//  g.get(node_traxel()).set(n01, tr01);
//  g.get(node_traxel()).set(n11, tr11);
//
//  stringstream ss;
//  write_lgf(g, ss, true);
//  cout << "with traxels\n";
//  cout << ss.str();
//  cout << "\n";
//
//  HypothesesGraph g_new;
//  read_lgf(g_new, ss, true);
//  cout << "restored from lgf\n";
//  stringstream ss2;
//  write_lgf(g_new, ss2, true);
//  cout << ss2.str();
//  cout << "\n";
//
//  BOOST_CHECK_EQUAL(ss.str(), ss2.str());
//}



// EOF
