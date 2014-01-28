#define BOOST_TEST_MODULE multi_hypotheses_test

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

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
#include "pgmlink/pgm_multi_hypotheses.h"
#include "pgmlink/reasoner_multi_hypotheses.h"

using namespace pgmlink;
using namespace std;
using namespace boost;

typedef MultiHypothesesGraph::Node Node;
typedef MultiHypothesesGraph::Arc Arc;
typedef pgm::multihypotheses::Model Model;

#if 0>1 // ignore this test for now

BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_serialize ) {

	MultiHypothesesTraxelStore ts;

	std::vector<Traxel>& t1 = ts.map[0][1];
	std::vector<Traxel>& t2 = ts.map[1][1];

	std::vector<std::vector<unsigned> >& c1 = ts.conflicts_by_timestep[0][1];
	std::vector<std::vector<unsigned> >& c2 = ts.conflicts_by_timestep[1][1];

	feature_array com(3, 1.);
	feature_array count(1, 0.1);
	feature_array var(1, 0.1);

	t1.push_back(Traxel(1, 0));
	float conflict_arr11[] = { 2., 3., 4., 5. };
	(t1.end() - 1)->features["conflicts"] = feature_array(conflict_arr11,
			conflict_arr11 + 4);
	(t1.end() - 1)->features["level"].push_back(0.);
	(t1.end() - 1)->features["com"] = com;
	(t1.end() - 1)->features["count_prediction"] = count;
	(t1.end() - 1)->features["Variance"] = var; // add a feature which will be removed while serialization

	c1.push_back(std::vector<unsigned>());
	c1.rbegin()->push_back(1);

	t2.push_back(Traxel(1, 1));
	float conflict_arr21[] = { 3., 4. };
	std::fill(com.begin(), com.end(), 0.);
	com[0] = 1;
	(t2.end() - 1)->features["conflicts"] = feature_array(conflict_arr21,
			conflict_arr21 + 2);
	(t2.end() - 1)->features["level"].push_back(0.);
	(t2.end() - 1)->features["com"] = com;
	(t2.end() - 1)->features["count_prediction"] = count;
	(t2.end() - 1)->features["Variance"] = count; // add a feature which will be removed while serialization

	t2.push_back(Traxel(3, 1));
	float conflict_arr23[] = { 1. };
	com[0] = 0;
	(t2.end() - 1)->features["conflicts"] = feature_array(conflict_arr23,
			conflict_arr23 + 1);
	(t2.end() - 1)->features["level"].push_back(1.);
	(t2.end() - 1)->features["com"] = com;

	t2.push_back(Traxel(4, 1));
	com[0] = 2;
	float conflict_arr24[] = { 1. };
	(t2.end() - 1)->features["conflicts"] = feature_array(conflict_arr24,
			conflict_arr24 + 1);
	(t2.end() - 1)->features["level"].push_back(1.);
	(t2.end() - 1)->features["com"] = com;

	c2.push_back(std::vector<unsigned>());
	c2.rbegin()->push_back(1);
	c2.rbegin()->push_back(3);
	c2.push_back(std::vector<unsigned>());
	c2.rbegin()->push_back(1);
	c2.rbegin()->push_back(4);

	std::cout << " -> workflow: building MultiHypothesesGraph" << std::endl;
	MultiHypothesesGraphBuilder::Options options(2, 1000000);
	MultiHypothesesGraphBuilder graph_builder(options);
	MultiHypothesesGraphPtr graph = graph_builder.build(ts);

	MultiHypothesesGraph& g = *graph;
	g.add(node_move_features());
	MultiHypothesesGraph::MoveFeatureMap& moves = g.get(node_move_features());
	const MultiHypothesesGraph::TraxelMap& traxels = g.get(node_traxel());
	MultiHypothesesGraph::Node n;
	for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
			it != traxels.endValue(); ++it) {
		if (*it == t1[0]) {
			n = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
			break;
		}
	}
	std::map<Traxel, feature_array>& move_probs = moves.get_value(n)[t1[0]];
	move_probs[t2[0]].push_back(0.3);
	move_probs[t2[0]].push_back(0.7);
	move_probs[t2[1]].push_back(0.99);
	move_probs[t2[1]].push_back(0.01);
	move_probs[t2[2]].push_back(0.05);
	move_probs[t2[2]].push_back(0.95);


	g.add(node_division_features());
	MultiHypothesesGraph::DivisionFeatureMap& divs = g.get(node_division_features());
//	for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
//			it != traxels.endValue(); ++it) {
//		if (*it == t1[0]) {
//			n = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
//			break;
//		}
//	}
	std::map<std::pair<Traxel, Traxel>, feature_array>& div_probs = divs.get_value(n)[t1[0]];
	std::vector<float> probs1, probs2;
	probs1.push_back(0.1);
	probs1.push_back(0.9);
	div_probs[std::make_pair(t2[0], t2[1])] = feature_array(probs1.size());
	std::copy(probs1.begin(),
			  probs1.end(),
			  div_probs[std::make_pair(t2[0], t2[1])].begin());
	probs2.push_back(0.1);
	probs2.push_back(0.9);
	div_probs[std::make_pair(t2[0], t2[2])] = feature_array(probs2.size());
	std::copy(probs2.begin(),
			  probs2.end(),
			  div_probs[std::make_pair(t2[0], t2[2])].begin());

//  // DEPRECATED
//	g.add(node_count_features());
//	MultiHypothesesGraph::CountFeatureMap& counts = g.get(node_count_features());
//	for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
//			it != traxels.endValue(); ++it) {
//		if (*it == t1[0]) {
//			n = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
//			break;
//		}
//	}
//	counts.get_value(n)[t1[0]] = 0.11;


//	g.add(node_conflict_sets());
//	MultiHypothesesGraph::ConflictSetMap& conflicts = g.get(node_conflict_sets());
//	for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
//			it != traxels.endValue(); ++it) {
//		if (*it == t1[0]) {
//			n = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
//			break;
//		}
//	}
//	MultiHypothesesGraph::ConflictSetMap::Value& conflict_sets = conflicts.get_value(n);
//	std::vector<unsigned> c_t1;
//	c_t1.push_back(5);
//	c_t1.push_back(2);
//	c_t1.push_back(1);
//	conflict_sets.insert(conflict_sets.end(), c_t1.begin(), c_t1.end());

	graph->remove_traxel_features();

	// save to string
	string s;
	{
		stringstream ss;
		boost::archive::text_oarchive oa(ss);
		oa & g;
		s = ss.str();
	}

	//cout << s << endl;

	// load from string and compare
	MultiHypothesesGraph loaded;
	{
		stringstream ss(s);
		boost::archive::text_iarchive ia(ss);
		ia & loaded;
	}

	loaded.get(node_traxel()); // shouldn't throw an exception
	loaded.get(node_conflict_sets()); // shouldn't throw an exception
	loaded.get(node_regions_in_component()); // shouldn't throw an exception

	BOOST_CHECK_EQUAL_COLLECTIONS(g.timesteps().begin(),
				g.timesteps().end(),
				loaded.timesteps().begin(),
				loaded.timesteps().end());

	BOOST_CHECK_EQUAL(countNodes(loaded), countNodes(g));
	BOOST_CHECK_EQUAL(countArcs(loaded), countArcs(g));

	MultiHypothesesGraph::Node n1;
	for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
			it != traxels.endValue(); ++it) {
		if (*it == t1[0]) {
			n1 = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
			break;
		}
	}
	std::vector<Traxel> loaded_traxels_n1 = loaded.get(node_regions_in_component())[n1];
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n1.begin(), loaded_traxels_n1.end(), t1.begin(), t1.end());
	string name = "conflicts";
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n1.front().features[name].begin(), loaded_traxels_n1.front().features[name].end(), t1.front().features[name].begin(), t1.front().features[name].end());
	name = "com";
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n1.front().features[name].begin(), loaded_traxels_n1.front().features[name].end(), t1.front().features[name].begin(), t1.front().features[name].end());
	name = "count_prediction";
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n1.front().features[name].begin(), loaded_traxels_n1.front().features[name].end(), t1.front().features[name].begin(), t1.front().features[name].end());
	BOOST_CHECK(loaded_traxels_n1.front().features.find("Variance") == loaded_traxels_n1.front().features.end());

	MultiHypothesesGraph::MoveFeatureMap& loaded_moves = loaded.get(node_move_features());
	std::map<Traxel, feature_array>& loaded_move_probs = loaded_moves.get_value(n1)[t1[0]];
	BOOST_CHECK_EQUAL_COLLECTIONS(move_probs[t2[0]].begin(), move_probs[t2[0]].end(), loaded_move_probs[t2[0]].begin(), loaded_move_probs[t2[0]].end());
	BOOST_CHECK_EQUAL_COLLECTIONS(move_probs[t2[1]].begin(), move_probs[t2[1]].end(), loaded_move_probs[t2[1]].begin(), loaded_move_probs[t2[1]].end());
	BOOST_CHECK_EQUAL_COLLECTIONS(move_probs[t2[2]].begin(), move_probs[t2[2]].end(), loaded_move_probs[t2[2]].begin(), loaded_move_probs[t2[2]].end());

	MultiHypothesesGraph::DivisionFeatureMap& loaded_divs = loaded.get(node_division_features());
	std::map<std::pair<Traxel, Traxel>, feature_array> loaded_div_probs = loaded_divs.get_value(n1)[t1[0]];
	BOOST_CHECK_EQUAL_COLLECTIONS(probs1.begin(), probs1.end(), loaded_div_probs[std::make_pair(t2[0], t2[1])].begin(), loaded_div_probs[std::make_pair(t2[0], t2[1])].end());
	BOOST_CHECK_EQUAL_COLLECTIONS(probs2.begin(), probs2.end(), loaded_div_probs[std::make_pair(t2[0], t2[2])].begin(), loaded_div_probs[std::make_pair(t2[0], t2[2])].end());



	MultiHypothesesGraph::Node n2;
	for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
			it != traxels.endValue(); ++it) {
		if (*it == t2[0]) {
			n2 = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
			break;
		}
	}
	std::vector<Traxel> loaded_traxels_n2 = loaded.get(node_regions_in_component())[n2];
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n2.begin(), loaded_traxels_n2.end(), t2.begin(), t2.end());
	name = "conflicts";
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n2.front().features[name].begin(), loaded_traxels_n2.front().features[name].end(), t2.front().features[name].begin(), t2.front().features[name].end());
	name = "com";
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n2.front().features[name].begin(), loaded_traxels_n2.front().features[name].end(), t2.front().features[name].begin(), t2.front().features[name].end());
	name = "count_prediction";
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n2.front().features[name].begin(), loaded_traxels_n2.front().features[name].end(), t2.front().features[name].begin(), t2.front().features[name].end());
}

BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_serialize_to_file ) {
	MultiHypothesesTraxelStore ts;

    size_t max_timesteps = 5;//1000000;
    size_t max_traxels_at = 5;//1000000;
    size_t shift = 1;
    for (size_t t = 0; t < max_timesteps; ++t) {
        for (size_t i = shift; i < max_traxels_at; ++i) {
        	std::vector<Traxel>& t1 = ts.map[t][shift + i];
        	std::vector<Traxel>& t2 = ts.map[t+1][i];

            std::vector<std::vector<unsigned> >& c1 = ts.conflicts_by_timestep[t][shift + i];
            std::vector<std::vector<unsigned> >& c2 = ts.conflicts_by_timestep[t+1][i];

            feature_array com(3, 1.);
            feature_array count(1, 0.1);

        	t1.push_back(Traxel(shift+i, t));
            float conflict_arr11[] = { 2., 3., 4., 5. };
            (t1.end() - 1)->features["conflicts"] = feature_array(conflict_arr11,
                    conflict_arr11 + 4);
            (t1.end() - 1)->features["level"].push_back(0.);
            (t1.end() - 1)->features["com"] = com;
            (t1.end() - 1)->features["count_prediction"] = count;

            c1.push_back(std::vector<unsigned>());
            c1.rbegin()->push_back(shift+i);

            t2.push_back(Traxel(i, t+1));
            float conflict_arr21[] = { 3., 4. };
            std::fill(com.begin(), com.end(), 0.);
            com[0] = 1;
            (t2.end() - 1)->features["conflicts"] = feature_array(conflict_arr21,
                    conflict_arr21 + 2);
            (t2.end() - 1)->features["level"].push_back(0.);
            (t2.end() - 1)->features["com"] = com;
            (t2.end() - 1)->features["count_prediction"] = count;

            t2.push_back(Traxel(i+1,t+1));
            float conflict_arr23[] = { 1. };
            com[0] = 0;
            (t2.end() - 1)->features["conflicts"] = feature_array(conflict_arr23,
                    conflict_arr23 + 1);
            (t2.end() - 1)->features["level"].push_back(1.);
            (t2.end() - 1)->features["com"] = com;

            t2.push_back(Traxel(i+2, t+1));
            com[0] = 2;
            float conflict_arr24[] = { 1. };
            (t2.end() - 1)->features["conflicts"] = feature_array(conflict_arr24,
                    conflict_arr24 + 1);
            (t2.end() - 1)->features["level"].push_back(1.);
            (t2.end() - 1)->features["com"] = com;

            c2.push_back(std::vector<unsigned>());
            c2.rbegin()->push_back(i);
            c2.rbegin()->push_back(i+1);
            c2.push_back(std::vector<unsigned>());
            c2.rbegin()->push_back(i);
            c2.rbegin()->push_back(i+2);

        }
        shift = max_traxels_at + 1;
    }
    std::cout << " -> workflow: building MultiHypothesesGraph" << std::endl;
    MultiHypothesesGraphBuilder::Options options(2, 1000000);
    MultiHypothesesGraphBuilder graph_builder(options);
    MultiHypothesesGraphPtr graph = graph_builder.build(ts);

    MultiHypothesesGraph& g = *graph;
    /*
    g.add(node_move_features());
    MultiHypothesesGraph::MoveFeatureMap& moves = g.get(node_move_features());
    const MultiHypothesesGraph::TraxelMap& traxels = g.get(node_traxel());
    MultiHypothesesGraph::Node n;
    for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
            it != traxels.endValue(); ++it) {
        if (*it == t1[0]) {
            n = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
            break;
        }
    }
    std::map<Traxel, feature_array>& move_probs = moves.get_value(n)[t1[0]];
    move_probs[t2[0]].push_back(0.3);
    move_probs[t2[0]].push_back(0.7);
    move_probs[t2[1]].push_back(0.99);
    move_probs[t2[1]].push_back(0.01);
    move_probs[t2[2]].push_back(0.05);
    move_probs[t2[2]].push_back(0.95);


    g.add(node_division_features());
    MultiHypothesesGraph::DivisionFeatureMap& divs = g.get(node_division_features());
    std::map<std::pair<Traxel, Traxel>, feature_array>& div_probs = divs.get_value(n)[t1[0]];
    std::vector<float> probs1, probs2;
    probs1.push_back(0.1);
    probs1.push_back(0.9);
    div_probs[std::make_pair(t2[0], t2[1])] = feature_array(probs1.size());
    std::copy(probs1.begin(),
      probs1.end(),
              div_probs[std::make_pair(t2[0], t2[1])].begin());
    probs2.push_back(0.1);
    probs2.push_back(0.9);
    div_probs[std::make_pair(t2[0], t2[2])] = feature_array(probs2.size());
    std::copy(probs2.begin(),
              probs2.end(),
              div_probs[std::make_pair(t2[0], t2[2])].begin());
              */

    std::string serialize_fn = "./serializationTest.dat";
    std::cout << "serialize to file" << std::endl;
    std::ofstream ofs(serialize_fn.c_str(), std::fstream::out | std::fstream::binary);
    // save data to archive
    {
          boost::archive::text_oarchive oa(ofs);
          // write class instance to archive
          oa << *graph;
        // archive and stream closed when destructors are called
    }


    std::cout << "deserialize from file" << std::endl;
    MultiHypothesesGraph loaded;
    {
        // create and open an archive for input
        std::ifstream ifs(serialize_fn.c_str(), std::fstream::binary | std::fstream::in);
        // read class state from archive
        boost::archive::text_iarchive ia(ifs);
        ia >> loaded;
        // archive and stream closed when destructors are called
    }


	loaded.get(node_traxel()); // shouldn't throw an exception
	loaded.get(node_conflict_sets()); // shouldn't throw an exception
	loaded.get(node_regions_in_component()); // shouldn't throw an exception

	BOOST_CHECK_EQUAL_COLLECTIONS(g.timesteps().begin(),
				g.timesteps().end(),
				loaded.timesteps().begin(),
				loaded.timesteps().end());

	BOOST_CHECK_EQUAL(countNodes(loaded), countNodes(g));
	BOOST_CHECK_EQUAL(countArcs(loaded), countArcs(g));

/*	MultiHypothesesGraph::Node n1;
	for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
			it != traxels.endValue(); ++it) {
		if (*it == t1[0]) {
			n1 = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
			break;
		}
	}
	std::vector<Traxel> loaded_traxels_n1 = loaded.get(node_regions_in_component())[n1];
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n1.begin(), loaded_traxels_n1.end(), t1.begin(), t1.end());
	string name = "conflicts";
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n1.front().features[name].begin(), loaded_traxels_n1.front().features[name].end(), t1.front().features[name].begin(), t1.front().features[name].end());
	name = "com";
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n1.front().features[name].begin(), loaded_traxels_n1.front().features[name].end(), t1.front().features[name].begin(), t1.front().features[name].end());
	name = "count_prediction";
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n1.front().features[name].begin(), loaded_traxels_n1.front().features[name].end(), t1.front().features[name].begin(), t1.front().features[name].end());

	MultiHypothesesGraph::MoveFeatureMap& loaded_moves = loaded.get(node_move_features());
	std::map<Traxel, feature_array>& loaded_move_probs = loaded_moves.get_value(n1)[t1[0]];
	BOOST_CHECK_EQUAL_COLLECTIONS(move_probs[t2[0]].begin(), move_probs[t2[0]].end(), loaded_move_probs[t2[0]].begin(), loaded_move_probs[t2[0]].end());
	BOOST_CHECK_EQUAL_COLLECTIONS(move_probs[t2[1]].begin(), move_probs[t2[1]].end(), loaded_move_probs[t2[1]].begin(), loaded_move_probs[t2[1]].end());
	BOOST_CHECK_EQUAL_COLLECTIONS(move_probs[t2[2]].begin(), move_probs[t2[2]].end(), loaded_move_probs[t2[2]].begin(), loaded_move_probs[t2[2]].end());

	MultiHypothesesGraph::DivisionFeatureMap& loaded_divs = loaded.get(node_division_features());
	std::map<std::pair<Traxel, Traxel>, feature_array> loaded_div_probs = loaded_divs.get_value(n1)[t1[0]];
	BOOST_CHECK_EQUAL_COLLECTIONS(probs1.begin(), probs1.end(), loaded_div_probs[std::make_pair(t2[0], t2[1])].begin(), loaded_div_probs[std::make_pair(t2[0], t2[1])].end());
	BOOST_CHECK_EQUAL_COLLECTIONS(probs2.begin(), probs2.end(), loaded_div_probs[std::make_pair(t2[0], t2[2])].begin(), loaded_div_probs[std::make_pair(t2[0], t2[2])].end());



	MultiHypothesesGraph::Node n2;
	for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
			it != traxels.endValue(); ++it) {
		if (*it == t2[0]) {
			n2 = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
			break;
		}
	}
	std::vector<Traxel> loaded_traxels_n2 = loaded.get(node_regions_in_component())[n2];
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n2.begin(), loaded_traxels_n2.end(), t2.begin(), t2.end());
	name = "conflicts";
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n2.front().features[name].begin(), loaded_traxels_n2.front().features[name].end(), t2.front().features[name].begin(), t2.front().features[name].end());
	name = "com";
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n2.front().features[name].begin(), loaded_traxels_n2.front().features[name].end(), t2.front().features[name].begin(), t2.front().features[name].end());
	name = "count_prediction";
	BOOST_CHECK_EQUAL_COLLECTIONS(loaded_traxels_n2.front().features[name].begin(), loaded_traxels_n2.front().features[name].end(), t2.front().features[name].begin(), t2.front().features[name].end());
*/
}

BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_pruning_nearest_neighbor ) {
  std::cout << "MultiHypothesesGraph_pruning_nearest_neighbor" << std::endl;
  MultiHypothesesGraph g;
  MultiHypothesesGraph::ContainedRegionsMap& regions = g.get(node_regions_in_component());

  const Node& n1 = g.add_node(1);
  const Node& n2 = g.add_node(1);
  const Node& n3 = g.add_node(2);
  const Node& n4 = g.add_node(2);

  g.addArc(n1, n3);
  g.addArc(n1, n4);
  g.addArc(n2, n3);
  g.addArc(n2, n4);

  feature_type com1[] = {2., 0., 0.};
  feature_type com2[] = {1., 0., 0.};
  feature_type com3[] = {3., 0., 0.};
  feature_type com4[] = {6., 0., 0.};
  feature_type com5[] = {5., 0., 0.};
  feature_type com6[] = {7., 0., 0.};

  feature_type count[] = {0.2, 0.7, 0.1};

  std::vector<Traxel>& t1 = regions.get_value(n1);
  std::vector<Traxel>& t2 = regions.get_value(n2);
  std::vector<Traxel>& t3 = regions.get_value(n3);
  std::vector<Traxel>& t4 = regions.get_value(n4);

  // generate traxels
  // t = 1
  t1.push_back(Traxel(1, 1));
  t1.rbegin()->features["com"] = feature_array(com1, com1 + 3);
  t1.rbegin()->features["count_prediction"] = feature_array(count, count + 3);
  t1.rbegin()->features["cardinality"] = feature_array(1, 2.);
  t1.push_back(Traxel(3, 1));
  t1.rbegin()->features["com"] = feature_array(com2, com2 + 3);
  t1.rbegin()->features["cardinality"] = feature_array(1, 1.);
  t1.push_back(Traxel(4, 1));
  t1.rbegin()->features["com"] = feature_array(com3, com3 + 3);
  t1.rbegin()->features["cardinality"] = feature_array(1, 1.);

  t2.push_back(Traxel(2, 1));
  t2.rbegin()->features["com"] = feature_array(com4, com4 + 3);
  t2.rbegin()->features["count_prediction"] = feature_array(count, count + 3);
  t2.rbegin()->features["cardinality"] = feature_array(1, 2.);
  t2.push_back(Traxel(5, 1));
  t2.rbegin()->features["com"] = feature_array(com5, com5 + 3);
  t2.rbegin()->features["cardinality"] = feature_array(1, 1.);
  t2.push_back(Traxel(6, 1));
  t2.rbegin()->features["com"] = feature_array(com6, com6 + 3);
  t2.rbegin()->features["cardinality"] = feature_array(1, 1.);

  // t = 2
  t3.push_back(Traxel(1, 2));
  t3.rbegin()->features["com"] = feature_array(com1, com1 + 3);
  t3.rbegin()->features["count_prediction"] = feature_array(count, count + 3);
  t3.rbegin()->features["cardinality"] = feature_array(1, 2.);
  t3.push_back(Traxel(3, 2));
  t3.rbegin()->features["com"] = feature_array(com2, com2 + 3);
  t3.rbegin()->features["cardinality"] = feature_array(1, 1.);
  t3.push_back(Traxel(4, 2));
  t3.rbegin()->features["com"] = feature_array(com3, com3 + 3);
  t3.rbegin()->features["cardinality"] = feature_array(1, 1.);

  t4.push_back(Traxel(2, 2));
  t4.rbegin()->features["com"] = feature_array(com4, com4 + 3);
  t4.rbegin()->features["count_prediction"] = feature_array(count, count + 3);
  t4.rbegin()->features["cardinality"] = feature_array(1, 2.);
  t4.push_back(Traxel(5, 2));
  t4.rbegin()->features["com"] = feature_array(com5, com5 + 3);
  t4.rbegin()->features["cardinality"] = feature_array(1, 1.);
  t4.push_back(Traxel(6, 2));
  t4.rbegin()->features["com"] = feature_array(com6, com6 + 3);
  t4.rbegin()->features["cardinality"] = feature_array(1, 1.);

  // arcs after pruning
  std::vector<TraxelArc> arcs;
  arcs.push_back(TraxelArc(Traxel(1,1), Traxel(1,2)));
  arcs.push_back(TraxelArc(Traxel(1,1), Traxel(5,2)));

  arcs.push_back(TraxelArc(Traxel(3,1), Traxel(3,2)));
  arcs.push_back(TraxelArc(Traxel(3,1), Traxel(5,2)));

  arcs.push_back(TraxelArc(Traxel(4,1), Traxel(4,2)));
  arcs.push_back(TraxelArc(Traxel(4,1), Traxel(5,2)));

  arcs.push_back(TraxelArc(Traxel(2,1), Traxel(2,2)));
  arcs.push_back(TraxelArc(Traxel(2,1), Traxel(4,2)));

  arcs.push_back(TraxelArc(Traxel(5,1), Traxel(5,2)));
  arcs.push_back(TraxelArc(Traxel(5,1), Traxel(4,2)));

  arcs.push_back(TraxelArc(Traxel(6,1), Traxel(6,2)));
  arcs.push_back(TraxelArc(Traxel(6,1), Traxel(4,2)));

  std::vector<bool> checks(arcs.size());

  ConstantFeature det(10);
  ConstantFeature mis(1000);
  ConstantFeature div(2000);
  pgm::multihypotheses::CVPR2014ModelBuilder builder( ConstantFeature(1000), // appearance
                                                      ConstantFeature(1000), // disappearance
                                                      SquaredDistance(), // move
                                                      ConstantFeature(0), // count
                                                      0, // forbidden_cost
                                                      0, // opportunity
                                                      50, // max_division_level
                                                      3 // max_count
                                                      );
  builder
      .with_detection_vars(det, mis)
      .with_divisions(div)
      .with_maximal_conflict_cliques(true)
      .with_maximum_arcs(1);
      ;

  MultiHypotheses reasoner(builder,
                           true, // with_constraints
                           0. // ep_gap
                           );

  std::cout << " -> workflow: formulating model" << std::endl;
  reasoner.formulate( g );
  const pgm::OpengmModel* model = reasoner.get_graphical_model();
  std::cout << "Checking the topology of the graphical model...\n";
  std::cout << "Number of variables: " << model->numberOfVariables() << '\n';
  BOOST_CHECK_EQUAL(model->numberOfVariables(), 24);

  const Model::arc_var_map& arc_var = reasoner.get_arc_map();
  for (Model::arc_var_map::const_iterator it = arc_var.begin(); it != arc_var.end(); ++it) {
    size_t idx = 0;
    for (std::vector<TraxelArc>::const_iterator a = arcs.begin(); a != arcs.end(); ++a, ++idx) {
      if (it->first.first == a->first && it->first.second == a->second) {
        break;
      }
    }
    // the required arc exists
    BOOST_REQUIRE(idx < checks.size());
    // the required arc has not been checked yet
    BOOST_REQUIRE(checks[idx] == false);
    checks[idx] = true;
  }
  // all arcs have been matched
  for (std::vector<bool>::const_iterator c = checks.begin(); c != checks.end(); ++c) {
    BOOST_CHECK(*c == true);
  }
}

#endif // ignore this test for now


// EOF
