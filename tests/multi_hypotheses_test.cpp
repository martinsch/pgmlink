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
#include "pgmlink/hypotheses.h"
#include "pgmlink/classifier_auxiliary.h"

using namespace pgmlink;
using namespace std;
using namespace boost;

typedef MultiHypothesesGraph::Node Node;
typedef MultiHypothesesGraph::Arc Arc;
typedef pgm::multihypotheses::Model Model;

BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_cardinalities ) {
  std::cout << "MultiHypothesesGraph_pruning_cardinalities" << std::endl;
  TraxelStore ts;
  feature_array com(3, 0.0);
  for (unsigned i = 1; i <= 5; ++i) {
    Traxel trax(i, 1);
    trax.features["com"] = com;
    add(ts, trax);
  }
  boost::shared_ptr<ConflictMap > conflicts(new ConflictMap);
  
  (*conflicts)[1].push_back(ConflictSet());
  (*conflicts)[1].rbegin()->push_back(1);
  (*conflicts)[1].rbegin()->push_back(5);
  
  (*conflicts)[1].push_back(ConflictSet());
  (*conflicts)[1].rbegin()->push_back(2);
  (*conflicts)[1].rbegin()->push_back(4);
  (*conflicts)[1].rbegin()->push_back(5);
  
  (*conflicts)[1].push_back(ConflictSet());
  (*conflicts)[1].rbegin()->push_back(3);
  (*conflicts)[1].rbegin()->push_back(4);
  (*conflicts)[1].rbegin()->push_back(5);

  std::map<Traxel, feature_type> cardinalities;
  cardinalities[Traxel(1, 1)] = 1.;
  cardinalities[Traxel(2, 1)] = 1.;
  cardinalities[Traxel(3, 1)] = 1.;
  cardinalities[Traxel(4, 1)] = 2.;
  cardinalities[Traxel(5, 1)] = 3.;

  SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(2, // max_nn
                                                               10000, // max_distance
                                                               true, // forward_backward
                                                               true, // consider_divisions
                                                               0.5 //division_threshold
                                                               );
  SingleTimestepTraxel_MultiHypothesesBuilder builder(&ts, builder_opts);
  MultiHypothesesGraphPtr g = builder.build_multi_hypotheses_graph();

  g->add_conflicts(conflicts);
  g->add_cardinalities();

  MultiHypothesesGraph::TraxelMap& traxel_map = g->get(node_traxel());
  for (MultiHypothesesGraph::NodeIt n(*g); n != lemon::INVALID; ++n) {
    const Traxel& trax = traxel_map[n]; 
    BOOST_CHECK_CLOSE(trax.features.find("cardinality")->second[0], cardinalities[trax], 0.001);
  }
}


BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_limit_arcs ) {
  std::cout << "MultiHypothesesGraph_pruning_limit_arcs" << std::endl;
  MultiHypothesesGraph g;
  MultiHypothesesGraph::TraxelMap& traxel_map = g.get(node_traxel());
  MultiHypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());
  ClassifierConstant mov(0.5, "move");
  ClassifierConstant div(0.5, "division");
  ClassifierConstant cnt(0.5, "count");
  ClassifierConstant det(0.5, "detProb");

  g.add_node(1);
  g.add_node(1);

  g.add_node(2);
  g.add_node(2);
  g.add_node(2);

  {
    unsigned id = 1;
    for (MultiHypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n, ++id) {
      Traxel trax = traxel_map[n];
      trax.Id = id;
      trax.Timestep = timestep_map[n];
      traxel_map.set(n, trax);
      if (trax.Timestep == 1) {
        for (MultiHypothesesGraph::node_timestep_map::ItemIt target(timestep_map, 2);
             target != lemon::INVALID;
             ++target) {
          g.addArc(n, target);
          std::cout << "ADDING ARC FROM " << traxel_map[n] << " TO " << traxel_map[target] << std::endl;
        }
      }
    }
  }

  g.add_classifier_features(&mov, &div, &cnt, &det);

  
  g.limit_arcs(2);

  {
    int count = 0;
    for (MultiHypothesesGraph::ArcIt arc(g); arc != lemon::INVALID; ++arc) ++count;
    BOOST_CHECK_EQUAL(count, 4);
    std::cout << "Total number of arcs in graph with outgoing arcs limited to 2: " << count << std::endl;
  }

  
  g.limit_arcs(1);

  {
    int count = 0;
    for (MultiHypothesesGraph::ArcIt arc(g); arc != lemon::INVALID; ++arc) ++count;
    BOOST_CHECK_EQUAL(count, 2);
    std::cout << "Total number of arcs in graph with outgoing arcs limited to 1: " << count << std::endl;
  }
}

// EOF
