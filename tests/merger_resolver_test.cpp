#define BOOST_TEST_MODULE merger_resolver_test

#include <stdexcept>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <set>
#include <vector>
#include <iterator>

#include <boost/test/unit_test.hpp>

#include <vigra/multi_array.hxx>
#include <vigra/tinyvector.hxx>

#include "pgmlink/hypotheses.h"
#include "pgmlink/traxels.h"
// enable tests of private class members
#define private public
#include "pgmlink/merger_resolving.h"
#undef private

#include <armadillo>
#include <mlpack/core.hpp>
#include <mlpack/methods/kmeans/kmeans.hpp>
#include <lemon/maps.h>
#include <lemon/concepts/digraph.h>

using namespace pgmlink;
using namespace std;
using namespace boost;


BOOST_AUTO_TEST_CASE( MergerResolver_no_mergers ) {
  HypothesesGraph src;
  HypothesesGraph dest;

  src.add(node_active())
    .add(arc_active())
    .add(node_active2())
    .add(node_traxel())
    .add(node_originated_from())
    .add(node_resolution_candidate())
    .add(arc_resolution_candidate())
    .add(arc_distance());

  MergerResolver m(&src);
  FeatureExtractorMCOMsFromGMM extractor(2);
  DistanceFromCOMs distance;
  FeatureHandlerFromTraxels handler(extractor, distance);
  m.resolve_mergers(handler);


  HypothesesGraph::Node n1 = src.add_node(1);
  HypothesesGraph::Node n2 = src.add_node(2);

  src.addArc(n1, n2);

  resolve_graph(src,
                dest,
                NegLnTransition(1), // weight 1
                0.05, // ep gap
                true, // with_tracklets
                100, // transition_parameter
                true // with_constraint
                );

  int n_arcs = lemon::countArcs(dest);
  int n_node = lemon::countNodes(dest);

  BOOST_CHECK_EQUAL(n_arcs, 0);
  BOOST_CHECK_EQUAL(n_node, 0);

}


BOOST_AUTO_TEST_CASE( MergerResolver_subgraph ) {
  HypothesesGraph g1, g2;
  g1.add(node_active()).add(arc_active()).add(node_active2()).add(node_traxel()).add(node_originated_from());
  HypothesesGraph::Node n1 = g1.add_node(1);
  HypothesesGraph::Node n2 = g1.add_node(2);
  HypothesesGraph::Node n3 = g1.add_node(2);
  property_map<node_active, HypothesesGraph::base_graph>::type& na_map = g1.get(node_active());
  property_map<node_active2, HypothesesGraph::base_graph>::type& na2_map = g1.get(node_active2());
  na_map.set(n1,true);
  na_map.set(n2,false);
  na_map.set(n3,true);
  na2_map.set(n1, 2);
  na2_map.set(n1, 3);
  na2_map.set(n1, 1);
  std::map<HypothesesGraph::Node, HypothesesGraph::Node> nr;
  std::map<HypothesesGraph::Arc, HypothesesGraph::Arc> ar;
  std::map<HypothesesGraph::Node, HypothesesGraph::Node> ncr;
  std::map<HypothesesGraph::Arc, HypothesesGraph::Arc> acr;
  copy_hypotheses_graph_subset<node_active, arc_active>(g1, g2, nr, ar, ncr, acr);
  g2.add(node_active()).add(node_active2());
  translate_property_bool_map<node_active, HypothesesGraph::Node>(g1, g2, nr);
  translate_property_value_map<node_active2, HypothesesGraph::Node>(g1, g2, nr);
  property_map<node_active, HypothesesGraph::base_graph>::type& n2a_map = g2.get(node_active());
  property_map<node_active2, HypothesesGraph::base_graph>::type& n2a2_map = g2.get(node_active2());

  std::map<HypothesesGraph::Node, HypothesesGraph::Node>::iterator it;
  for (it = nr.begin(); it != nr.end(); ++it) {
    std::cout << "Id (old): " << g1.id(it->first) <<  "(," <<  na_map[it->first] << "," << na2_map[it->first]
              << "), Id (new): " << g2.id(it->second) << "," << n2a_map[it->second] << "," << n2a2_map[it->second] << "\n";
  }

  for (it = ncr.begin(); it != ncr.end(); ++it) {
    std::cout << "Id (old): " << g1.id(it->second) << "(," << na_map[it->second] << "," << na2_map[it->second]
              << "), Id (new): " << g2.id(it->first) << ",(" << n2a_map[it->first] << "," << n2a2_map[it->first] << ")\n";
  }

  
  
}


BOOST_AUTO_TEST_CASE( MergerResolver_constructor ) {
  HypothesesGraph g;
  BOOST_CHECK_THROW(MergerResolver m(&g), std::runtime_error);
  g.add(node_active2());
  BOOST_CHECK_THROW(MergerResolver m(&g), std::runtime_error);
  g.add(arc_active());
  BOOST_CHECK_THROW(MergerResolver m(&g), std::runtime_error);
  g.add(arc_distance());
  MergerResolver m(&g);
  BOOST_CHECK_EQUAL(m.g_, &g);
  // check that merger_resolved_to property has been added
  BOOST_CHECK(m.g_->has_property(merger_resolved_to()));
  BOOST_CHECK(m.g_->has_property(node_originated_from()));

  // check exception on intialization with null pointer
  HypothesesGraph* G = 0; // = hyp_builder.build();
  BOOST_CHECK_THROW(MergerResolver M(G), std::runtime_error);
}


BOOST_AUTO_TEST_CASE( MergerResolver_resolve_mergers_3 ) {
  LOG(logINFO) << "Starting test MergerResolver_resolve_mergers_3";
  //  t=1       2
  //       --- (2)
  //     |
  //   (3) --- (1)

  // -> each of the nodes in timesteps t has a possible arc to all nodes in t+1
  //    o ----- o
  //     |     |
  //      -----
  //     | | | |
  //    o --x-- o
  //     | | | |
  //      -----
  //     |     |
  //    o ----- o

  HypothesesGraph g;
  g.add(node_traxel()).add(arc_distance()).add(arc_active()).add(node_active2());
  feature_array com(3,0);
  feature_array pCOM(6*3, 0);

  pCOM[0]  = 3;
  pCOM[3]  = 1;
  pCOM[6]  = 6;
  pCOM[9]  = 1;
  pCOM[12] = 3;
  pCOM[15] = 6;

  Traxel t11;
  t11.Timestep = 1;
  t11.Id = 11;
  com[0] = 3; t11.features["com"] = com;
  t11.features["possibleCOMs"] = pCOM;

  Traxel t21;
  t21.Timestep = 2;
  t21.Id = 21;
  com[0] = 1.5; t21.features["com"] = com;
  t21.features["possibleCOMs"] = pCOM;

  Traxel t22;
  t22.Timestep = 2;
  t22.Id = 22;
  com[0] = 3; t22.features["com"] = com;

  HypothesesGraph::Node n11 = g.add_node(1);
  HypothesesGraph::Node n21 = g.add_node(2);
  HypothesesGraph::Node n22 = g.add_node(2);

  HypothesesGraph::Arc a11_21 = g.addArc(n11, n21);
  HypothesesGraph::Arc a11_22 = g.addArc(n11, n22);



  property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
  traxel_map.set(n11, t11);
  traxel_map.set(n21, t21);
  traxel_map.set(n22, t22);

  property_map<arc_active, HypothesesGraph::base_graph>::type& arc_map = g.get(arc_active());
  arc_map.set(a11_21, true);
  arc_map.set(a11_22, true);

  property_map<node_active2, HypothesesGraph::base_graph>::type& active_map = g.get(node_active2());
  active_map.set(n11, 3);
  active_map.set(n21, 2);
  active_map.set(n22, 1);

  property_map<arc_distance, HypothesesGraph::base_graph>::type& dist_map = g.get(arc_distance());

  MergerResolver m(&g);
  FeatureExtractorMCOMsFromPCOMs extractor;
  DistanceFromCOMs distance;
  FeatureHandlerFromTraxels handler(extractor, distance);
  m.resolve_mergers(handler);
  prune_inactive(g);
  

  // check that arcs and nodes have been deactivated

  BOOST_CHECK(!g.valid(a11_22) || (g.source(a11_22) != n11 || g.target(a11_22) != n22));
  BOOST_CHECK(!g.valid(a11_21) || (g.source(a11_21) != n11 || g.target(a11_21) != n21));

  BOOST_CHECK(!g.valid(n11));
  BOOST_CHECK(!g.valid(n21));


  // check that traxel ids have been set correctly
  set<int> trx_ids_1;
  set<int> trx_ids_2;
  trx_ids_1.insert(12);
  trx_ids_1.insert(13);
  trx_ids_1.insert(14);
  trx_ids_1.insert(22);
  trx_ids_1.insert(23);
  trx_ids_1.insert(24);
  int trx_count = 0;
  property_map<node_active2, HypothesesGraph::base_graph>::type::ItemIt active_it(active_map, 1);
  for (; active_it != lemon::INVALID; ++active_it, ++trx_count) {
    trx_ids_2.insert(traxel_map[active_it].Id);

    HypothesesGraph::InArcIt IAIT(g, active_it);
    int count_n = 0;
    for(; IAIT != lemon::INVALID; ++IAIT, ++count_n) {
      BOOST_CHECK(arc_map[IAIT]);
    }
    HypothesesGraph::OutArcIt OAIT(g, active_it);
    for(; OAIT != lemon::INVALID; ++OAIT, ++count_n) {
      BOOST_CHECK(arc_map[OAIT]);
    }
    BOOST_CHECK_EQUAL(count_n, 3);
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(trx_ids_1.begin(), trx_ids_1.end(), trx_ids_2.begin(), trx_ids_2.end());
  BOOST_CHECK_EQUAL(trx_count, 6);

  
  // check that distances are calculated correctly and active arcs are valid
  set<double> arc_dist_1;
  set<double> arc_dist_2;
  arc_dist_1.insert(0);
  arc_dist_1.insert(2);
  arc_dist_1.insert(3);
  arc_dist_1.insert(5);
  int arc_count = 0;
  
  property_map<arc_active, HypothesesGraph::base_graph>::type::ItemIt arc_it(arc_map, true);
  for(; arc_it != lemon::INVALID; ++arc_it, ++arc_count) {
    arc_dist_2.insert(dist_map[arc_it]);
    BOOST_CHECK(g.valid(arc_it));
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(arc_dist_1.begin(), arc_dist_1.end(), arc_dist_2.begin(), arc_dist_2.end());
  BOOST_CHECK_EQUAL(arc_count, 9);


  

  // check that deactivated nodes are pruned, i.e. ItemIt(active_map, 0) should be equal  to lemon::INVALID
  property_map<node_active2, HypothesesGraph::base_graph>::type::ItemIt deactive_it(active_map, 0);
  BOOST_CHECK(!(deactive_it != lemon::INVALID));

  // check that deactivated arcs are pruned, i.e. FalseIt should be equal to lemon::INVALID
  property_map<arc_active, HypothesesGraph::base_graph>::type::FalseIt f_it(arc_map);
  BOOST_CHECK(!(f_it != lemon::INVALID));
      
  

  HypothesesGraph g_res;
  resolve_graph(g, g_res, NegLnTransition(1), 0.05, false);
  prune_inactive(g);

  vector<vector<Event> > ev = *(events(g));
  unsigned resolve_count = 0;
  for (vector<vector<Event> >::iterator t_it = ev.begin(); t_it != ev.end(); ++t_it) {
    for (vector<Event>::iterator e_it = t_it->begin(); e_it != t_it->end(); ++ e_it) {
      cout << *e_it << "\n";
      if (e_it->type == Event::ResolvedTo) ++resolve_count;
    }
  }
  BOOST_CHECK_EQUAL(resolve_count, 0);

}


BOOST_AUTO_TEST_CASE( MergerResolver_resolve_mergers_2 ) {
  LOG(logINFO) << "Starting test MergerResolver_resolve_mergers_2";
  
  //  t=1      2      3
  //    o ----    ----o
  //          | |
  //           O
  //          | |
  //    o ----    ----o

  // ->
  //    o ---- o ---- o
  //     |    | |    |
  //       --     --
  //     |    | |    |
  //    o ---- o ---- o

  HypothesesGraph g;
  g.add(node_traxel()).add(arc_distance()).add(arc_active()).add(node_active2()).add(merger_resolved_to());

  feature_array com(3,0);
  feature_array pCOM(6*3, 0);

  pCOM[0]  = 3;
  pCOM[3]  = 1;
  pCOM[6]  = 6;
  pCOM[9]  = 1;
  pCOM[12] = 3;
  pCOM[16] = 6;

  Traxel t11;
  t11.Timestep = 1;
  t11.Id = 11;
  com[0] = 1; t11.features["com"] = com;

  Traxel t12 = t11;
  t12.Id = 12;
  t12.features["com"][0] = 6;

  Traxel t21;
  t21.Timestep = 2;
  t21.Id = 21;
  com[0] = 3; t21.features["com"] = com;
  t21.features["possibleCOMs"] = pCOM;

  Traxel t31;
  t31.Timestep = 3;
  t31.Id = 31;
  com[0] = 1; t31.features["com"] = com;
  
  Traxel t32 = t31;
  t32.Id = 32;
  t32.features["com"][0] = 6;

  
  HypothesesGraph::Node n11 = g.add_node(1);
  HypothesesGraph::Node n12 = g.add_node(1);
  HypothesesGraph::Node n21 = g.add_node(2);
  HypothesesGraph::Node n31 = g.add_node(3);
  HypothesesGraph::Node n32 = g.add_node(3);

  HypothesesGraph::Arc a11_21 = g.addArc(n11, n21);
  HypothesesGraph::Arc a12_21 = g.addArc(n12, n21);
  HypothesesGraph::Arc a21_31 = g.addArc(n21, n31);
  HypothesesGraph::Arc a21_32 = g.addArc(n21, n32);

  property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
  traxel_map.set(n11, t11);
  traxel_map.set(n12, t12);
  traxel_map.set(n21, t21);
  traxel_map.set(n31, t31);
  traxel_map.set(n32, t32);

  property_map<arc_active, HypothesesGraph::base_graph>::type& arc_map = g.get(arc_active());
  arc_map.set(a11_21, true);
  arc_map.set(a12_21, true);
  arc_map.set(a21_31, true);
  arc_map.set(a21_32, true);

  

  property_map<node_active2, HypothesesGraph::base_graph>::type& active_map = g.get(node_active2());
  active_map.set(n11, 1);
  active_map.set(n12, 1);
  active_map.set(n21, 2);
  active_map.set(n31, 1);
  active_map.set(n32, 1);

  MergerResolver m(&g);
  FeatureExtractorMCOMsFromPCOMs extractor;
  DistanceFromCOMs distance;
  FeatureHandlerFromTraxels handler(extractor, distance);
  m.resolve_mergers(handler);
  prune_inactive(g);
  // property_map<merger_resolved_to, HypothesesGraph::base_graph>::type& resolved_map = g.get(merger_resolved_to());
  property_map<node_timestep, HypothesesGraph::base_graph>::type& time_map = g.get(node_timestep());
  property_map<node_timestep, HypothesesGraph::base_graph>::type::ItemIt IT(time_map, 2);

  // deactivated arcs from and to merger node
  BOOST_CHECK(!g.valid(a11_21) || g.source(a11_21) != n11 || g.target(a11_21) != n21);
  BOOST_CHECK(!g.valid(a12_21) || g.source(a12_21) != n12 || g.target(a12_21) != n21);
  BOOST_CHECK(!g.valid(a21_31) || g.source(a21_31) != n21 || g.target(a21_31) != n31);
  BOOST_CHECK(!g.valid(a21_32) || g.source(a21_32) != n21 || g.target(a21_32) != n32);
  
  vector<vector<Event> > ev = *(events(g));
  unsigned time = 0;
  unsigned resolve_count = 0;
  // LOG(logINFO) << "Detected the following events:";
  for (vector<vector<Event> >::iterator it = ev.begin(); it != ev.end(); ++it, ++time) {
    for (vector<Event>::iterator It = it->begin(); It != it->end(); ++It) {
      // LOG(logINFO) << " " << time << ": " << *It;
      if (It->type == Event::ResolvedTo) ++resolve_count;
    }
  } 
  BOOST_CHECK_EQUAL(resolve_count, 1);

  HypothesesGraph g_res;
  resolve_graph(g, g_res, NegLnTransition(1), 0.05, false);
  prune_inactive(g);
  vector<vector<Event> > evt = *(events(g));
  for (vector<vector<Event> >::iterator t_it = evt.begin(); t_it != evt.end(); ++t_it) {
    for (vector<Event>::iterator e_it = t_it->begin(); e_it != t_it->end(); ++ e_it) {
      cout << "t=" << t_it - evt.begin() + 1 << ": " << *e_it << "\n";
    }
  }

  vector<vector<Event> > evts = *(multi_frame_move_events(g));
  for (vector<vector<Event> >::iterator t_it = evts.begin(); t_it != evts.end(); ++t_it) {
    for (vector<Event>::iterator e_it = t_it->begin(); e_it != t_it->end(); ++ e_it) {
      cout << "t=" << t_it - evts.begin() + 1 << ": " << *e_it << "\n";
    }
  }

  vector<vector<Event> > merged_ev = *merge_event_vectors(evt, evts);
  for (vector<vector<Event> >::iterator t_it = merged_ev.begin(); t_it != merged_ev.end(); ++t_it) {
    size_t pos = t_it - merged_ev.begin();
    BOOST_CHECK_EQUAL(t_it->size(), (evt.begin()+pos)->size() + (evts.begin()+pos)->size());
  }
}


BOOST_AUTO_TEST_CASE( MergerResolver_resolve_mergers ) {
  LOG(logINFO) << "Starting test MergerResolver_resolve_mergers";
  HypothesesGraph g;
  g.add(node_traxel()).add(arc_distance()).add(arc_active()).add(node_active2()).add(merger_resolved_to());
  //  t=1      2
  //    o ----  
  //          | 
  //           O
  //          | 
  //    o ----   

  // ->
  //    o ---- o
  //     |    | 
  //       --   
  //     |    | 
  //    o ---- o

  feature_array com(3,0);
  feature_array pCOM(6*3, 0);

  pCOM[0]  = 3;
  pCOM[3]  = 1;
  pCOM[6]  = 6;
  pCOM[9]  = 1;
  pCOM[12] = 3;
  pCOM[16] = 6;

  Traxel t11;
  t11.Timestep = 1;
  t11.Id = 11;
  com[0] = 1; t11.features["com"] = com;

  Traxel t12 = t11;
  t12.Id = 12;
  t12.features["com"][0] = 6;

  Traxel t21;
  t21.Timestep = 2;
  t21.Id = 21;
  com[0] = 3; t21.features["com"] = com;
  t21.features["possibleCOMs"] = pCOM;

  
  HypothesesGraph::Node n11 = g.add_node(1);
  HypothesesGraph::Node n12 = g.add_node(1);
  HypothesesGraph::Node n21 = g.add_node(2);

  HypothesesGraph::Arc a11_21 = g.addArc(n11, n21);
  HypothesesGraph::Arc a12_21 = g.addArc(n12, n21);

  property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
  traxel_map.set(n11, t11);
  traxel_map.set(n12, t12);
  traxel_map.set(n21, t21);

  property_map<arc_active, HypothesesGraph::base_graph>::type& arc_map = g.get(arc_active());
  arc_map.set(a11_21, true);
  arc_map.set(a12_21, true);

  

  property_map<node_active2, HypothesesGraph::base_graph>::type& active_map = g.get(node_active2());
  active_map.set(n11, 1);
  active_map.set(n12, 1);
  active_map.set(n21, 2);

  MergerResolver m(&g);
  FeatureExtractorMCOMsFromPCOMs extractor;
  DistanceFromCOMs distance;
  FeatureHandlerFromTraxels handler(extractor, distance);
  m.resolve_mergers(handler);
  prune_inactive(g);

  // setup tests
  property_map<node_active2, HypothesesGraph::base_graph>::type::ValueIt active_valueIt = active_map.beginValue();
  property_map<arc_distance, HypothesesGraph::base_graph>::type& distance_map = g.get(arc_distance());
  std::set<int> active_values;
  std::set<int> values_found;
  active_values.insert(1);
  // active_values.insert(2);
  std::set<double> distances;
  distances.insert(0);
  distances.insert(5);
  std::set<int> new_ids;
  new_ids.insert(22);
  new_ids.insert(23);
  std::set<int> set_ids;
  
  for (; active_valueIt != active_map.endValue(); ++active_valueIt) {
    property_map<node_active2, HypothesesGraph::base_graph>::type::ItemIt active_itemIt(active_map, *active_valueIt);
    property_map<merger_resolved_to, HypothesesGraph::base_graph>::type& to_map = g.get(merger_resolved_to());
    values_found.insert(*active_valueIt);
    for (; active_itemIt != lemon::INVALID; ++active_itemIt) {
      int count = 0;
      std::set<double> measured_distances;
      for (HypothesesGraph::InArcIt it(g, active_itemIt); it != lemon::INVALID; ++it) {
	++count;
	measured_distances.insert(distance_map[it]);
      }
      for (HypothesesGraph::OutArcIt it(g, active_itemIt); it != lemon::INVALID; ++it) {
	++count;
	measured_distances.insert(distance_map[it]);
      }
      if (*active_valueIt == 1) {
	BOOST_CHECK_EQUAL(count, 2);
	BOOST_CHECK_EQUAL_COLLECTIONS(distances.begin(), distances.end(), measured_distances.begin(), measured_distances.end());
      } else if (*active_valueIt == 0) {
	BOOST_CHECK_EQUAL(count, 0);
	BOOST_CHECK_EQUAL(traxel_map[active_itemIt].Id, 21);
	set_ids = std::set<int>(to_map[active_itemIt].begin(), to_map[active_itemIt].end());
      }
    }
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(values_found.begin(), values_found.end(), active_values.begin(), active_values.end());
  // BOOST_CHECK_EQUAL_COLLECTIONS(set_ids.begin(), set_ids.end(), new_ids.begin(), new_ids.end());

  vector<vector<Event> > ev = *(events(g));
  unsigned time = 0;
  LOG(logINFO) << "Detected the following events:";
  for (vector<vector<Event> >::iterator it = ev.begin(); it != ev.end(); ++it, ++time) {
    for (vector<Event>::iterator It = it->begin(); It != it->end(); ++It) {
      LOG(logINFO) << " " << time << ": " << *It;
    }
  }
  
}


BOOST_AUTO_TEST_CASE( MergerResolver_refine_node ) {
  LOG(logINFO) << "Starting test MergerResolver_refine_node";
  // MergerResolver::refine_node(HypothesesGraph::Node node, std::size_t nMerger)

  //  t=1      2      3 
  //    o ----   --- o
  //          | |
  //           O
  //          | |
  //    o ----   --- o

  // ->
  //    o ---- o ---- o
  //     |    | |    |
  //       --     --
  //     |    | |    |
  //    o ---- o ---- o

  
  HypothesesGraph g;
  // g.add(arc_distance()).add(tracklet_intern_dist()).add(node_tracklet()).add(tracklet_intern_arc_ids()).add(traxel_arc_id());
  g.add(node_traxel()).add(arc_distance()).add(arc_active()).add(node_active2()).add(merger_resolved_to());
  
  feature_array com(3,0);
  feature_array pCOM(6*3, 0);

  pCOM[0]  = 3;
  pCOM[3]  = 1;
  pCOM[6]  = 6;
  pCOM[9]  = 1;
  pCOM[12] = 3;
  pCOM[16] = 6;

  feature_array mCOM(pCOM.begin()+3, pCOM.begin()+9);
  feature_array com1(mCOM.begin(), mCOM.begin()+3);
  feature_array com2(mCOM.begin()+3, mCOM.end());
  
  std::set<double> dist;
  dist.insert(0);
  dist.insert(5);

  std::vector<unsigned int> new_ids;
  new_ids.push_back(22);
  new_ids.push_back(23);
    
  Traxel t11;
  t11.Timestep = 1;
  t11.Id = 11;
  com[0] = 1; t11.features["com"] = com;

  Traxel t12 = t11;
  t12.Id = 12;
  t12.features["com"][0] = 6;

  Traxel t21;
  t21.Timestep = 2;
  t21.Id = 21;
  com[0] = 3; t21.features["com"] = com;
  t21.features["mergerCOMs"] = mCOM;

  
  Traxel t31;
  t31.Timestep = 3;
  t31.Id = 31;
  com[0] = 1; t31.features["com"] = com;
  
  Traxel t32 = t31;
  t32.Id = 32;
  t32.features["com"][0] = 6;
  
  HypothesesGraph::Node n11 = g.add_node(1);
  HypothesesGraph::Node n12 = g.add_node(1);
  HypothesesGraph::Node n21 = g.add_node(2);
  HypothesesGraph::Node n31 = g.add_node(3);
  HypothesesGraph::Node n32 = g.add_node(3);

  HypothesesGraph::Arc a11_21 = g.addArc(n11, n21);
  HypothesesGraph::Arc a12_21 = g.addArc(n12, n21);
  HypothesesGraph::Arc a21_31 = g.addArc(n21, n31);
  HypothesesGraph::Arc a21_32 = g.addArc(n21, n32);

  property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
  traxel_map.set(n11, t11);
  traxel_map.set(n12, t12);
  traxel_map.set(n21, t21);
  traxel_map.set(n31, t31);
  traxel_map.set(n32, t32);

  property_map<arc_active, HypothesesGraph::base_graph>::type& arc_map = g.get(arc_active());
  arc_map.set(a11_21, true);
  arc_map.set(a12_21, true);
  arc_map.set(a21_31, true);
  arc_map.set(a21_32, true);

  property_map<node_active2, HypothesesGraph::base_graph>::type& active_map = g.get(node_active2());
  active_map.set(n11, 1);
  active_map.set(n12, 1);
  active_map.set(n21, 2);
  active_map.set(n31, 1);
  active_map.set(n32, 1);

  property_map<node_timestep, HypothesesGraph::base_graph>::type& timestep_map = g.get(node_timestep());
  property_map<arc_distance, HypothesesGraph::base_graph>::type& distance_map = g.get(arc_distance());
  

  MergerResolver m(&g);
  FeatureExtractorMCOMsFromMCOMs extractor;
  DistanceFromCOMs distance;
  FeatureHandlerFromTraxels handler(extractor, distance);
  m.refine_node(n21, 2, handler);
  prune_inactive(g);
  
  // deactivated arcs from and to merger node
  BOOST_CHECK(!g.valid(a11_21) || g.source(a11_21) != n11 || g.target(a11_21) != n21);
  BOOST_CHECK(!g.valid(a12_21) || g.source(a12_21) != n12 || g.target(a12_21) != n21);
  BOOST_CHECK(!g.valid(a21_31) || g.source(a21_31) != n21 || g.target(a21_31) != n31);
  BOOST_CHECK(!g.valid(a21_32) || g.source(a21_32) != n21 || g.target(a21_32) != n32);

  int count = 0;
  property_map<node_timestep, HypothesesGraph::base_graph>::type::ItemIt timeIt(timestep_map, 2);
  for (; timeIt != lemon::INVALID; ++timeIt) {
    HypothesesGraph::base_graph::InArcIt inIt(g, timeIt);
    HypothesesGraph::base_graph::OutArcIt outIt(g, timeIt);
    if (timeIt == n21) {
      
      // merger node still active, will be deactivated in MergerResolver::resolve_mergers
      BOOST_CHECK_EQUAL(active_map[timeIt], 2);
    } else {

      // if correct merger COM exists in traxel, count++
      Traxel trax = traxel_map[timeIt];
      com = trax.features["com"];
      if ( (std::equal(com1.begin(), com1.end(), com.begin()) ||
	    std::equal(com2.begin(), com2.end(), com.begin()) ) &&
	   (trax.Id == 22 || trax.Id == 23) ) {
	++count;
      }
      
      
      std::set<double> inDist;
      for (; inIt != lemon::INVALID; ++inIt) {
	// save distances of arcs
	inDist.insert(distance_map[inIt]);
	// check if activated
        BOOST_CHECK(arc_map[inIt]);
      }
      // check distances of new arcs
      BOOST_CHECK_EQUAL_COLLECTIONS(inDist.begin(), inDist.end(), dist.begin(), dist.end());
      
      std::set<double> outDist;
      for (; outIt != lemon::INVALID; ++outIt) {
	// save distances of arcs
	outDist.insert(distance_map[outIt]);
	// check if activated
        BOOST_CHECK(arc_map[outIt]);	
      }
      // check distances of new arcs
      BOOST_CHECK_EQUAL_COLLECTIONS(outDist.begin(), outDist.end(), dist.begin(), dist.end());
      
      // activated merger replacement nodes
      BOOST_CHECK_EQUAL(active_map[timeIt], 1);
    }
  }
  // check number of newly created nodes
  BOOST_CHECK_EQUAL(count, 2);
  // check that merger_resolved_to property is set properly
  std::vector<unsigned int> property = g.get(merger_resolved_to())[n21];
  BOOST_CHECK_EQUAL_COLLECTIONS(property.begin(), property.end(), new_ids.begin(), new_ids.end());
  
}


BOOST_AUTO_TEST_CASE( MergerResolver_deactivate_arcs ) {
  // deactivate_arcs(std::vector<HypothesesGraph::base_graph::Arc>)
  
  //  t=1      2      3 
  //    o ---- o ---- o
             
  // ->
  //    o      o      o

  
  HypothesesGraph g;
  g.add(arc_active()).add(node_active2()).add(arc_distance()).add(node_traxel());
  HypothesesGraph::Node n1 = g.add_node(1);
  HypothesesGraph::Node n2 = g.add_node(2);
  HypothesesGraph::Node n3 = g.add_node(3);
  HypothesesGraph::Arc a12 = g.addArc(n1, n2);
  HypothesesGraph::Arc a23 = g.addArc(n2, n3);

  g.get(arc_active()).set(a12, true);
  g.get(arc_active()).set(a23, true);
  
  std::vector<HypothesesGraph::base_graph::Arc> arcs;
  arcs.push_back(a12);
  arcs.push_back(a23);
  MergerResolver m(&g);
  m.deactivate_arcs(arcs);
  for (std::vector<HypothesesGraph::base_graph::Arc>::iterator it = arcs.begin(); it != arcs.end(); ++it) {
    // arc is deactivated
    BOOST_CHECK_EQUAL(g.get(arc_active())[*it], false);
  }
}


BOOST_AUTO_TEST_CASE( MergerResolver_deactivate_nodes ) {
  LOG(logINFO) << "Starting test MergerResolver_deactivate_nodes";
  // deactivate_nodes(std::vector<HypothesesGraph::Node> nodes)

  // (1) -> (0)

  
  HypothesesGraph g;
  g.add(node_active2()).add(arc_active()).add(arc_distance());
  HypothesesGraph::Node n = g.add_node(1);
  g.get(node_active2()).set(n, 1);
  set<int> active_values;
  set<int> active_values_from_properties;
  active_values.insert(0);
  
  // node is active with 1 object
  BOOST_CHECK_EQUAL(g.get(node_active2())[n], 1);

  std::vector<HypothesesGraph::Node> nodes;
  nodes.push_back(n);
  MergerResolver m(&g);
  m.deactivate_nodes(nodes);
  property_map<node_active2, HypothesesGraph::base_graph>::type& active_map = g.get(node_active2());
  property_map<node_active2, HypothesesGraph::base_graph>::type::ValueIt active_It=active_map.beginValue();
  for (; active_It != active_map.endValue(); ++active_It) active_values_from_properties.insert(*active_It);
  BOOST_CHECK_EQUAL_COLLECTIONS(active_values.begin(), active_values.end(), active_values_from_properties.begin(), active_values_from_properties.end());
  
  // node is not active (0) objects
  BOOST_CHECK_EQUAL(g.get(node_active2())[n], 0);
}


BOOST_AUTO_TEST_CASE( MergerResolver_get_max_id ) {
  // get_max_id(int ts)

  // ids in timesteps:
  //  t=1      2
  //    1      2
  //    2
  
  HypothesesGraph g;
  g.add(node_traxel()).add(node_active2()).add(arc_active()).add(arc_distance());
  Traxel t11, t12, t21;
  
  t11.Timestep = 1;
  t11.Id = 1;

  t12.Timestep = 1;
  t12.Id = 2;

  t21.Timestep = 2;
  t21.Id = 2;

  HypothesesGraph::Node n11 = g.add_node(1);
  HypothesesGraph::Node n12 = g.add_node(1);
  HypothesesGraph::Node n21 = g.add_node(2);

  property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());

  traxel_map.set(n11, t11);
  traxel_map.set(n12, t12);
  traxel_map.set(n21, t21);

  MergerResolver m(&g);
  BOOST_CHECK_EQUAL(m.get_max_id(1), 2);
  BOOST_CHECK_EQUAL(m.get_max_id(2), 2);  
}


BOOST_AUTO_TEST_CASE( MergerResolver_add_arcs_for_replacement_node ) {
  // BOOST_CHECK(false); // FIXME: fix this test!
//  // MergerResolver::add_arcs_for_replacement_node(HypothesesGraph::Node node,
//  //						   Traxel trax,
//  // 						   std::vector<HypothesesGraph::base_graph::Arc> src,
//  //						   std::vector<HypothesesGraph::base_graph::Arc> dest);
//
//  //  t=1      2      3
//  //    o ---- o ---- o
//
//  // ->
//  //    o      o      o
//  //     \           /
//  //      ---- o ----
//
//
//  feature_array COM(3, 0.0);
//  TraxelStore ts;
//  Traxel t1, t2, t3, t22;
//
//  t1.Id = 11;
//  t1.Timestep = 1;
//  t1.features["com"] = COM;
//
//  t2.Id = 21;
//  t2.Timestep = 2;
//  t2.features["com"] = COM;
//
//  t3.Id = 31;
//  t3.Timestep = 3;
//  t3.features["com"] = COM;
//
//  COM[0] = 100;
//  t22.Id = 22;
//  t22.Timestep = 2;
//  t22.features["com"] = COM;
//
//  add(ts, t1);
//  add(ts, t2);
//  add(ts, t3);
//  SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(1, // max NN
//							       1, // max Dist
//							       false, // forward_backward
//							       false, //  consider_divisions
//							       100 // division_threshold
//							       );
//  SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
//  HypothesesGraph* g = hyp_builder.build();
//  HypothesesGraph& G = *g;
//  G.add(arc_distance()).add(arc_active()).add(node_active2());
//  MergerResolver m(g);
//  DistanceFromCOMs distance;
//
//  property_map<node_timestep, HypothesesGraph::base_graph>::type& time_map = g->get(node_timestep());
//  property_map<node_timestep, HypothesesGraph::base_graph>::type::ItemIt time_it(time_map, 2);
//  property_map<node_traxel, HypothesesGraph::base_graph>::type& trax_map = g->get(node_traxel());
//  HypothesesGraph::Node newNode = g->add_node(t22.Timestep);
//  trax_map.set(newNode, t22);
//  for (; time_it != lemon::INVALID; ++time_it) {
//    HypothesesGraph::Node node = time_it;
//    std::vector<HypothesesGraph::base_graph::Arc> sources;
//    std::vector<HypothesesGraph::base_graph::Arc> targets;
//    m.collect_arcs(HypothesesGraph::base_graph::InArcIt(*g, node), sources);
//    m.collect_arcs(HypothesesGraph::base_graph::OutArcIt(*g, node), targets);
//    if (sources.size() > 0 && targets.size() > 0) {
//      m.add_arcs_for_replacement_node(newNode, t22, sources, targets, distance);
//    }
//  }
//  property_map<arc_distance, HypothesesGraph::base_graph>::type& distance_map = g->get(arc_distance());
//  property_map<arc_active, HypothesesGraph::base_graph>::type& active_map = g->get(arc_active());
//  HypothesesGraph::base_graph::InArcIt inIt(G, newNode);
//  int k = 0;
//  for (; inIt != lemon::INVALID; ++k, ++inIt) {
//    BOOST_CHECK_EQUAL(distance_map[inIt], 100);
//    BOOST_CHECK(active_map[inIt]);
//  }
//  BOOST_CHECK_EQUAL(k, 1);
//
//  HypothesesGraph::base_graph::OutArcIt outIt(G, newNode);
//  k = 0;
//  for (; outIt != lemon::INVALID; ++k, ++outIt) {
//    BOOST_CHECK_EQUAL(distance_map[outIt], 100);
//    BOOST_CHECK(active_map[outIt]);
//  }
//  BOOST_CHECK_EQUAL(k, 1);
}


BOOST_AUTO_TEST_CASE( MergerResolver_collect_arcs ) {
  // MergerResolver::collect_arcs(ArcIterator arcIt, std::vector<HypothesesGraph::base_graph::Arc>& res)

  //  t=1      2      3 
  //    o ---- o ---- o 

  
  feature_array COM(3, 0.0);
  TraxelStore ts;
  Traxel t1, t2, t3;
  
  t1.Id = 11;
  t1.Timestep = 1;
  t1.features["com"] = COM;
  
  t2.Id = 21;
  t2.Timestep = 2;
  t2.features["com"] = COM;
  
  t3.Id = 31;
  t3.Timestep = 3;
  t3.features["com"] = COM;
  
  add(ts, t1);
  add(ts, t2);
  add(ts, t3);
  SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(1, // max NN
							       100, // max Dist
							       false, // forward_backward
							       false, // consider_divisions
							       100 // divison_threshold
							       );
  SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
  HypothesesGraph* g = hyp_builder.build();
  g->add(node_active2()).add(arc_active()).add(arc_distance());
  MergerResolver m(g);
  
  std::vector<HypothesesGraph::base_graph::Arc> sources;
  std::vector<HypothesesGraph::base_graph::Arc> targets;
  
  property_map<node_timestep, HypothesesGraph::base_graph>::type& time_map = g->get(node_timestep());
  property_map<node_timestep, HypothesesGraph::base_graph>::type::ItemIt time_it(time_map, 2);
  property_map<arc_active, HypothesesGraph::base_graph>::type& arc_active_map = g->get(arc_active());
  for (HypothesesGraph::base_graph::ArcIt arcIt(*g); arcIt != lemon::INVALID; ++arcIt) {
    arc_active_map.set(arcIt, true);
  }
  HypothesesGraph::Node node = time_it;
  m.collect_arcs(HypothesesGraph::base_graph::InArcIt(*g, node), sources);
  m.collect_arcs(HypothesesGraph::base_graph::OutArcIt(*g, node), targets);

  // There should be one arc from t1 to t2 and one from t2 to t3
  BOOST_CHECK_EQUAL(sources.size(), 1);
  BOOST_CHECK_EQUAL(targets.size(), 1);

  property_map<node_traxel, HypothesesGraph::base_graph>::type& trax_map = g->get(node_traxel());
  Traxel src_from = trax_map[g->source(sources[0])];
  Traxel src_to = trax_map[g->target(sources[0])];
  Traxel tar_from = trax_map[g->source(targets[0])];
  Traxel tar_to = trax_map[g->target(targets[0])];

  BOOST_CHECK_EQUAL(src_from.Id, t1.Id);
  BOOST_CHECK_EQUAL(src_to.Id, t2.Id);
  BOOST_CHECK_EQUAL(tar_from.Id, t2.Id);
  BOOST_CHECK_EQUAL(tar_to.Id, t3.Id);
}


BOOST_AUTO_TEST_CASE( MergerResolver_kmeans) {
  // KMeans(int k, const feature_array& data, feature_array& centers)
  // KMeans::operator()
  float arr[] = {-6, 0, 0, -5, 0, 0, -4, 0, 0, 6, 0, 0, 5, 0, 0, 4, 0, 0};
  float arr_res[] = {5, 0, 0, -5, 0, 0};
  feature_array data(arr, arr + sizeof(arr)/sizeof(arr[0]));
  KMeans kMeans(2, data);
  feature_array centers = kMeans();
  BOOST_CHECK_EQUAL_COLLECTIONS(centers.begin(), centers.end(), arr_res, arr_res+sizeof(arr_res)/sizeof(arr_res[0]));
}


BOOST_AUTO_TEST_CASE( MergerResolver_extract_coordinates ) {
  std::cout << "MergerResolver_extract_coordinates" << std::endl;
  vigra::MultiArray<2, unsigned> label_image2D(vigra::Shape2(10,10), 0u);
  vigra::MultiArray<3, unsigned> label_image3D(vigra::Shape3(10,10,3), 0u);

  int timestep = 5;

  TimestepIdCoordinateMap coordinate_base2D;
  TimestepIdCoordinateMap coordinate_base3D;

  TimestepIdCoordinateMapPtr coordinate_result2D(new TimestepIdCoordinateMap);
  TimestepIdCoordinateMapPtr coordinate_result3D(new TimestepIdCoordinateMap);

  long int lower1, lower2, upper1, upper2, offset1 = 2, offset2 = 0;
  unsigned label;
  vigra::TinyVector<long int, 2> offset2D(offset1, offset2);
  vigra::TinyVector<long int, 3> offset3D(offset1, offset2, 0);


  std::cout << "MergerResolver_extract_coordinates -- create label images and coordinate lists" << std::endl;

  lower1 = 0;
  upper1 = 6;
  lower2 = 0;
  upper2 = 4;
  label  = 3;
  Traxel trax1(label, timestep);
  trax1.features["count"] = feature_array(1, (upper1-lower1)*(upper2-lower2));

  coordinate_base2D[std::make_pair(timestep, label)] = arma::mat(2, (upper1-lower1)*(upper2-lower2));
  coordinate_base3D[std::make_pair(timestep, label)] = arma::mat(3, (upper1-lower1)*(upper2-lower2));
  {
    int count = 0;
    for (long idx2 = lower2; idx2 < upper2; ++idx2) {
      for (long idx1 = lower1; idx1 < upper1; ++idx1, ++count) {
        label_image2D(idx1, idx2) = label;
        label_image3D(idx1, idx2, 1) = label;
        coordinate_base2D[std::make_pair(timestep, label)](0, count) = idx1 + offset1;
        coordinate_base2D[std::make_pair(timestep, label)](1, count) = idx2 + offset2;

        coordinate_base3D[std::make_pair(timestep, label)](0, count) = idx1 + offset1;
        coordinate_base3D[std::make_pair(timestep, label)](1, count) = idx2 + offset2;
        coordinate_base3D[std::make_pair(timestep, label)](2, count) = 1;
      }
    }
  }

  lower1 = 1;
  upper1 = 4;
  lower2 = 7;
  upper2 = 10;
  label  = 15;
  Traxel trax2(label, timestep);
  trax2.features["count"] = feature_array(1, (upper1-lower1)*(upper2-lower2));
  
  coordinate_base2D[std::make_pair(timestep, label)] = arma::mat(2, (upper1-lower1)*(upper2-lower2));
  coordinate_base3D[std::make_pair(timestep, label)] = arma::mat(3, (upper1-lower1)*(upper2-lower2));
  {
    int count = 0;
    for (long idx2 = lower2; idx2 < upper2; ++idx2) {
      for (long idx1 = lower1; idx1 < upper1; ++idx1, ++count) {
        label_image2D(idx1, idx2) = label;
        label_image3D(idx1, idx2, 1) = label;
        coordinate_base2D[std::make_pair(timestep, label)](0, count) = idx1 + offset1;
        coordinate_base2D[std::make_pair(timestep, label)](1, count) = idx2 + offset2;
      
        coordinate_base3D[std::make_pair(timestep, label)](0, count) = idx1 + offset1;
        coordinate_base3D[std::make_pair(timestep, label)](1, count) = idx2 + offset2;
        coordinate_base3D[std::make_pair(timestep, label)](2, count) = 1;
      }
    }
  }

  lower1 = 6;
  upper1 = 9;
  lower2 = 1;
  upper2 = 9;
  label  = 7;
  Traxel trax3(label, timestep);
  trax3.features["count"] = feature_array(1, (upper1-lower1)*(upper2-lower2));
  
  coordinate_base2D[std::make_pair(timestep, label)] = arma::mat(2, (upper1-lower1)*(upper2-lower2));
  coordinate_base3D[std::make_pair(timestep, label)] = arma::mat(3, (upper1-lower1)*(upper2-lower2));
  {
    int count = 0;
    for (long idx2 = lower2; idx2 < upper2; ++idx2) {
      for (long idx1 = lower1; idx1 < upper1; ++idx1, ++count) {
        label_image2D(idx1, idx2) = label;
        label_image3D(idx1, idx2, 1) = label;
        coordinate_base2D[std::make_pair(timestep, label)](0, count) = idx1 + offset1;
        coordinate_base2D[std::make_pair(timestep, label)](1, count) = idx2 + offset2;
      
        coordinate_base3D[std::make_pair(timestep, label)](0, count) = idx1 + offset1;
        coordinate_base3D[std::make_pair(timestep, label)](1, count) = idx2 + offset2;
        coordinate_base3D[std::make_pair(timestep, label)](2, count) = 1;
      }
    }
  }

  std::cout << "MergerResolver_extract_coordinates -- extract coordinates" << std::endl;

  extract_coordinates<2, unsigned>(coordinate_result2D, label_image2D, offset2D, trax1);
  extract_coordinates<2, unsigned>(coordinate_result2D, label_image2D, offset2D, trax2);
  extract_coordinates<2, unsigned>(coordinate_result2D, label_image2D, offset2D, trax3);

  extract_coordinates<3, unsigned>(coordinate_result3D, label_image3D, offset3D, trax1);
  extract_coordinates<3, unsigned>(coordinate_result3D, label_image3D, offset3D, trax2);
  extract_coordinates<3, unsigned>(coordinate_result3D, label_image3D, offset3D, trax3);

  std::cout << "MergerResolver_extract_coordinates -- compare results" << std::endl;

  BOOST_REQUIRE_EQUAL(coordinate_result2D->size(), coordinate_base2D.size());
  BOOST_REQUIRE_EQUAL(coordinate_result3D->size(), coordinate_base3D.size());
  

  TimestepIdCoordinateMap::const_iterator result2D_it = coordinate_result2D->begin();
  TimestepIdCoordinateMap::const_iterator base2D_it = coordinate_base2D.begin();

  for (; result2D_it != coordinate_result2D->end(); ++result2D_it, ++base2D_it) {
    BOOST_REQUIRE_EQUAL(result2D_it->first, base2D_it->first);
    BOOST_REQUIRE_EQUAL(result2D_it->second.n_rows, base2D_it->second.n_rows);
    BOOST_REQUIRE_EQUAL(result2D_it->second.n_cols, base2D_it->second.n_cols);
    BOOST_CHECK_EQUAL_COLLECTIONS(result2D_it->second.begin(),
                                  result2D_it->second.end(),
                                  base2D_it->second.begin(),
                                  base2D_it->second.end()
                                  );
  }

  TimestepIdCoordinateMap::const_iterator result3D_it = coordinate_result3D->begin();
  TimestepIdCoordinateMap::const_iterator base3D_it = coordinate_base3D.begin();

  for (; result3D_it != coordinate_result3D->end(); ++result3D_it, ++base3D_it) {
    BOOST_REQUIRE_EQUAL(result3D_it->first, base3D_it->first);
    BOOST_REQUIRE_EQUAL(result3D_it->second.n_rows, base3D_it->second.n_rows);
    BOOST_REQUIRE_EQUAL(result3D_it->second.n_cols, base3D_it->second.n_cols);
    BOOST_CHECK_EQUAL_COLLECTIONS(result3D_it->second.begin(),
                                  result3D_it->second.end(),
                                  base3D_it->second.begin(),
                                  base3D_it->second.end()
                                  );
  }
  
  
}


//BOOST_AUTO_TEST_CASE( MergerResolver_helper_functions_vector_to_mat ) {
//  // void feature_array_to_arma_mat(const feature_array& in, arma::mat& out);
//
//  float arr[] = {1, 2, 0, 3, 1, 2};
//  feature_array ft(arr, arr + sizeof(arr)/sizeof(arr[0]));
//  arma::fmat m1(6,1);
//  arma::fmat m2(3,2);
//  arma::fmat m3(6,2);
//  feature_array::iterator ftIt;
//
//  feature_array_to_arma_mat(ft, m1);
//  ftIt = ft.begin();
//  for (unsigned c = 0; c < m1.n_cols; ++c, ftIt += m1.n_rows) {
//    arma::Col<float> vc = m1.col(c);
//    BOOST_CHECK_EQUAL_COLLECTIONS(ftIt, ftIt + m1.n_rows, vc.begin(), vc.end());
//  }
//
//  feature_array_to_arma_mat(ft, m2);
//  ftIt = ft.begin();
//  for (unsigned c = 0; c < m2.n_cols; ++c, ftIt += m2.n_rows) {
//    arma::Col<float> vc = m2.col(c);
//    BOOST_CHECK_EQUAL_COLLECTIONS(ftIt, ftIt + m2.n_rows, vc.begin(), vc.end());
//  }
//
//  BOOST_CHECK_THROW(feature_array_to_arma_mat(ft, m3), std::range_error);
//}
//
//
//BOOST_AUTO_TEST_CASE( MergerResolver_helper_functions_get_centers ) {
//  // void get_centers(const arma::mat& data, const arma::Col<size_t> labels, arma::mat& centers, int k);
//
//  float arr[] = {1, 2, 0, 3, 1, 2, -1, -2, 0, 1, 1, 1};
//  size_t ass[] = {0, 2, 0, 1};
//  std::vector<std::vector<float> > c;
//  for (int i = 0; i < 3; ++i) c.push_back(std::vector<float>(3));
//  c[0][0] = 0; c[0][1] = 0; c[0][2] = 0;
//  c[1][0] = 1; c[1][1] = 1; c[1][2] = 1;
//  c[2][0] = 3; c[2][1] = 1; c[2][2] = 2;
//  int k = 3;
//  feature_array ft(arr, arr + sizeof(arr)/sizeof(arr[0]));
//  std::vector<size_t> asgn(ass, ass + sizeof(ass)/sizeof(ass[0]));
//  arma::mat data(3,4);
//  arma::mat centers(3,3);
//  feature_array_to_arma_mat(ft, data);
//  arma::Col<size_t> assign(asgn);
//  get_centers(data, assign, centers, k);
//  for (unsigned n = 0; n < centers.n_cols; ++n) {
//    arma::Col<double> fvec = centers.col(n);
//    BOOST_CHECK_EQUAL_COLLECTIONS(c[n].begin(), c[n].end(), fvec.begin(), fvec.end());
//  }
//}


/* BOOST_AUTO_TEST_CASE( Merger_Resolver_GMM ) {
  float arr[] = {-6, -5, -4, 6, 5, 4};
  float arr_res[] = {-5, 5};
  feature_array data(arr, arr + sizeof(arr)/sizeof(arr[0]));
  GMM gmm(2, 1, data);
  feature_array centers = gmm();
  cout << "GMM score: " << gmm.score() << "\n";
  BOOST_CHECK_EQUAL_COLLECTIONS(centers.begin(), centers.end(), arr_res, arr_res+sizeof(arr_res)/sizeof(arr_res[0]));
} */

// EOF
