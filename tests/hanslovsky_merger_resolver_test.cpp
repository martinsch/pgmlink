#define BOOST_TEST_MODULE hanslovsky_merger_resolver_test

#include <stdexcept>
#include <cstring>
#include <iostream>
#include <algorithm>

#include <boost/test/unit_test.hpp>

#include "pgmlink/hypotheses.h"
#include "pgmlink/traxels.h"
// enable tests of private class members
#define private public
#include "pgmlink/hanslovsky.h"
#undef private

#include <armadillo>
#include <mlpack/core.hpp>
#include <mlpack/methods/kmeans/kmeans.hpp>
#include <lemon/maps.h>

using namespace pgmlink;
using namespace std;
using namespace boost;


BOOST_AUTO_TEST_CASE( MergerResolver_constructor ) {
  HypothesesGraph g;
  MergerResolver m(&g);
  BOOST_CHECK_EQUAL(m.g_, &g);
}


BOOST_AUTO_TEST_CASE( MergerResolver_resolve_mergers ) {
  HypothesesGraph g;
  MergerResolver m(&g);
}


BOOST_AUTO_TEST_CASE( MergerResolver_refine_node) {
  // MergerResolver::refine_node(HypothesesGraph::Node node, std::size_t nMerger)

  //  t=1      2      3 
  //    o ----   --- o
  //          \ /
  //           O
  //          / \
  //    o ----   --- o

  // ->
  //    o ---- o ---- o
  //     \    / \    /
  //       --     --
  //     /    \ /    \
  //    o ---- o ---- o

  
  HypothesesGraph g;
  // g.add(arc_distance()).add(tracklet_intern_dist()).add(node_tracklet()).add(tracklet_intern_arc_ids()).add(traxel_arc_id());
  g.add(node_traxel()).add(arc_distance()).add(arc_active()).add(node_active2());
  
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
  property_map<node_timestep, HypothesesGraph::base_graph>::type::ItemIt timeIt(timestep_map, 2);

  MergerResolver m(&g);
  m.refine_node(n21, 2);
  
  // deactivated arcs from and to merger node
  BOOST_CHECK_EQUAL(arc_map[a11_21], false);
  BOOST_CHECK_EQUAL(arc_map[a12_21], false);
  BOOST_CHECK_EQUAL(arc_map[a21_31], false);
  BOOST_CHECK_EQUAL(arc_map[a21_32], false);

  int count = 0;
  for (; timeIt != lemon::INVALID; ++timeIt) {
    if (timeIt == n21) {
      // deactivated merger node
      BOOST_CHECK_EQUAL(active_map[timeIt], 2);
    } else {

      // if correct merger COM exists in traxel, count++
      Traxel trax = traxel_map[timeIt];
      com = trax.features["com"];
      if (std::equal(com1.begin(), com1.end(), com.begin()) ||
	  std::equal(com2.begin(), com2.end(), com.begin())) {
	++count;
      }
      // activated merger replacement nodes
      BOOST_CHECK_EQUAL(active_map[timeIt], 1);
    }
  }
  BOOST_CHECK_EQUAL(count, 2);
  
}


BOOST_AUTO_TEST_CASE( MergerResolver_deactivate_arcs ) {
  // deactivate_arcs(std::vector<HypothesesGraph::base_graph::Arc>)
  
  //  t=1      2      3 
  //    o ---- o ---- o
             
  // ->
  //    o      o      o

  
  HypothesesGraph g;
  g.add(arc_active());
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
  // deactivate_nodes(std::vector<HypothesesGraph::Node> nodes)

  // (1) -> (0)

  
  HypothesesGraph g;
  g.add(node_active2());
  HypothesesGraph::Node n = g.add_node(1);
  g.get(node_active2()).set(n, 1);
  
  // node is active with 1 object
  BOOST_CHECK_EQUAL(g.get(node_active2())[n], 1);

  std::vector<HypothesesGraph::Node> nodes;
  nodes.push_back(n);
  MergerResolver m(&g);
  m.deactivate_nodes(nodes);
  
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
  g.add(node_traxel());
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
  // MergerResolver::add_arcs_for_replacement_node(HypothesesGraph::Node node,
  //						   Traxel trax,
  // 						   std::vector<HypothesesGraph::base_graph::Arc> src,
  //						   std::vector<HypothesesGraph::base_graph::Arc> dest);
  
  //  t=1      2      3 
  //    o ---- o ---- o
             
  // ->
  //    o      o      o
  //     \           /
  //      ---- o ----


  feature_array COM(3, 0.0);
  TraxelStore ts;
  Traxel t1, t2, t3, t22;
  
  t1.Id = 11;
  t1.Timestep = 1;
  t1.features["com"] = COM;
  
  t2.Id = 21;
  t2.Timestep = 2;
  t2.features["com"] = COM;
  
  t3.Id = 31;
  t3.Timestep = 3;
  t3.features["com"] = COM;

  COM[0] = 100;
  t22.Id = 22;
  t22.Timestep = 2;
  t22.features["com"] = COM;
  
  add(ts, t1);
  add(ts, t2);
  add(ts, t3);
  SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(1, // max NN
							       1, // max Dist
							       false, // forward_backward
							       false, //  consider_divisions
							       100 // division_threshold
							       );
  SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
  HypothesesGraph* g = hyp_builder.build();
  HypothesesGraph& G = *g;
  G.add(arc_distance()).add(arc_active());
  MergerResolver m(g);

  property_map<node_timestep, HypothesesGraph::base_graph>::type& time_map = g->get(node_timestep());
  property_map<node_timestep, HypothesesGraph::base_graph>::type::ItemIt time_it(time_map, 2);
  property_map<node_traxel, HypothesesGraph::base_graph>::type& trax_map = g->get(node_traxel());
  HypothesesGraph::Node newNode = g->add_node(t22.Timestep);
  trax_map.set(newNode, t22);
  for (; time_it != lemon::INVALID; ++time_it) {
    HypothesesGraph::Node node = time_it;
    std::vector<HypothesesGraph::base_graph::Arc> sources;
    std::vector<HypothesesGraph::base_graph::Arc> targets;
    m.collect_arcs(HypothesesGraph::base_graph::InArcIt(*g, node), sources);
    m.collect_arcs(HypothesesGraph::base_graph::OutArcIt(*g, node), targets);
    if (sources.size() > 0 && targets.size() > 0) {
      m.add_arcs_for_replacement_node(newNode, t22, sources, targets);
    }
  }
  property_map<arc_distance, HypothesesGraph::base_graph>::type& distance_map = g->get(arc_distance());
  property_map<arc_active, HypothesesGraph::base_graph>::type& active_map = g->get(arc_active());
  HypothesesGraph::base_graph::InArcIt inIt(G, newNode);
  int k = 0;
  for (; inIt != lemon::INVALID; ++k, ++inIt) {
    BOOST_CHECK_EQUAL(distance_map[inIt], 100);
    BOOST_CHECK(active_map[inIt]);
  }
  BOOST_CHECK_EQUAL(k, 1);
  
  HypothesesGraph::base_graph::OutArcIt outIt(G, newNode);
  k = 0;
  for (; outIt != lemon::INVALID; ++k, ++outIt) {
    BOOST_CHECK_EQUAL(distance_map[outIt], 100);
    BOOST_CHECK(active_map[outIt]);
  }
  BOOST_CHECK_EQUAL(k, 1);
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
  MergerResolver m(g);
  
  std::vector<HypothesesGraph::base_graph::Arc> sources;
  std::vector<HypothesesGraph::base_graph::Arc> targets;
  
  property_map<node_timestep, HypothesesGraph::base_graph>::type& time_map = g->get(node_timestep());
  property_map<node_timestep, HypothesesGraph::base_graph>::type::ItemIt time_it(time_map, 2);
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


BOOST_AUTO_TEST_CASE( MergerResolver_helper_functions_vector_to_mat ) {
  // void feature_array_to_arma_mat(const feature_array& in, arma::mat& out);
  
  float arr[] = {1, 2, 0, 3, 1, 2};
  feature_array ft(arr, arr + sizeof(arr)/sizeof(arr[0]));
  arma::fmat m1(6,1);
  arma::fmat m2(3,2);
  arma::fmat m3(6,2);
  feature_array::iterator ftIt;

  feature_array_to_arma_mat(ft, m1);
  ftIt = ft.begin();
  for (unsigned c = 0; c < m1.n_cols; ++c, ftIt += m1.n_rows) {
    arma::Col<float> vc = m1.col(c);
    BOOST_CHECK_EQUAL_COLLECTIONS(ftIt, ftIt + m1.n_rows, vc.begin(), vc.end());
  }
				
  feature_array_to_arma_mat(ft, m2);
  ftIt = ft.begin();
  for (unsigned c = 0; c < m2.n_cols; ++c, ftIt += m2.n_rows) {
    arma::Col<float> vc = m2.col(c);
    BOOST_CHECK_EQUAL_COLLECTIONS(ftIt, ftIt + m2.n_rows, vc.begin(), vc.end());
  }

  BOOST_CHECK_THROW(feature_array_to_arma_mat(ft, m3), std::range_error);
}


BOOST_AUTO_TEST_CASE( MergerResolver_helper_functions_get_centers ) {
  // void get_centers(const arma::mat& data, const arma::Col<size_t> labels, arma::mat& centers, int k);
  
  float arr[] = {1, 2, 0, 3, 1, 2, -1, -2, 0, 1, 1, 1};
  size_t ass[] = {0, 2, 0, 1};
  std::vector<std::vector<float> > c;
  for (int i = 0; i < 3; ++i) c.push_back(std::vector<float>(3));
  c[0][0] = 0; c[0][1] = 0; c[0][2] = 0;
  c[1][0] = 1; c[1][1] = 1; c[1][2] = 1;
  c[2][0] = 3; c[2][1] = 1; c[2][2] = 2;
  int k = 3;
  feature_array ft(arr, arr + sizeof(arr)/sizeof(arr[0]));
  std::vector<size_t> asgn(ass, ass + sizeof(ass)/sizeof(ass[0]));
  arma::mat data(3,4);
  arma::mat centers(3,3);
  feature_array_to_arma_mat(ft, data);
  arma::Col<size_t> assign(asgn);
  get_centers(data, assign, centers, k);
  for (unsigned n = 0; n < centers.n_cols; ++n) {
    arma::Col<double> fvec = centers.col(n);
    BOOST_CHECK_EQUAL_COLLECTIONS(c[n].begin(), c[n].end(), fvec.begin(), fvec.end());
  }
}


// EOF
