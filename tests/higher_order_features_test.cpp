#define BOOST_TEST_MODULE higher_order_features_auxiliary_test

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "pgmlink/log.h"
#include "pgmlink/higher_order_features.h"

using namespace pgmlink;
// using namespace boost;
typedef typename
  property_map<node_timestep, HypothesesGraph::base_graph>::type
  node_timestep_type;
typedef typename
  property_map<node_active_count, HypothesesGraph::base_graph>::type
  node_active_map_type;
typedef typename
  property_map<arc_active_count, HypothesesGraph::base_graph>::type
  arc_active_map_type;
typedef typename
  property_map<division_active_count, HypothesesGraph::base_graph>::type
  div_active_map_type;
typedef typename
  property_map<node_traxel, HypothesesGraph::base_graph>::type
  node_traxel_type;
typedef typename
  property_map<node_tracklet, HypothesesGraph::base_graph>::type
  node_tracklet_type;
typedef typename HypothesesGraph::NodeIt NodeIt;
typedef typename HypothesesGraph::InArcIt InArcIt;
typedef typename HypothesesGraph::OutArcIt OutArcIt;

void get_graph(HypothesesGraph& graph) {
  HypothesesGraph::Node n00 = graph.add_node(1);
  HypothesesGraph::Node n01 = graph.add_node(1);
  HypothesesGraph::Node n10 = graph.add_node(2);
  HypothesesGraph::Node n11 = graph.add_node(2);
  HypothesesGraph::Node n20 = graph.add_node(3);
  HypothesesGraph::Node n21 = graph.add_node(3);

  HypothesesGraph::Arc a13 = graph.addArc(n00, n10);
  HypothesesGraph::Arc a35 = graph.addArc(n10, n20);
  HypothesesGraph::Arc a36 = graph.addArc(n10, n21);
  HypothesesGraph::Arc a24 = graph.addArc(n01, n11);

  // corresponding traxels
  Traxel t00, t01, t10, t11, t20, t21;
  t00.Id = 1;
  t00.Timestep = 1;
  t00.features["com"].push_back(0.);
  t00.features["com"].push_back(0.);
  t00.features["id"].push_back(1.);
  t01.Id = 2;
  t01.Timestep = 1;
  t01.features["com"].push_back(0.);
  t01.features["com"].push_back(1.);
  t01.features["id"].push_back(2.);
  t10.Id = 3;
  t10.Timestep = 2;
  t10.features["com"].push_back(1.);
  t10.features["com"].push_back(0.);
  t10.features["id"].push_back(3.);
  t11.Id = 4;
  t11.Timestep = 2;
  t11.features["com"].push_back(1.);
  t11.features["com"].push_back(1.);
  t11.features["id"].push_back(4.);
  t20.Id = 5;
  t20.Timestep = 3;
  t20.features["com"].push_back(2.);
  t20.features["com"].push_back(0.);
  t20.features["id"].push_back(5.);
  t21.Id = 6;
  t21.Timestep = 3;
  t21.features["com"].push_back(2.);
  t21.features["com"].push_back(1.);
  t21.features["id"].push_back(6.);

  graph.add(node_traxel());
  graph.get(node_traxel()).set(n00, t00);
  graph.get(node_traxel()).set(n01, t01);
  graph.get(node_traxel()).set(n10, t10);
  graph.get(node_traxel()).set(n11, t11);
  graph.get(node_traxel()).set(n20, t20);
  graph.get(node_traxel()).set(n21, t21);

  graph.add(node_active_count());
  graph.get(node_active_count()).set(n00, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n01, std::vector<size_t>(1,0));
  graph.get(node_active_count()).set(n10, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n11, std::vector<size_t>(1,0));
  graph.get(node_active_count()).set(n20, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n21, std::vector<size_t>(1,1));

  graph.add(arc_active_count());
  graph.get(arc_active_count()).set(a13, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a35, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a36, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a24, std::vector<bool>(1,false));

  graph.add(division_active_count());
  graph.get(division_active_count()).set(n00, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n01, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n10, std::vector<bool>(1,true));
  graph.get(division_active_count()).set(n11, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n20, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n21, std::vector<bool>(1,false));
}

void get_tracklet_graph(HypothesesGraph& graph) {
  HypothesesGraph::Node n00 = graph.add_node(1);
  HypothesesGraph::Node n01 = graph.add_node(1);
  HypothesesGraph::Node n10 = graph.add_node(2);
  HypothesesGraph::Node n11 = graph.add_node(2);

  HypothesesGraph::Arc a13 = graph.addArc(n00, n10);
  HypothesesGraph::Arc a14 = graph.addArc(n00, n11);

  // corresponding traxels
  Traxel t00, t01, t10, t11, t20, t21;
  t00.Id = 1;
  t00.Timestep = 1;
  t00.features["com"].push_back(0.);
  t00.features["com"].push_back(0.);
  t00.features["id"].push_back(1.);
  t01.Id = 2;
  t01.Timestep = 1;
  t01.features["com"].push_back(0.);
  t01.features["com"].push_back(1.);
  t01.features["id"].push_back(2.);
  t10.Id = 3;
  t10.Timestep = 2;
  t10.features["com"].push_back(1.);
  t10.features["com"].push_back(0.);
  t10.features["id"].push_back(3.);
  t11.Id = 4;
  t11.Timestep = 2;
  t11.features["com"].push_back(1.);
  t11.features["com"].push_back(1.);
  t11.features["id"].push_back(4.);
  t20.Id = 5;
  t20.Timestep = 3;
  t20.features["com"].push_back(2.);
  t20.features["com"].push_back(0.);
  t20.features["id"].push_back(5.);
  t21.Id = 6;
  t21.Timestep = 3;
  t21.features["com"].push_back(2.);
  t21.features["com"].push_back(1.);
  t21.features["id"].push_back(6.);

  graph.add(node_tracklet());
  std::vector<Traxel> t1; t1.push_back(t00); t1.push_back(t10);
  graph.get(node_tracklet()).set(n00, t1);
  std::vector<Traxel> t2; t2.push_back(t01); t2.push_back(t11);
  graph.get(node_tracklet()).set(n01, t2);
  graph.get(node_tracklet()).set(n10, std::vector<Traxel>(1,t20));
  graph.get(node_tracklet()).set(n11, std::vector<Traxel>(1,t21));

  graph.add(node_active_count());
  graph.get(node_active_count()).set(n00, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n01, std::vector<size_t>(1,0));
  graph.get(node_active_count()).set(n10, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n11, std::vector<size_t>(1,1));

  graph.add(arc_active_count());
  graph.get(arc_active_count()).set(a13, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a14, std::vector<bool>(1,true));

  graph.add(division_active_count());
  graph.get(division_active_count()).set(n00, std::vector<bool>(1,true));
  graph.get(division_active_count()).set(n01, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n10, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n11, std::vector<bool>(1,false));
}

BOOST_AUTO_TEST_CASE( SubsetFeaturesIdentity_extract_matrix_traxelgraph ) {
  LOG(logINFO) << "test case: SubsetFeaturesIdentity_operator_traxelgraph";

  // set up the graph
  HypothesesGraph graph;
  get_graph(graph);

  // get the subset of all traxels in the order of their traxel ids
  ConstTraxelRefVector subset(6);
  for (NodeIt n_it(graph); n_it != lemon::INVALID; ++n_it) {
    const Traxel* tref = &(graph.get(node_traxel())[n_it]);
    subset[tref->Id-1] = tref;
  }

  set_solution(graph, 0);

  LOG(logINFO) << "  test \"SubsetFeaturesIdentity\" with string as constructor argument";
  {
    SubsetFeaturesIdentity identity("com");

    FeatureMatrix com_matrix = identity.extract_matrix(subset);
    
    BOOST_CHECK_EQUAL(com_matrix.shape(0), 6);
    BOOST_CHECK_EQUAL(com_matrix.shape(1), 2);
    for (size_t i = 0; i != 3; ++i) {
      size_t x = 2*i;
      BOOST_CHECK_EQUAL(com_matrix(x  ,0), static_cast<feature_type>(i));
      BOOST_CHECK_EQUAL(com_matrix(x  ,1), 0);
      BOOST_CHECK_EQUAL(com_matrix(x+1,0), static_cast<feature_type>(i));
      BOOST_CHECK_EQUAL(com_matrix(x+1,1), 1);
    }
  }

  LOG(logINFO) << "  test \"SubsetFeaturesIdentity\" with length one vector as constructor argument";
  {
    std::vector<std::string> feature_names(1, "com");
    SubsetFeaturesIdentity identity(feature_names);

    FeatureMatrix com_matrix = identity.extract_matrix(subset);
    
    BOOST_CHECK_EQUAL(com_matrix.shape(0), 6);
    BOOST_CHECK_EQUAL(com_matrix.shape(1), 2);
    for (size_t i = 0; i != 3; ++i) {
      size_t x = 2*i;
      BOOST_CHECK_EQUAL(com_matrix(x  ,0), static_cast<feature_type>(i));
      BOOST_CHECK_EQUAL(com_matrix(x  ,1), 0);
      BOOST_CHECK_EQUAL(com_matrix(x+1,0), static_cast<feature_type>(i));
      BOOST_CHECK_EQUAL(com_matrix(x+1,1), 1);
    }
  }

  LOG(logINFO) << "  test \"SubsetFeaturesIdentity\" with length two vector as constructor argument";
  {
    std::vector<std::string> feature_names;
    feature_names.push_back("com");
    feature_names.push_back("id");

    SubsetFeaturesIdentity identity(feature_names);

    FeatureMatrix com_matrix = identity.extract_matrix(subset);
    
    BOOST_CHECK_EQUAL(com_matrix.shape(0), 6);
    BOOST_CHECK_EQUAL(com_matrix.shape(1), 3);
    for (size_t i = 0; i != 3; ++i) {
      size_t x = 2*i;
      BOOST_CHECK_EQUAL(com_matrix(x  ,0), static_cast<feature_type>(i));
      BOOST_CHECK_EQUAL(com_matrix(x  ,1), 0);
      BOOST_CHECK_EQUAL(com_matrix(x  ,2), static_cast<feature_type>(2*i+1));
      BOOST_CHECK_EQUAL(com_matrix(x+1,0), static_cast<feature_type>(i));
      BOOST_CHECK_EQUAL(com_matrix(x+1,1), 1);
      BOOST_CHECK_EQUAL(com_matrix(x+1,2), static_cast<feature_type>(2*i+2));
    }
  }
}

BOOST_AUTO_TEST_CASE( TrackSubsets_operator_traxelgraph ) {
  
  LOG(logINFO) << "test case: SubsetFeaturesIdentity_operator";

  // set up the graph
  HypothesesGraph graph;
  get_graph(graph);

  // set the solution index
  set_solution(graph, 0);
  
  // get the track subsets
  TrackSubsets get_track_subsets;
  std::vector<ConstTraxelRefVector> track_subsets = get_track_subsets(graph);

  LOG(logINFO) << "  there are " << track_subsets.size() << " tracks";
  for (
    std::vector<ConstTraxelRefVector>::iterator tvec_it = track_subsets.begin();
    tvec_it != track_subsets.end();
    tvec_it++
  ) {
    LOG(logINFO) << "    count of nodes in track: " << tvec_it->size();
    std::stringstream sstream;
    sstream << "    the ids are: ";
    for (
      ConstTraxelRefVector::iterator tref_it = tvec_it->begin();
      tref_it != tvec_it->end();
      tref_it++
    ) {
      sstream << (*tref_it)->Id << " ";
    }
    LOG(logINFO) << sstream.str();
  }
  LOG(logINFO) << "  test sizes";
  BOOST_CHECK_EQUAL(track_subsets.size(), 3);
  BOOST_CHECK_EQUAL(track_subsets[0].size(), 2);
  BOOST_CHECK_EQUAL(track_subsets[1].size(), 1);
  BOOST_CHECK_EQUAL(track_subsets[2].size(), 1);
}

BOOST_AUTO_TEST_CASE( MVNOutlierCalculator_calculate ) {
  LOG(logINFO) << "test case: MVNOutlierCalculator_calculate";
  // Set up the test data with outlier x6
  feature_type x_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {9., 8.}
  };
  feature_arrays x;
  for(size_t i = 0; i < 7; i++) {
    feature_array y(x_array[i], x_array[i]+2);
    x.push_back(y);
  }

  // Create outlier detection
  LOG(logINFO) << "  set up the outlier calculator";
  MVNOutlierCalculator mvnoutlier;
  LOG(logINFO) << "  test outliers";
  std::vector<size_t> outlier(mvnoutlier.calculate(x));
  BOOST_CHECK_EQUAL(outlier.size(), 1);
  BOOST_CHECK_EQUAL(outlier[0], 6);
}

BOOST_AUTO_TEST_CASE( OutlierCountAggregator_calculate ) {
  LOG(logINFO) << "test case: OutlierCountAggregator_calculate";
  // Set up the test data with outlier x6
  feature_type x_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {9., 8.}
  };
  feature_arrays x;
  for(size_t i = 0; i < 7; i++) {
    feature_array y(x_array[i], x_array[i]+2);
    x.push_back(y);
  }

  // Create outlier detection
  LOG(logINFO) << "  Set up the outlier count aggregator";
  OutlierCountAggregator outliercount;
  LOG(logINFO) << "  Calculate the outlier count";
  size_t outlier = outliercount(x);
  LOG(logINFO) << "  outlier count: " << outlier;

  BOOST_CHECK_EQUAL(outlier, 1);
}

BOOST_AUTO_TEST_CASE( OutlierBadnessAggregator_calculate ) {
  LOG(logINFO) << "test case: OutlierBadnessAggregator_calculate";
  // Set up the test data with outlier x6
  feature_type x_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {9., 8.}
  };
  feature_arrays x;
  for(size_t i = 0; i < 7; i++) {
    feature_array y(x_array[i], x_array[i]+2);
    x.push_back(y);
  }

  // Create outlier detection
  LOG(logINFO) << "  Set up the outlier badness aggregator";
  OutlierBadnessAggregator outlierbadness;
  LOG(logINFO) << "  Calculate the outlier badness";
  feature_type value = outlierbadness(x);
  BOOST_CHECK( value > 3 );
}
