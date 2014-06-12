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
  HypothesesGraph::Node n30 = graph.add_node(4);
  HypothesesGraph::Node n31 = graph.add_node(4);

  HypothesesGraph::Arc a13 = graph.addArc(n00, n10);
  HypothesesGraph::Arc a24 = graph.addArc(n01, n11);
  HypothesesGraph::Arc a35 = graph.addArc(n10, n20);
  HypothesesGraph::Arc a36 = graph.addArc(n10, n21);
  HypothesesGraph::Arc a57 = graph.addArc(n20, n30);
  HypothesesGraph::Arc a68 = graph.addArc(n21, n31);


  // corresponding traxels
  Traxel t00, t01, t10, t11, t20, t21, t30, t31;
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
  t30.Id = 7;
  t30.Timestep = 4;
  t30.features["com"].push_back(3.);
  t30.features["com"].push_back(0.);
  t30.features["id"].push_back(7.);
  t31.Id = 8;
  t31.Timestep = 4;
  t31.features["com"].push_back(3.);
  t31.features["com"].push_back(1.);
  t31.features["id"].push_back(8.);

  graph.add(node_traxel());
  graph.get(node_traxel()).set(n00, t00);
  graph.get(node_traxel()).set(n01, t01);
  graph.get(node_traxel()).set(n10, t10);
  graph.get(node_traxel()).set(n11, t11);
  graph.get(node_traxel()).set(n20, t20);
  graph.get(node_traxel()).set(n21, t21);
  graph.get(node_traxel()).set(n30, t30);
  graph.get(node_traxel()).set(n31, t31);

  graph.add(node_active_count());
  graph.get(node_active_count()).set(n00, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n01, std::vector<size_t>(1,0));
  graph.get(node_active_count()).set(n10, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n11, std::vector<size_t>(1,0));
  graph.get(node_active_count()).set(n20, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n21, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n30, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n31, std::vector<size_t>(1,1));

  graph.add(arc_active_count());
  graph.get(arc_active_count()).set(a13, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a24, std::vector<bool>(1,false));
  graph.get(arc_active_count()).set(a35, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a36, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a57, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a68, std::vector<bool>(1,true));

  graph.add(division_active_count());
  graph.get(division_active_count()).set(n00, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n01, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n10, std::vector<bool>(1,true));
  graph.get(division_active_count()).set(n11, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n20, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n21, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n30, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n31, std::vector<bool>(1,false));
}

void get_tracklet_graph(HypothesesGraph& graph) {
  HypothesesGraph::Node n00 = graph.add_node(1);
  HypothesesGraph::Node n01 = graph.add_node(1);
  HypothesesGraph::Node n10 = graph.add_node(3);
  HypothesesGraph::Node n11 = graph.add_node(3);
  HypothesesGraph::Node n21 = graph.add_node(4);

  HypothesesGraph::Arc a13 = graph.addArc(n00, n10);
  HypothesesGraph::Arc a14 = graph.addArc(n00, n11);
  HypothesesGraph::Arc a45 = graph.addArc(n11, n21);

  // corresponding traxels
  Traxel t00, t01, t10, t11, t20, t21, t30, t31;
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
  t30.Id = 7;
  t30.Timestep = 4;
  t30.features["com"].push_back(3.);
  t30.features["com"].push_back(0.);
  t30.features["id"].push_back(7.);
  t31.Id = 8;
  t31.Timestep = 4;
  t31.features["com"].push_back(3.);
  t31.features["com"].push_back(1.);
  t31.features["id"].push_back(8.);

  graph.add(node_tracklet());
  std::vector<Traxel> t1; t1.push_back(t00); t1.push_back(t10);
  graph.get(node_tracklet()).set(n00, t1);
  std::vector<Traxel> t2; t2.push_back(t01); t2.push_back(t11);
  graph.get(node_tracklet()).set(n01, t2);
  std::vector<Traxel> t3; t3.push_back(t20); t3.push_back(t30);
  graph.get(node_tracklet()).set(n10, t3);
  graph.get(node_tracklet()).set(n11, std::vector<Traxel>(1,t21));
  graph.get(node_tracklet()).set(n21, std::vector<Traxel>(1,t31));

  graph.add(node_active_count());
  graph.get(node_active_count()).set(n00, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n01, std::vector<size_t>(1,0));
  graph.get(node_active_count()).set(n10, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n11, std::vector<size_t>(1,1));
  graph.get(node_active_count()).set(n21, std::vector<size_t>(1,1));

  graph.add(arc_active_count());
  graph.get(arc_active_count()).set(a13, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a14, std::vector<bool>(1,true));
  graph.get(arc_active_count()).set(a45, std::vector<bool>(1,true));

  graph.add(division_active_count());
  graph.get(division_active_count()).set(n00, std::vector<bool>(1,true));
  graph.get(division_active_count()).set(n01, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n10, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n11, std::vector<bool>(1,false));
  graph.get(division_active_count()).set(n21, std::vector<bool>(1,false));
}

void get_feature_matrix(FeatureMatrix& x) {
  feature_type x_array[7][2] = {
    {3., 3.},
    {4., 3.},
    {3., 4.},
    {4., 4.},
    {5., 4.},
    {5., 5.},
    {9., 8.}
  };

  x.reshape(vigra::Shape2(7,2));
  for(size_t i = 0; i < 7; i++) {
    for(size_t j = 0; j < 2; j++) {
      x(i, j) = x_array[i][j];
    }
  }
}

BOOST_AUTO_TEST_CASE( SubsetFeaturesIdentity_extract ) {
  LOG(logINFO) << "test case: SubsetFeaturesIdentity_operator";

  // set up the graph
  HypothesesGraph graph;
  get_graph(graph);

  // get the subset of all traxels in the order of their traxel ids
  ConstTraxelRefVector subset(8);
  for (NodeIt n_it(graph); n_it != lemon::INVALID; ++n_it) {
    const Traxel* tref = &(graph.get(node_traxel())[n_it]);
    subset[tref->Id-1] = tref;
  }

  LOG(logINFO) << "  test \"SubsetFeaturesIdentity\" with string as constructor argument";
  {
    SubsetFeaturesIdentity identity("com");

    FeatureMatrix com_matrix;
    identity.extract(subset, com_matrix);
    
    BOOST_CHECK_EQUAL(com_matrix.shape(0), 8);
    BOOST_CHECK_EQUAL(com_matrix.shape(1), 2);
    for (size_t i = 0; i != 4; ++i) {
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

    FeatureMatrix com_matrix;
    identity.extract(subset, com_matrix);
    
    BOOST_CHECK_EQUAL(com_matrix.shape(0), 8);
    BOOST_CHECK_EQUAL(com_matrix.shape(1), 2);
    for (size_t i = 0; i != 4; ++i) {
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

    FeatureMatrix com_matrix;
    identity.extract(subset, com_matrix);
    
    BOOST_CHECK_EQUAL(com_matrix.shape(0), 8);
    BOOST_CHECK_EQUAL(com_matrix.shape(1), 3);
    for (size_t i = 0; i != 4; ++i) {
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
  LOG(logINFO) << "test case: TrackSubsets_operator_traxelgraph";

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
  BOOST_CHECK_EQUAL(track_subsets[1].size(), 2);
  BOOST_CHECK_EQUAL(track_subsets[2].size(), 2);
}

BOOST_AUTO_TEST_CASE( TrackSubsets_operator_trackletgraph ) {
  LOG(logINFO) << "test case: TrackSubsets_operator_trackletgraph";

  // set up the graph
  HypothesesGraph graph;
  get_tracklet_graph(graph);

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
  BOOST_CHECK_EQUAL(track_subsets[1].size(), 2);
  BOOST_CHECK_EQUAL(track_subsets[2].size(), 2);
}

BOOST_AUTO_TEST_CASE( DivisionSubsets_operator_traxelgraph ) {
  LOG(logINFO) << "test case: DivisionSubsets_operator_traxelgraph";

  // set up the graph
  HypothesesGraph graph;
  get_graph(graph);

  // set the solution index
  set_solution(graph, 0);

  // get the division subsets to depth 1
  LOG(logINFO) << "  get division subsets to depth 1";
  {
    DivisionSubsets get_div_subsets;
    std::vector<ConstTraxelRefVector> div_subsets = get_div_subsets(graph);

    LOG(logINFO) << "  count of division subsets: " << div_subsets.size();
    LOG(logINFO) << "  check size and indices";
    BOOST_CHECK_EQUAL(div_subsets.size(), 1);
    ConstTraxelRefVector traxelrefs = div_subsets[0];
    BOOST_CHECK_EQUAL(traxelrefs.size(), 3);
    BOOST_CHECK_EQUAL(traxelrefs[0]->Id, 3);
    BOOST_CHECK((traxelrefs[1]->Id == 5) or (traxelrefs[1]->Id == 6));
    BOOST_CHECK((traxelrefs[2]->Id == 5) or (traxelrefs[2]->Id == 6));
  }

  // get the division subsets to depth 2
  LOG(logINFO) << "  get division subsets to depth 2";
  {
    DivisionSubsets get_div_subsets(2);
    std::vector<ConstTraxelRefVector> div_subsets = get_div_subsets(graph);

    LOG(logINFO) << "  count of division subsets: " << div_subsets.size();
    LOG(logINFO) << "  check size and indices";
    BOOST_CHECK_EQUAL(div_subsets.size(), 1);

    ConstTraxelRefVector traxelrefs = div_subsets[0];
    std::stringstream sstream;
    sstream << "    indices in division: ";
    for (
      ConstTraxelRefVector::iterator tref = traxelrefs.begin();
      tref != traxelrefs.end();
      tref++
    ) {
      sstream << (*tref)->Id << " ";
    }
    LOG(logINFO) << sstream.str();
    BOOST_CHECK_EQUAL(traxelrefs.size(), 6);
    BOOST_CHECK_EQUAL(traxelrefs[0]->Id, 3);
    BOOST_CHECK_EQUAL(traxelrefs[1]->Id, 1);
    BOOST_CHECK((traxelrefs[2]->Id == 5) or (traxelrefs[2]->Id == 6));
    BOOST_CHECK_EQUAL(traxelrefs[3]->Id - traxelrefs[2]->Id, 2);
    BOOST_CHECK((traxelrefs[4]->Id == 5) or (traxelrefs[4]->Id == 6));
    BOOST_CHECK_EQUAL(traxelrefs[5]->Id - traxelrefs[4]->Id, 2);
  }
}

BOOST_AUTO_TEST_CASE( DivisionSubsets_operator_trackletgraph ) {
  LOG(logINFO) << "test case: DivisionSubsets_operator_trackletgraph";

  // set up the graph
  HypothesesGraph graph;
  get_tracklet_graph(graph);

  // set the solution index
  set_solution(graph, 0);

  // get the division subsets to depth 1
  LOG(logINFO) << "  get division subsets to depth 1";
  {
    DivisionSubsets get_div_subsets;
    std::vector<ConstTraxelRefVector> div_subsets = get_div_subsets(graph);

    LOG(logINFO) << "  count of division subsets: " << div_subsets.size();
    LOG(logINFO) << "  check size and indices";
    BOOST_CHECK_EQUAL(div_subsets.size(), 1);
    ConstTraxelRefVector traxelrefs = div_subsets[0];
    BOOST_CHECK_EQUAL(traxelrefs.size(), 3);
    BOOST_CHECK_EQUAL(traxelrefs[0]->Id, 3);
    BOOST_CHECK((traxelrefs[1]->Id == 5) or (traxelrefs[1]->Id == 6));
    BOOST_CHECK((traxelrefs[2]->Id == 5) or (traxelrefs[2]->Id == 6));
  }

  // get the division subsets to depth 2
  LOG(logINFO) << "  get division subsets to depth 2";
  {
    DivisionSubsets get_div_subsets(2);
    std::vector<ConstTraxelRefVector> div_subsets = get_div_subsets(graph);

    LOG(logINFO) << "  count of division subsets: " << div_subsets.size();
    LOG(logINFO) << "  check size and indices";
    BOOST_CHECK_EQUAL(div_subsets.size(), 1);

    ConstTraxelRefVector traxelrefs = div_subsets[0];
    std::stringstream sstream;
    sstream << "    indices in division: ";
    for (
      ConstTraxelRefVector::iterator tref = traxelrefs.begin();
      tref != traxelrefs.end();
      tref++
    ) {
      sstream << (*tref)->Id << " ";
    }
    LOG(logINFO) << sstream.str();
    BOOST_CHECK_EQUAL(traxelrefs.size(), 6);
    BOOST_CHECK_EQUAL(traxelrefs[0]->Id, 3);
    BOOST_CHECK_EQUAL(traxelrefs[1]->Id, 1);
    BOOST_CHECK((traxelrefs[2]->Id == 5) or (traxelrefs[2]->Id == 6));
    BOOST_CHECK_EQUAL(traxelrefs[3]->Id - traxelrefs[2]->Id, 2);
    BOOST_CHECK((traxelrefs[4]->Id == 5) or (traxelrefs[4]->Id == 6));
    BOOST_CHECK_EQUAL(traxelrefs[5]->Id - traxelrefs[4]->Id, 2);
  }
}

BOOST_AUTO_TEST_CASE( SumCalculator_test ) {
  LOG(logINFO) << "test case: SumCalculator_test";

  // get the test data
  FeatureMatrix x(vigra::Shape2(0,0));
  get_feature_matrix(x);

  // set up the reference data
  FeatureVector x_sum_vec_ref(vigra::Shape1(2));
  x_sum_vec_ref(0) = 33;
  x_sum_vec_ref(1) = 31;

  // set up the SumCalculator
  SumCalculator sum_calculator;

  FeatureMatrix x_sum_vec;
  sum_calculator.calculate(x, x_sum_vec);
  BOOST_CHECK_EQUAL(x_sum_vec.shape(0), 1);
  BOOST_CHECK_EQUAL(x_sum_vec.shape(1), 2);
  BOOST_CHECK_EQUAL(x_sum_vec(0, 0), x_sum_vec_ref(0));
  BOOST_CHECK_EQUAL(x_sum_vec(0, 1), x_sum_vec_ref(1));

  // with an empty feature matrix
  FeatureMatrix y(vigra::Shape2(0,0));
  FeatureMatrix y_sum_vec;
  sum_calculator.calculate(y, y_sum_vec);
  BOOST_CHECK_EQUAL(y_sum_vec.shape(0), 1);
  BOOST_CHECK_EQUAL(y_sum_vec.shape(1), 1);
}

BOOST_AUTO_TEST_CASE( DiffCalculator_test ) {
  LOG(logINFO) << "test case: DiffCalculator_test";

  // get the test data
  FeatureMatrix x;
  get_feature_matrix(x);

  // set up the reference data
  FeatureMatrix ref_matrix(vigra::Shape2(6,2));
  ref_matrix(0, 0) =  1.; ref_matrix(0, 1) =  0.;
  ref_matrix(1, 0) = -1.; ref_matrix(1, 1) =  1.;
  ref_matrix(2, 0) =  1.; ref_matrix(2, 1) =  0.;
  ref_matrix(3, 0) =  1.; ref_matrix(3, 1) =  0.;
  ref_matrix(4, 0) =  0.; ref_matrix(4, 1) =  1.;
  ref_matrix(5, 0) =  4.; ref_matrix(5, 1) =  3.;

  // set up the DiffCalculator and run the calculation
  DiffCalculator difference_calculator;
  FeatureMatrix matrix;
  difference_calculator.calculate(x, matrix);
  
  BOOST_CHECK_EQUAL(matrix.shape(0), 6);
  BOOST_CHECK_EQUAL(matrix.shape(1), 2);
  for(size_t i = 0; i < 6; i++) {
    for(size_t j = 0; j < 2; j++) {
      BOOST_CHECK_EQUAL(ref_matrix(i,j), matrix(i,j));
    }
  }

  // with an empty feature matrix
  FeatureMatrix y(vigra::Shape2(0,0));
  difference_calculator.calculate(y, matrix);
  BOOST_CHECK_EQUAL(matrix.shape(0), 1);
  BOOST_CHECK_EQUAL(matrix.shape(1), 0);
}

BOOST_AUTO_TEST_CASE( CurveCalculator_test ) {
  LOG(logINFO) << "test case: CurveCalculator_test";

  // get the test data
  FeatureMatrix x;
  get_feature_matrix(x);

  // set up the reference data
  FeatureMatrix ref_matrix(vigra::Shape2(5,2));
  ref_matrix(0, 0) = -2.; ref_matrix(0, 1) =  1.;
  ref_matrix(1, 0) =  2.; ref_matrix(1, 1) = -1.;
  ref_matrix(2, 0) =  0.; ref_matrix(2, 1) =  0.;
  ref_matrix(3, 0) = -1.; ref_matrix(3, 1) =  1.;
  ref_matrix(4, 0) =  4.; ref_matrix(4, 1) =  2.;

  // set up the CurveCalculator and run the calculation
  CurveCalculator curvature_calculator;
  FeatureMatrix matrix;
  curvature_calculator.calculate(x, matrix);

  BOOST_CHECK_EQUAL(matrix.shape(0), 5);
  BOOST_CHECK_EQUAL(matrix.shape(1), 2);
  for(size_t i = 0; i < 5; i++) {
    for(size_t j = 0; j < 2; j++) {
      BOOST_CHECK_EQUAL(ref_matrix(i,j), matrix(i,j));
    }
  }
  
  // with an empty feature matrix
  FeatureMatrix y(vigra::Shape2(0,0));
  curvature_calculator.calculate(y, matrix);
  BOOST_CHECK_EQUAL(matrix.shape(0), 1);
  BOOST_CHECK_EQUAL(matrix.shape(1), 0);
}

BOOST_AUTO_TEST_CASE( MaxCalculator_test ) {
  LOG(logINFO) << "test case: MaxCalculator_test";

  // get the test data
  FeatureMatrix x;
  get_feature_matrix(x);

  // set up the CurveCalculator and run the calculation
  MaxCalculator max_calculator;
  FeatureMatrix matrix;
  max_calculator.calculate(x, matrix);

  BOOST_CHECK_EQUAL(matrix.shape(0), 1);
  BOOST_CHECK_EQUAL(matrix.shape(1), 2);
  BOOST_CHECK_EQUAL(matrix(0, 0), 9.);
  BOOST_CHECK_EQUAL(matrix(0, 1), 8.);
}

BOOST_AUTO_TEST_CASE( ChildParentDiffCalculator_test ) {
  LOG(logINFO) << "test case: ChildParentDiffCalculator_test";
  // get the test data
  FeatureMatrix x;
  get_feature_matrix(x);
  
  // reshape
  FeatureMatrixView y = x.subarray(vigra::Shape2(0,0), vigra::Shape2(6,2));

  // set up the reference data
  FeatureMatrix ref_matrix(vigra::Shape2(2,2));
  ref_matrix(0,0) =  0.; ref_matrix(0,1) =  1.;
  ref_matrix(1,0) =  2.; ref_matrix(1,1) =  1.;

  // set up ChildParentDiffCalculator and calculate the result
  ChildParentDiffCalculator cp_diff;
  FeatureMatrix matrix;
  cp_diff.calculate(y, matrix);
  BOOST_CHECK_EQUAL(matrix.shape(0), 2);
  BOOST_CHECK_EQUAL(matrix.shape(1), 2);
  BOOST_CHECK_EQUAL(matrix(0,0), ref_matrix(0,0));
  BOOST_CHECK_EQUAL(matrix(0,1), ref_matrix(0,1));
  BOOST_CHECK_EQUAL(matrix(1,0), ref_matrix(1,0));
  BOOST_CHECK_EQUAL(matrix(1,1), ref_matrix(1,1));

  // calculate with wrong dimension of x
  // should return 2x2 matrix filled with zeros
  cp_diff.calculate(x, matrix);
  BOOST_CHECK_EQUAL(matrix.shape(0), 2);
  BOOST_CHECK_EQUAL(matrix.shape(1), 2);
  BOOST_CHECK_EQUAL(matrix(0,0), 0.);
  BOOST_CHECK_EQUAL(matrix(0,1), 0.);
  BOOST_CHECK_EQUAL(matrix(1,0), 0.);
  BOOST_CHECK_EQUAL(matrix(1,1), 0.);
}

BOOST_AUTO_TEST_CASE( SquaredDiffCalculator_test) {
  LOG(logINFO) << "test case: SquaredDiffCalculator_test";

  // get the test data
  FeatureMatrix x;
  get_feature_matrix(x);

  // set up calculator
  SquaredDiffCalculator squared_diff_calculator;
  FeatureMatrix sq_diff;
  squared_diff_calculator.calculate(x, sq_diff);

  BOOST_CHECK_EQUAL(sq_diff.shape(0), 6);
  BOOST_CHECK_EQUAL(sq_diff.shape(1), 1);
  BOOST_CHECK_EQUAL(sq_diff(0, 0), 1.);
  BOOST_CHECK_EQUAL(sq_diff(1, 0), 2.);
  BOOST_CHECK_EQUAL(sq_diff(2, 0), 1.);
  BOOST_CHECK_EQUAL(sq_diff(3, 0), 1.);
  BOOST_CHECK_EQUAL(sq_diff(4, 0), 1.);
  BOOST_CHECK_EQUAL(sq_diff(5, 0),25.);
}

BOOST_AUTO_TEST_CASE( DiffusionCalculator_test ) {
  LOG(logINFO) << "test case: DiffusionCalculator_test";

  // get the test data
  FeatureMatrix x;
  get_feature_matrix(x);

  // set up calculator
  DiffusionCalculator diffusion_calculator;
  FeatureMatrix diffusion;
  diffusion_calculator.calculate(x, diffusion);

  BOOST_CHECK_EQUAL(diffusion.shape(0), 1);
  BOOST_CHECK_EQUAL(diffusion.shape(1), 1);
  BOOST_CHECK_EQUAL(diffusion(0, 0), static_cast<FeatureScalar>(31./6.));
}

BOOST_AUTO_TEST_CASE( MVNOutlierCalculator_calculate_inverse_covariance_matrix ) {
  LOG(logINFO) << "test case: MVNOutlierCalculator_calculate_inverse_covariance_matrix";

  // get the test data with outlier in x6
  FeatureMatrix x;
  get_feature_matrix(x);

  // Create outlier detection
  LOG(logINFO) << "  set up the outlier calculator";
  MVNOutlierCalculator mvnoutlier;
  LOG(logINFO) << "  calculate inverse covariance matrix:";
  FeatureMatrix inv_cov;
  mvnoutlier.calculate_inverse_covariance_matrix(x, inv_cov);
  
  LOG(logINFO) << "  " << inv_cov(0, 0) << "\t" << inv_cov(1, 0);
  LOG(logINFO) << "  " << inv_cov(0, 1) << "\t" << inv_cov(1, 1);
  BOOST_CHECK(( 1.89 < inv_cov(0,0)) and (inv_cov(0,0) <  1.90));
  BOOST_CHECK((-2.13 < inv_cov(1,0)) and (inv_cov(1,0) < -2.12));
  BOOST_CHECK((-2.13 < inv_cov(0,1)) and (inv_cov(0,1) < -2.12));
  BOOST_CHECK(( 2.71 < inv_cov(1,1)) and (inv_cov(1,1) <  2.72));
}

BOOST_AUTO_TEST_CASE( MVNOutlierCalculator_calculate_outlier_badness ) {
  LOG(logINFO) << "test case: MVNOutlierCalculator_calculate_outlier_badness";

  // get the test data with outlier in x6
  FeatureMatrix x;
  get_feature_matrix(x);

  // Create outlier detection
  LOG(logINFO) << "  set up the outlier calculator";
  MVNOutlierCalculator mvnoutlier;
  LOG(logINFO) << "  calculate vector";
  FeatureMatrix matrix;
  mvnoutlier.calculate_outlier_badness(x, matrix);

  LOG(logINFO) << "  vector of outlier badnesses:";
  for(size_t i = 0; i < 7; i++) {
    LOG(logINFO) << "  " << i << " : " << matrix(i, 0);
  }
  BOOST_CHECK_EQUAL(matrix.shape(0), 7);
  BOOST_CHECK_EQUAL(matrix.shape(1), 1);
  BOOST_CHECK(( 0.7153 < matrix(0, 0)) and (matrix(0, 0) <  0.7154));
  BOOST_CHECK(( 2.1810 < matrix(1, 0)) and (matrix(1, 0) <  2.1811));
  BOOST_CHECK(( 2.9443 < matrix(2, 0)) and (matrix(2, 0) <  2.9444));
  BOOST_CHECK(( 0.1657 < matrix(3, 0)) and (matrix(3, 0) <  0.1658));
  BOOST_CHECK(( 1.1733 < matrix(4, 0)) and (matrix(4, 0) <  1.1734));
  BOOST_CHECK(( 0.3489 < matrix(5, 0)) and (matrix(5, 0) <  0.3490));
  BOOST_CHECK(( 4.4711 < matrix(6, 0)) and (matrix(6, 0) <  4.4712));
}

BOOST_AUTO_TEST_CASE( MVNOutlierCalculator_calculate ) {
  LOG(logINFO) << "test case: MVNOutlierCalculator_calculate";

  // get the test data with outlier in x6
  FeatureMatrix x;
  get_feature_matrix(x);

  // Create outlier detection
  LOG(logINFO) << "  set up the outlier calculator";
  MVNOutlierCalculator mvnoutlier;
  LOG(logINFO) << "  calculate count of outliers normalized to the length";
  FeatureMatrix s3;
  mvnoutlier.calculate(x, s3);

  FeatureMatrix s4;
  mvnoutlier.calculate(x, s4, 4.0);

  BOOST_CHECK_EQUAL(s3.shape(0), 1);
  BOOST_CHECK_EQUAL(s3.shape(1), 1);
  BOOST_CHECK_EQUAL(s4.shape(0), 1);
  BOOST_CHECK_EQUAL(s4.shape(1), 1);
  BOOST_CHECK_EQUAL(s3(0, 0), static_cast<FeatureScalar>(1./7.));
  BOOST_CHECK_EQUAL(s4(0, 0), static_cast<FeatureScalar>(1./7.));
}

// BOOST_AUTO_TEST_CASE( GraphFeatureCalculator_calculate_vector ) {
//   LOG(logINFO) << "test case: GraphFeatureCalculator_calculate_vector";
//   
//   // set up the graph and set the solution index
//   HypothesesGraph graph;
//   get_graph(graph);
//   set_solution(graph, 0);
// 
//   // set up the calculators
//   boost::shared_ptr<SubsetsOfInterest> track_subsets_extractor_ptr(
//     new TrackSubsets
//   );
//   boost::shared_ptr<SubsetFeaturesIdentity> features_identity_extractor_ptr(
//     new SubsetFeaturesIdentity("id")
//   );
//   boost::shared_ptr<SubsetFeatureCalculator> sum_calculator_ptr(
//     new SumCalculator
//   );
//   GraphFeatureCalculator sum_of_ids_in_track(
//     track_subsets_extractor_ptr,
//     features_identity_extractor_ptr,
//     sum_calculator_ptr
//   );
// 
//   // calculate the vector
//   FeatureVector result = sum_of_ids_in_track.calculate_vector(graph);
//   BOOST_CHECK_EQUAL(result.shape(0), 3);
//   BOOST_CHECK_EQUAL(result(0),  4.);
//   BOOST_CHECK_EQUAL(result(1), 12.);
//   BOOST_CHECK_EQUAL(result(2), 14.);
// }
