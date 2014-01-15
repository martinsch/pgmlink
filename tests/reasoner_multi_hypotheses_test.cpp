#define BOOST_TEST_MODULE reasoner_multi_hypotheses_test

// stl
#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <stdexcept>


// boost
#include <boost/test/unit_test.hpp>
#include <boost/shared_ptr.hpp>

// pgmlink
#include "pgmlink/graph.h"
#include "pgmlink/multi_hypotheses_graph.h"
#include "pgmlink/feature.h"
#include "pgmlink/reasoner_multi_hypotheses.h"
#include "pgmlink/traxels.h"

using namespace pgmlink;
typedef MultiHypothesesGraph::Node Node;
typedef MultiHypothesesGraph::Arc Arc;



BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_build_hyp4 ) {
  std::cout << "MultiHypothesesGraph_build_hyp4" << std::endl;
  // like 4, but use conflict sets
  std::cout << "Constructing MultiHypothesesGraph" << std::endl;
  std::cout << std::endl;

  MultiHypothesesGraph g;
  g.add(node_conflict_sets());
  MultiHypothesesGraph::ContainedRegionsMap& regions = g.get(node_regions_in_component());
  MultiHypothesesGraph::ConflictSetMap& conflict_sets = g.get(node_conflict_sets());


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

  std::vector<std::vector<unsigned> >& c1 = conflict_sets.get_value(n1);
  std::vector<std::vector<unsigned> >& c2 = conflict_sets.get_value(n2);
  std::vector<std::vector<unsigned> >& c3 = conflict_sets.get_value(n3);

  feature_array com(3,1.);
  feature_array count(1, 0.1);





  t1.push_back(Traxel(1, 0));
  (t1.end()-1)->features["conflicts"] = feature_array();
  (t1.end()-1)->features["level"].push_back(0.);
  (t1.end()-1)->features["com"] = com;
  (t1.end()-1)->features["count_prediction"] = count;

  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);


  t2.push_back(Traxel(1, 1));
  float conflict_arr21[] = {3., 4.};
  std::fill(com.begin(), com.end(), 0.); com[0] = 1;
  (t2.end()-1)->features["conflicts"] = feature_array(conflict_arr21, conflict_arr21 + 2);
  (t2.end()-1)->features["level"].push_back(0);
  (t2.end()-1)->features["com"] = com;
  (t2.end()-1)->features["count_prediction"] = count;

  t2.push_back(Traxel(3, 1));
  float conflict_arr23[] = {1.};
  com[0] = 0;
  (t2.end()-1)->features["conflicts"] = feature_array(conflict_arr23, conflict_arr23 + 1);
  (t2.end()-1)->features["level"].push_back(1);
  (t2.end()-1)->features["com"] = com;

  t2.push_back(Traxel(4, 1));
  com[0] = 2;
  float conflict_arr24[] = {1.};
  (t2.end()-1)->features["conflicts"] = feature_array(conflict_arr24, conflict_arr24 + 1);
  (t2.end()-1)->features["level"].push_back(1);
  (t2.end()-1)->features["com"] = com;

  c2.push_back(std::vector<unsigned>());
  c2.rbegin()->push_back(1);
  c2.rbegin()->push_back(3);
  c2.push_back(std::vector<unsigned>());
  c2.rbegin()->push_back(1);
  c2.rbegin()->push_back(4);



  t3.push_back(Traxel(2, 1));
  std::fill(com.begin(), com.end(), 2.);
  (t3.end()-1)->features["conflicts"] = feature_array();
  (t3.end()-1)->features["level"].push_back(0);
  (t3.end()-1)->features["com"] = com;
  (t3.end()-1)->features["count_prediction"] = count;

  c3.push_back(std::vector<unsigned>());
  c3.rbegin()->push_back(2);


  // add cardinalities
  g.add_cardinalities(n1);
  g.add_cardinalities(n2);
  g.add_cardinalities(n3);

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
      .with_maximal_conflict_cliques(true);

  MultiHypotheses reasoner(builder,
                           true, // with_constraints
                           0. // ep_gap
                           );

  std::cout << " -> workflow: formulating model" << std::endl;
  reasoner.formulate( g );

  ////
  //// Topology
  ////
  const pgm::OpengmModel* model = reasoner.get_graphical_model();
  std::cout << "Checking the topology of the graphical model...\n";
  std::cout << "Number of variables: " << model->numberOfVariables() << '\n';
  std::cout << "Number of factors:   " << model->numberOfFactors() << '\n';

  BOOST_CHECK_EQUAL( model->numberOfVariables(), 9);


  std::cout << " -> workflow: infer" << std::endl;
//  double objective = reasoner.infer();
  reasoner.infer();

//  BOOST_CHECK_CLOSE(objective, 3021.4142, 0.0001);

  std::cout << " -> workflow: conclude" << std::endl;
  reasoner.conclude( g );

  for (MultiHypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
    std::vector<Traxel>& traxels = regions.get_value(n);
    std::cout << "Region " << traxels[0].Id << " at time " << traxels[0].Timestep << '\n';
    for (std::vector<Traxel>::iterator t = traxels.begin(); t != traxels.end(); ++t) {
      std::cout << *t << " is active? " << t->features["active"][0];
      if (t->features["active"][0] > 0.) {
        std::cout << "   descendants: ";
        std::ostream_iterator<feature_type> os_it(std::cout, ", ");
        std::copy(t->features["outgoing"].begin(),
                  t->features["outgoing"].end(),
                  os_it);
        std::cout << " parent: ";
        std::copy(t->features["parent"].begin(),
                  t->features["parent"].end(),
                  os_it);
      }
      std::cout << '\n';
    }
    std::cout << '\n';
  }


}


BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_build_hyp5 ) {
  std::cout << "MultiHypothesesGraph_build_hyp5" << std::endl;
  // like the first test case, but with maximal conflict cliques
  std::cout << "Constructing MultiHypothesesGraph" << std::endl;
  std::cout << std::endl;

  MultiHypothesesGraph g;
  g.add(node_conflict_sets());
  MultiHypothesesGraph::ContainedRegionsMap& regions = g.get(node_regions_in_component());
  MultiHypothesesGraph::ConflictSetMap& conflict_sets = g.get(node_conflict_sets());

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

  std::vector<std::vector<unsigned> >& c1 = conflict_sets.get_value(n1);
  std::vector<std::vector<unsigned> >& c2 = conflict_sets.get_value(n2);
  std::vector<std::vector<unsigned> >& c3 = conflict_sets.get_value(n3);

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

  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);
  c1.rbegin()->push_back(2);
  c1.rbegin()->push_back(4);
  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);
  c1.rbegin()->push_back(2);
  c1.rbegin()->push_back(5);
  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);
  c1.rbegin()->push_back(3);


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

  c2.push_back(std::vector<unsigned>());
  c2.rbegin()->push_back(1);
  c2.rbegin()->push_back(3);
  c2.push_back(std::vector<unsigned>());
  c2.rbegin()->push_back(1);
  c2.rbegin()->push_back(4);


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

  c3.push_back(std::vector<unsigned>());
  c3.rbegin()->push_back(2);
  c3.rbegin()->push_back(5);
  c3.push_back(std::vector<unsigned>());
  c3.rbegin()->push_back(2);
  c3.rbegin()->push_back(6);

  g.add_cardinalities(n1);
  g.add_cardinalities(n2);
  g.add_cardinalities(n3);

  ConstantFeature det(10);
  ConstantFeature mis(1000);
  ConstantFeature div(5);

  pgm::multihypotheses::CVPR2014ModelBuilder builder( ConstantFeature(10), // appearance
                                                      ConstantFeature(10), // disappearance
                                                      SquaredDistance(), // move
                                                      ConstantFeature(0), // count
                                                      0, // forbidden_cost
                                                      0., // opportunity cost
                                                      50, // max_division_level
                                                      3 // max_count
                                                      );
  builder
      .with_detection_vars(det, mis)
      .with_divisions(div)
      .with_maximal_conflict_cliques(true)
      ;

  MultiHypotheses reasoner(builder,
                           true, // with_constraints
                           0. // ep_gap
                           );

  std::cout << " -> workflow: formulating model" << std::endl;
  reasoner.formulate( g );

  ////
  //// Topology
  ////
  const pgm::OpengmModel* model = reasoner.get_graphical_model();
  std::cout << "Checking the topology of the graphical model...\n";
  std::cout << "Number of variables: " << model->numberOfVariables() << '\n';
  std::cout << "Number of factors:   " << model->numberOfFactors() << '\n';

  BOOST_CHECK_EQUAL( model->numberOfVariables(), 41);


  std::cout << " -> workflow: infer" << std::endl;
  double objective = reasoner.infer();

  BOOST_CHECK_CLOSE(objective, 4075., 0.0001);

  std::cout << " -> workflow: conclude" << std::endl;
  reasoner.conclude( g );

  for (MultiHypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
    std::vector<Traxel>& traxels = regions.get_value(n);
    std::cout << "Region " << traxels[0].Id << " at time " << traxels[0].Timestep << '\n';
    for (std::vector<Traxel>::iterator t = traxels.begin(); t != traxels.end(); ++t) {
      std::cout << *t << " is active? " << t->features["active"][0];
      if (t->features["active"][0] > 0.) {
        std::cout << "   descendants: ";
        std::ostream_iterator<feature_type> os_it(std::cout, ", ");
        std::copy(t->features["outgoing"].begin(),
                  t->features["outgoing"].end(),
                  os_it);
        std::cout << " parent: ";
        std::copy(t->features["parent"].begin(),
                  t->features["parent"].end(),
                  os_it);
      }
      std::cout << '\n';
    }
    std::cout << '\n';
  }
}





BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_build_hyp6 ) {
  std::cout << "MultiHypothesesGraph_build_hyp6" << std::endl;
  // like the first test case, but with maximal conflict cliques

  MultiHypothesesTraxelStore ts;

  std::vector<Traxel>& t1 = ts.map[0][1];
  std::vector<Traxel>& t2 = ts.map[1][1];
  std::vector<Traxel>& t3 = ts.map[1][2];

  std::vector<std::vector<unsigned> >& c1 = ts.conflicts_by_timestep[0][1];
  std::vector<std::vector<unsigned> >& c2 = ts.conflicts_by_timestep[1][1];
  std::vector<std::vector<unsigned> >& c3 = ts.conflicts_by_timestep[1][2];

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

  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);
  c1.rbegin()->push_back(2);
  c1.rbegin()->push_back(4);
  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);
  c1.rbegin()->push_back(2);
  c1.rbegin()->push_back(5);
  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);
  c1.rbegin()->push_back(3);


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

  c2.push_back(std::vector<unsigned>());
  c2.rbegin()->push_back(1);
  c2.rbegin()->push_back(3);
  c2.push_back(std::vector<unsigned>());
  c2.rbegin()->push_back(1);
  c2.rbegin()->push_back(4);


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

  c3.push_back(std::vector<unsigned>());
  c3.rbegin()->push_back(2);
  c3.rbegin()->push_back(5);
  c3.push_back(std::vector<unsigned>());
  c3.rbegin()->push_back(2);
  c3.rbegin()->push_back(6);


  std::cout << " -> workflow: building MultiHypothesesGraph" <<std::endl;
  MultiHypothesesGraphBuilder::Options options(2, 1000000);
  MultiHypothesesGraphBuilder graph_builder(options);
  MultiHypothesesGraphPtr graph = graph_builder.build(ts);
  MultiHypothesesGraph& g = *graph;




  ConstantFeature det(10);
  ConstantFeature mis(1000);
  ConstantFeature div(5);

  pgm::multihypotheses::CVPR2014ModelBuilder builder( ConstantFeature(10), // appearance
                                                      ConstantFeature(10), // disappearance
                                                      SquaredDistance(), // move
                                                      ConstantFeature(0), // count
                                                      0, // forbidden_cost
                                                      0., // opportunity cost
                                                      50, // max_division_level
                                                      3 // max_count
                                                      );
  builder
      .with_detection_vars(det, mis)
      .with_divisions(div)
      .with_maximal_conflict_cliques(true)
      ;

  MultiHypotheses reasoner(builder,
                           true, // with_constraints
                           0. // ep_gap
                           );

  std::cout << " -> workflow: formulating model" << std::endl;
  reasoner.formulate( g );

  ////
  //// Topology
  ////
  const pgm::OpengmModel* model = reasoner.get_graphical_model();
  std::cout << "Checking the topology of the graphical model...\n";
  std::cout << "Number of variables: " << model->numberOfVariables() << '\n';
  std::cout << "Number of factors:   " << model->numberOfFactors() << '\n';

  BOOST_CHECK_EQUAL( model->numberOfVariables(), 41);


  std::cout << " -> workflow: infer" << std::endl;
//  double objective = reasoner.infer();
  reasoner.infer();

//  BOOST_CHECK_CLOSE(objective, 4075., 0.0001);

  std::cout << " -> workflow: conclude" << std::endl;
  reasoner.conclude( g );

  MultiHypothesesGraph::ContainedRegionsMap& regions = g.get(node_regions_in_component());

  for (MultiHypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
    std::vector<Traxel>& traxels = regions.get_value(n);
    std::cout << "Region " << traxels[0].Id << " at time " << traxels[0].Timestep << '\n';
    for (std::vector<Traxel>::iterator t = traxels.begin(); t != traxels.end(); ++t) {
      std::cout << *t << " is active? " << t->features["active"][0];
      if (t->features["active"][0] > 0.) {
        std::cout << "   descendants: ";
        std::ostream_iterator<feature_type> os_it(std::cout, ", ");
        std::copy(t->features["outgoing"].begin(),
                  t->features["outgoing"].end(),
                  os_it);
        std::cout << " parent: ";
        std::copy(t->features["parent"].begin(),
                  t->features["parent"].end(),
                  os_it);
      }
      std::cout << '\n';
    }
    std::cout << '\n';
  }
}



BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_build_hierarchical_count_factor ) {
  std::cout << "MultiHypothesesGraph_hierarchical_count_factor" << std::endl;
  // like the first test case, but with maximal conflict cliques

  MultiHypothesesTraxelStore ts;

  std::vector<Traxel>& t1 = ts.map[0][1];
  std::vector<Traxel>& t2 = ts.map[1][1];
  std::vector<Traxel>& t3 = ts.map[1][2];

  std::vector<std::vector<unsigned> >& c1 = ts.conflicts_by_timestep[0][1];
  std::vector<std::vector<unsigned> >& c2 = ts.conflicts_by_timestep[1][1];
  std::vector<std::vector<unsigned> >& c3 = ts.conflicts_by_timestep[1][2];

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

  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);
  c1.rbegin()->push_back(2);
  c1.rbegin()->push_back(4);
  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);
  c1.rbegin()->push_back(2);
  c1.rbegin()->push_back(5);
  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);
  c1.rbegin()->push_back(3);


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

  c2.push_back(std::vector<unsigned>());
  c2.rbegin()->push_back(1);
  c2.rbegin()->push_back(3);
  c2.push_back(std::vector<unsigned>());
  c2.rbegin()->push_back(1);
  c2.rbegin()->push_back(4);


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

  c3.push_back(std::vector<unsigned>());
  c3.rbegin()->push_back(2);
  c3.rbegin()->push_back(5);
  c3.push_back(std::vector<unsigned>());
  c3.rbegin()->push_back(2);
  c3.rbegin()->push_back(6);


  std::cout << " -> workflow: building MultiHypothesesGraph" <<std::endl;
  MultiHypothesesGraphBuilder::Options options(2, 1000000);
  MultiHypothesesGraphBuilder graph_builder(options);
  MultiHypothesesGraphPtr graph = graph_builder.build(ts);
  MultiHypothesesGraph& g = *graph;




  ConstantFeature det(10);
  ConstantFeature mis(1000);
  ConstantFeature div(5);

  pgm::multihypotheses::CVPR2014ModelBuilder builder( ConstantFeature(10), // appearance
                                                      ConstantFeature(10), // disappearance
                                                      SquaredDistance(), // move
                                                      ConstantFeature(0), // count
                                                      0, // forbidden_cost
                                                      0., // opportunity cost
                                                      50, // max_division_level
                                                      3 // max_count
                                                      );
  builder
      .with_detection_vars(det, mis)
      .with_divisions(div)
      .with_maximal_conflict_cliques(true)
      .with_hierarchical_counting_factor(true)
      ;

  MultiHypotheses reasoner(builder,
                           true, // with_constraints
                           0. // ep_gap
                           );

  std::cout << " -> workflow: formulating model" << std::endl;
  reasoner.formulate( g );

  ////
  //// Topology
  ////
  const pgm::OpengmModel* model = reasoner.get_graphical_model();
  std::cout << "Checking the topology of the graphical model...\n";
  std::cout << "Number of variables: " << model->numberOfVariables() << '\n';
  std::cout << "Number of factors:   " << model->numberOfFactors() << '\n';



  std::cout << " -> workflow: infer" << std::endl;
//  double objective = reasoner.infer();
  reasoner.infer();

//  BOOST_CHECK_CLOSE(objective, 4075., 0.0001);
  // need to figure out objective deviation!

  std::cout << " -> workflow: conclude" << std::endl;
  reasoner.conclude( g );

  MultiHypothesesGraph::ContainedRegionsMap& regions = g.get(node_regions_in_component());

  for (MultiHypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
    std::vector<Traxel>& traxels = regions.get_value(n);
    std::cout << "Region " << traxels[0].Id << " at time " << traxels[0].Timestep << '\n';
    for (std::vector<Traxel>::iterator t = traxels.begin(); t != traxels.end(); ++t) {
      std::cout << *t << " is active? " << t->features["active"][0];
      if (t->features["active"][0] > 0.) {
        std::cout << "   descendants: ";
        std::ostream_iterator<feature_type> os_it(std::cout, ", ");
        std::copy(t->features["outgoing"].begin(),
                  t->features["outgoing"].end(),
                  os_it);
        std::cout << " parent: ";
        std::copy(t->features["parent"].begin(),
                  t->features["parent"].end(),
                  os_it);
      }
      std::cout << '\n';
    }
    std::cout << '\n';
  }
}



BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_build_division_segments ) {
  std::cout << "MultiHypothesesGraph_build_division_segments" << std::endl;

  MultiHypothesesTraxelStore ts;

  //	   --o
  //      |
  // o|o <
  //      |
  //       --o


  std::vector<Traxel>& t1 = ts.map[0][3]; // [timestep][connected_comp]
  std::vector<Traxel>& t2 = ts.map[1][4]; // [timestep][connected_comp]
  std::vector<Traxel>& t3 = ts.map[1][5]; // [timestep][connected_comp]

  std::vector<std::vector<unsigned> >& c1 = ts.conflicts_by_timestep[0][3];
  std::vector<std::vector<unsigned> >& c2 = ts.conflicts_by_timestep[1][4];
  std::vector<std::vector<unsigned> >& c3 = ts.conflicts_by_timestep[1][5];

  feature_array com(3,1.);
  feature_array count;
  feature_array detProb(2, 0.);

  // connected component 3, t=0
  t1.push_back(Traxel(3, 0));
  float conflict_arr13[] = {1., 2.};
  (t1.end()-1)->features["conflicts"] = feature_array(conflict_arr13, conflict_arr13 + 2);
  (t1.end()-1)->features["level"].push_back(0.);  // connected component
  (t1.end()-1)->features["com"] = com;
  detProb[0] = 0.03;
  detProb[1] = 0.97;
  (t1.end()-1)->features["detProb"] = detProb;
  count.clear();
  count.push_back(0.02); // 0
  count.push_back(0.8);  // 1
  count.push_back(0.18); // 2
  (t1.end()-1)->features["count_prediction"] = count;

  // segment 2 of connected component 3, t=0
  t1.push_back(Traxel(2, 0));
  float conflict_arr12[] = {3.};
  std::fill(com.begin(), com.end(), 2.);
  (t1.end()-1)->features["conflicts"] = feature_array(conflict_arr12, conflict_arr12 + 1);
  (t1.end()-1)->features["level"].push_back(1.);
  (t1.end()-1)->features["com"] = com;
  detProb[0] = 0.14;
  detProb[1] = 0.86;
  (t1.end()-1)->features["detProb"] = detProb;

  // segment 1 of connected component 3, t=0
  t1.push_back(Traxel(1, 0));
  float conflict_arr11[] = {1.};
  std::fill(com.begin(), com.end(), 0.);
  (t1.end()-1)->features["conflicts"] = feature_array(conflict_arr11, conflict_arr11 + 1);
  (t1.end()-1)->features["level"].push_back(1.);
  (t1.end()-1)->features["com"] = com;
  detProb[0] = 0.11;
  detProb[1] = 0.89;
  (t1.end()-1)->features["detProb"] = detProb;

  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);
  c1.rbegin()->push_back(3);
  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(2);
  c1.rbegin()->push_back(3);


  // connected component 4, t=1 (first daughter)
  t2.push_back(Traxel(4, 1));
  std::fill(com.begin(), com.end(), 0.); com[0] = 1;
  (t2.end()-1)->features["conflicts"] = feature_array();
  (t2.end()-1)->features["level"].push_back(0.);
  (t2.end()-1)->features["com"] = com;
  count.clear();
  count.push_back(0.02); // 0
  count.push_back(0.98);  // 1
  (t2.end()-1)->features["count_prediction"] = count;
  detProb[0] = 0.02;
  detProb[1] = 0.98;
  (t2.end()-1)->features["detProb"] = detProb;

  c2.push_back(std::vector<unsigned>());
  c2.rbegin()->push_back(4);

  // connected component 5, t=1 (second daughter)
  t3.push_back(Traxel(5, 1));
  std::fill(com.begin(), com.end(), 2.);
  (t3.end()-1)->features["conflicts"] = feature_array();
  (t3.end()-1)->features["level"].push_back(0);
  (t3.end()-1)->features["com"] = com;
  count.clear();
  count.push_back(0.02); // 0
  count.push_back(0.98);  // 1
  (t3.end()-1)->features["count_prediction"] = count;
  detProb[0] = 0.02;
  detProb[1] = 0.98;
  (t3.end()-1)->features["detProb"] = detProb;

  c3.push_back(std::vector<unsigned>());
  c3.rbegin()->push_back(5);

  std::cout << " -> workflow: building MultiHypothesesGraph" <<std::endl;
  MultiHypothesesGraphBuilder::Options options(1, 1000000, true); // 1 nearest neighbor, forward_backward
  MultiHypothesesGraphBuilder graph_builder(options);
  MultiHypothesesGraphPtr graph = graph_builder.build(ts);
  MultiHypothesesGraph& g = *graph;

//  MultiHypothesesGraph::ContainedRegionsMap& regions = g.get(node_regions_in_component());

  // add move probabilities
  g.add(node_move_features());
  MultiHypothesesGraph::MoveFeatureMap& moves = g.get(node_move_features());
  const MultiHypothesesGraph::TraxelMap& traxels = g.get(node_traxel());

  // moves from 3 (connected comp, t=0)
  MultiHypothesesGraph::Node n;
  for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
  			it != traxels.endValue(); ++it) {
  		if (*it == t1[0]) {
  			n = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
  			break;
  		}
  	}
  	std::map<Traxel, feature_array>& move_probs1 = moves.get_value(n)[t1[0]];
  	move_probs1[t2[0]].push_back(0.4);
  	move_probs1[t2[0]].push_back(0.6);
  	move_probs1[t3[0]].push_back(0.4);
  	move_probs1[t3[0]].push_back(0.6);

  	// moves from 2 (segment in 3)
  	  for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
  	  			it != traxels.endValue(); ++it) {
  	  		if (*it == t1[1]) {
  	  			n = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
  	  			break;
  	  		}
  	  	}
  	std::map<Traxel, feature_array>& move_probs2 = moves.get_value(n)[t1[1]];
	move_probs2[t2[0]].push_back(0.4);
	move_probs2[t2[0]].push_back(0.6);
	move_probs2[t3[0]].push_back(0.4);
	move_probs2[t3[0]].push_back(0.6);


  	// moves from 1 (segment in 3)
  	  for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
  	  			it != traxels.endValue(); ++it) {
  	  		if (*it == t1[2]) {
  	  			n = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
  	  			break;
  	  		}
  	  	}
  	std::map<Traxel, feature_array>& move_probs3 = moves.get_value(n)[t1[2]];
	move_probs3[t2[0]].push_back(0.4);
	move_probs3[t2[0]].push_back(0.6);
	move_probs3[t3[0]].push_back(0.4);
	move_probs3[t3[0]].push_back(0.6);


	// add division probabilities
  	g.add(node_division_features());
  	MultiHypothesesGraph::DivisionFeatureMap& divs = g.get(node_division_features());

  	// divs for connected component 3
    for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
  			it != traxels.endValue(); ++it) {
  		if (*it == t1[0]) {
  			n = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
  			break;
  		}
  	}
  	std::map<std::pair<Traxel, Traxel>, feature_array>& div_probs1 = divs.get_value(n)[t1[0]];
  	std::vector<float> probs;
  	probs.push_back(0.03);
  	probs.push_back(0.97);
  	div_probs1[std::make_pair(t3[0], t2[0])] = feature_array(probs.size());
  	std::copy(probs.begin(),
  			  probs.end(),
  			  div_probs1[std::make_pair(t3[0], t2[0])].begin());
  	div_probs1[std::make_pair(t2[0], t3[0])] = feature_array(probs.size());
	std::copy(probs.begin(),
  	  			  probs.end(),
  	  			  div_probs1[std::make_pair(t2[0], t3[0])].begin());

  	// divs for segment 2 (in connected comp. 3, t=0)
	for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
			it != traxels.endValue(); ++it) {
		if (*it == t1[1]) {
			n = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
			break;
		}
	}
	std::map<std::pair<Traxel, Traxel>, feature_array>& div_probs2 = divs.get_value(n)[t1[1]];
	probs.clear();
	probs.push_back(0.14);
	probs.push_back(0.86);
	div_probs2[std::make_pair(t3[0], t2[0])] = feature_array(probs.size());
	std::copy(probs.begin(),
			  probs.end(),
			  div_probs2[std::make_pair(t3[0], t2[0])].begin());
	div_probs2[std::make_pair(t2[0], t3[0])] = feature_array(probs.size());
	std::copy(probs.begin(),
				  probs.end(),
				  div_probs2[std::make_pair(t2[0], t3[0])].begin());

  	// divs for segment 1 (in connected comp. 3, t=0)
	for (MultiHypothesesGraph::TraxelMap::ValueIt it = traxels.beginValue();
			it != traxels.endValue(); ++it) {
		if (*it == t1[2]) {
			n = MultiHypothesesGraph::TraxelMap::ItemIt(traxels, *it);
			break;
		}
	}
	std::map<std::pair<Traxel, Traxel>, feature_array>&  div_probs3 = divs.get_value(n)[t1[2]];
	probs.clear();
	probs.push_back(0.11);
	probs.push_back(0.89);
	div_probs3[std::make_pair(t3[0], t2[0])] = feature_array(probs.size());
	std::copy(probs.begin(),
			  probs.end(),
			  div_probs3[std::make_pair(t3[0], t2[0])].begin());
	div_probs3[std::make_pair(t2[0], t3[0])] = feature_array(probs.size());
	std::copy(probs.begin(),
				  probs.end(),
				  div_probs3[std::make_pair(t2[0], t3[0])].begin());


  double maxDiv = 10;
  ConstantFeature tmp_mov(0);
  ConstantFeature tmp_cnt(0);
  ConstantFeature app(0);
  ConstantFeature dis(0);
  double opp = 500;
  double forbidden = 500;
  double outarcs = 0;

  pgm::multihypotheses::CVPR2014ModelBuilder builder( app, // appearance
                                                      dis, // disappearance
                                                      tmp_mov, // move
                                                      tmp_cnt, // count
                                                      forbidden, // forbidden_cost
                                                      opp, // opportunity cost
                                                      maxDiv, // max_division_level
                                                      outarcs // max_count
                                                      );

  function<double (const Traxel&, size_t)> det = NegLnCardinalityDetection(50);
  function<double (const Traxel&, const Traxel&, const Traxel&, feature_type)> div =
  		  NegLnDivision(1, maxDiv);
  function<double (const Traxel&, const Traxel&, feature_type)> mov = NegLnTransition(11);
  function<double (feature_type)> cnt = NegLn(50);

  builder
      .with_detection_vars(det, det)
      .with_divisions(div)
      .with_maximal_conflict_cliques(true)
      .with_hierarchical_counting_factor(true)
      .with_classifier_priors(mov, cnt)
	  .with_conflict_factors(det);
      ;

  MultiHypotheses reasoner(builder,
                           true, // with_constraints
                           0. // ep_gap
                           );

  std::cout << " -> workflow: formulating model" << std::endl;
  reasoner.formulate( g );

  ////
  //// Topology
  ////
  const pgm::OpengmModel* model = reasoner.get_graphical_model();
  std::cout << "Checking the topology of the graphical model...\n";
  std::cout << "Number of variables: " << model->numberOfVariables() << '\n';
  std::cout << "Number of factors:   " << model->numberOfFactors() << '\n';



  std::cout << " -> workflow: infer" << std::endl;
  reasoner.infer();

  std::cout << " -> workflow: conclude" << std::endl;
  reasoner.conclude( g );



  MultiHypothesesGraph::ContainedRegionsMap& regions = g.get(node_regions_in_component());

  for (MultiHypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
    std::vector<Traxel>& traxels = regions.get_value(n);
    std::cout << "Region " << traxels[0].Id << " at time " << traxels[0].Timestep << '\n';
    for (std::vector<Traxel>::iterator t = traxels.begin(); t != traxels.end(); ++t) {
      std::cout << *t << " is active? " << t->features["active"][0];
	  switch(t->Id) {
	  case 1: // fall through
	  case 2: BOOST_CHECK_EQUAL(t->features["outgoing"].size(), 0); break;
	  case 3: BOOST_CHECK_EQUAL(t->features["outgoing"].size(), 2); break;
	  case 4: // fall through
	  case 5: BOOST_CHECK_EQUAL(int(t->features["active"][0]), 1); break;
	  default: BOOST_CHECK(false); break;
      }

      if (t->features["active"][0] > 0.) {
        std::cout << "   descendants: ";
        std::ostream_iterator<feature_type> os_it(std::cout, ", ");
        std::copy(t->features["outgoing"].begin(),
                  t->features["outgoing"].end(),
                  os_it);
        std::cout << " parent: ";
        std::copy(t->features["parent"].begin(),
                  t->features["parent"].end(),
                  os_it);
      }
      std::cout << '\n';
    }
    std::cout << '\n';
  }

}


BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_build_incoming_counting_factor ) {
  std::cout << "MultiHypothesesGraph_build_incoming_counting_factor" << std::endl;
  // like the first test case, but with incoming counting factor

  MultiHypothesesTraxelStore ts;

  std::vector<Traxel>& t1 = ts.map[0][1];
  std::vector<Traxel>& t2 = ts.map[1][1];
  std::vector<Traxel>& t3 = ts.map[1][2];

  std::vector<std::vector<unsigned> >& c1 = ts.conflicts_by_timestep[0][1];
  std::vector<std::vector<unsigned> >& c2 = ts.conflicts_by_timestep[1][1];
  std::vector<std::vector<unsigned> >& c3 = ts.conflicts_by_timestep[1][2];

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

  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);
  c1.rbegin()->push_back(2);
  c1.rbegin()->push_back(4);
  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);
  c1.rbegin()->push_back(2);
  c1.rbegin()->push_back(5);
  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);
  c1.rbegin()->push_back(3);


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

  c2.push_back(std::vector<unsigned>());
  c2.rbegin()->push_back(1);
  c2.rbegin()->push_back(3);
  c2.push_back(std::vector<unsigned>());
  c2.rbegin()->push_back(1);
  c2.rbegin()->push_back(4);


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

  c3.push_back(std::vector<unsigned>());
  c3.rbegin()->push_back(2);
  c3.rbegin()->push_back(5);
  c3.push_back(std::vector<unsigned>());
  c3.rbegin()->push_back(2);
  c3.rbegin()->push_back(6);


  std::cout << " -> workflow: building MultiHypothesesGraph" <<std::endl;
  MultiHypothesesGraphBuilder::Options options(2, 1000000);
  MultiHypothesesGraphBuilder graph_builder(options);
  MultiHypothesesGraphPtr graph = graph_builder.build(ts);
  MultiHypothesesGraph& g = *graph;




  ConstantFeature det(10);
  ConstantFeature mis(1000);
  ConstantFeature div(5);

  pgm::multihypotheses::CVPR2014ModelBuilder builder( ConstantFeature(10), // appearance
                                                      ConstantFeature(10), // disappearance
                                                      SquaredDistance(), // move
                                                      ConstantFeature(0), // count
                                                      0, // forbidden_cost
                                                      0., // opportunity cost
                                                      50, // max_division_level
                                                      3 // max_count
                                                      );
  builder
      .with_detection_vars(det, mis)
      .with_divisions(div)
      .with_maximal_conflict_cliques(true)
      .with_counting_incoming_factor(true)
      ;

  MultiHypotheses reasoner(builder,
                           true, // with_constraints
                           0. // ep_gap
                           );

  std::cout << " -> workflow: formulating model" << std::endl;
  reasoner.formulate( g );

  ////
  //// Topology
  ////
  const pgm::OpengmModel* model = reasoner.get_graphical_model();
  std::cout << "Checking the topology of the graphical model...\n";
  std::cout << "Number of variables: " << model->numberOfVariables() << '\n';
  std::cout << "Number of factors:   " << model->numberOfFactors() << '\n';



  std::cout << " -> workflow: infer" << std::endl;
//  double objective = reasoner.infer();
  reasoner.infer();

//  BOOST_CHECK_CLOSE(objective, 4075., 0.0001);

  std::cout << " -> workflow: conclude" << std::endl;
  reasoner.conclude( g );

  MultiHypothesesGraph::ContainedRegionsMap& regions = g.get(node_regions_in_component());

  for (MultiHypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
    std::vector<Traxel>& traxels = regions.get_value(n);
    std::cout << "Region " << traxels[0].Id << " at time " << traxels[0].Timestep << '\n';
    for (std::vector<Traxel>::iterator t = traxels.begin(); t != traxels.end(); ++t) {
      std::cout << *t << " is active? " << t->features["active"][0];
      if (t->features["active"][0] > 0.) {
        std::cout << "   descendants: ";
        std::ostream_iterator<feature_type> os_it(std::cout, ", ");
        std::copy(t->features["outgoing"].begin(),
                  t->features["outgoing"].end(),
                  os_it);
        std::cout << " parent: ";
        std::copy(t->features["parent"].begin(),
                  t->features["parent"].end(),
                  os_it);
      }
      std::cout << '\n';
    }
    std::cout << '\n';
  }
}


BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_restrict_in_time ) {
  std::cout << "MultiHypothesesGraph_restrict_in_time" << std::endl;
  // like the first test case, but with maximal conflict cliques

  MultiHypothesesTraxelStore ts;

  std::vector<Traxel>& t1 = ts.map[0][1];
  std::vector<Traxel>& t2 = ts.map[1][1];

  std::vector<std::vector<unsigned> >& c1 = ts.conflicts_by_timestep[0][1];
  std::vector<std::vector<unsigned> >& c2 = ts.conflicts_by_timestep[1][1];

  feature_array com(3,1.);
  feature_array count(1, 0.1);

  t1.push_back(Traxel(1, 0));
  float conflict_arr11[] = {2., 3., 4., 5.};
  (t1.end()-1)->features["conflicts"] = feature_array(conflict_arr11, conflict_arr11 + 4);
  (t1.end()-1)->features["level"].push_back(0.);
  (t1.end()-1)->features["com"] = com;
  (t1.end()-1)->features["count_prediction"] = count;

  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);


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
       it != traxels.endValue();
       ++it) {
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




  ConstantFeature det(10);
  ConstantFeature mis(1000);
  ConstantFeature div(5);

  pgm::multihypotheses::CVPR2014ModelBuilder builder( ConstantFeature(10), // appearance
                                                      ConstantFeature(10), // disappearance
                                                      SquaredDistance(), // move
                                                      ConstantFeature(0), // count
                                                      0, // forbidden_cost
                                                      0., // opportunity cost
                                                      50, // max_division_level
                                                      3 // max_count
                                                      );
  builder
      .with_detection_vars(det, mis)
      .with_divisions(div)
      .with_maximal_conflict_cliques(true)
//      .with_maximum_arcs(1)
      .with_timestep_range(0, 0)
      ;

  MultiHypotheses reasoner(builder,
                           true, // with_constraints
                           0. // ep_gap
                           );

  std::cout << " -> workflow: formulating model" << std::endl;
  reasoner.formulate( g );


  ////
  //// Topology
  ////
  const pgm::OpengmModel* model = reasoner.get_graphical_model();
  std::cout << "Checking the topology of the graphical model...\n";
  std::cout << "Number of variables: " << model->numberOfVariables() << '\n';
  std::cout << "Number of factors:   " << model->numberOfFactors() << '\n';

  // 5 variables: 4*detection + 1*assignment
  BOOST_CHECK_EQUAL( model->numberOfVariables(), 1);
  // 14 factors: 2 count factors + 4*(1 outgoing, 1 incoming, 1 unary)
  BOOST_CHECK_EQUAL( model->numberOfFactors(), 4);


}


BOOST_AUTO_TEST_CASE( MultiHypothesesGraph_restrict_in_time2 ) {
  std::cout << "MultiHypothesesGraph_restrict_in_time2" << std::endl;
  // like the first test case, but with maximal conflict cliques

  MultiHypothesesTraxelStore ts;

  std::vector<Traxel>& t1 = ts.map[0][1];
  std::vector<Traxel>& t2 = ts.map[1][1];
  std::vector<Traxel>& t3 = ts.map[2][1];
  std::vector<Traxel>& t4 = ts.map[3][1];

  std::vector<std::vector<unsigned> >& c1 = ts.conflicts_by_timestep[0][1];
  std::vector<std::vector<unsigned> >& c2 = ts.conflicts_by_timestep[1][1];
  std::vector<std::vector<unsigned> >& c3 = ts.conflicts_by_timestep[2][1];
  std::vector<std::vector<unsigned> >& c4 = ts.conflicts_by_timestep[3][1];

  feature_array com(3,1.);
  feature_array count(1, 0.1);

  t1.push_back(Traxel(1, 0));
  float conflict_arr11[] = {2., 3., 4., 5.};
  (t1.end()-1)->features["conflicts"] = feature_array(conflict_arr11, conflict_arr11 + 4);
  (t1.end()-1)->features["level"].push_back(0.);
  (t1.end()-1)->features["com"] = com;
  (t1.end()-1)->features["count_prediction"] = count;

  c1.push_back(std::vector<unsigned>());
  c1.rbegin()->push_back(1);


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

  c2.push_back(std::vector<unsigned>());
  c2.rbegin()->push_back(1);
  c2.rbegin()->push_back(3);
  c2.push_back(std::vector<unsigned>());
  c2.rbegin()->push_back(1);
  c2.rbegin()->push_back(4);


  t3.push_back(Traxel(1, 2));
  std::fill(com.begin(), com.end(), 0);
  assert((t3.end()-1)->Id == 1);
  (t3.end()-1)->features["conflicts"] = feature_array(0);
  (t3.end()-1)->features["level"].push_back(0.);
  (t3.end()-1)->features["com"] = com;
  assert((t3.end()-1)->features.count("com") > 0);
  (t3.end()-1)->features["count_prediction"] = count;

  c3.push_back(std::vector<unsigned>());
  c3.rbegin()->push_back(1);


  t4.push_back(Traxel(1, 3));
  std::fill(com.begin(), com.end(), 0);
  (t4.end()-1)->features["conflicts"] = feature_array(0);
  (t4.end()-1)->features["level"].push_back(0.);
  (t4.end()-1)->features["com"] = com;
  (t4.end()-1)->features["count_prediction"] = count;

  c4.push_back(std::vector<unsigned>());
  c4.rbegin()->push_back(1);

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
       it != traxels.endValue();
       ++it) {
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




  ConstantFeature det(10);
  ConstantFeature mis(1000);
  ConstantFeature div(5);

  pgm::multihypotheses::CVPR2014ModelBuilder builder( ConstantFeature(10), // appearance
                                                      ConstantFeature(10), // disappearance
                                                      SquaredDistance(), // move
                                                      ConstantFeature(0), // count
                                                      0, // forbidden_cost
                                                      0., // opportunity cost
                                                      50, // max_division_level
                                                      3 // max_count
                                                      );
  builder
      .with_detection_vars(det, mis)
      .with_divisions(div)
      .with_maximal_conflict_cliques(true)
      .with_timestep_range(1, 2)
      ;

  MultiHypotheses reasoner(builder,
                           true, // with_constraints
                           0. // ep_gap
                           );

  std::cout << " -> workflow: formulating model" << std::endl;
  reasoner.formulate( g );


  ////
  //// Topology
  ////
  const pgm::OpengmModel* model = reasoner.get_graphical_model();
  std::cout << "Checking the topology of the graphical model...\n";
  std::cout << "Number of variables: " << model->numberOfVariables() << '\n';
  std::cout << "Number of factors:   " << model->numberOfFactors() << '\n';

  // 5 variables: 4*detection + 1*assignment
  BOOST_CHECK_EQUAL( model->numberOfVariables(), 7);
  // 14 factors: 2 count factors + 4*(1 outgoing, 1 incoming, 1 unary)
  BOOST_CHECK_EQUAL( model->numberOfFactors(), 14);


}
