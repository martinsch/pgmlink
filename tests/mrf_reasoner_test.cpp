#define BOOST_TEST_MODULE new_reasoner_test

#include <vector>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/bind.hpp>

#include <lemon/color.h>
#include <lemon/graph_to_eps.h>

#include <lemon/core.h>
#include <lemon/concepts/digraph.h>
#include <lemon/list_graph.h>
#include <lemon/maps.h>

#include "graph.h"
#include "hypotheses.h"
#include "energy.h"
#include "reasoning/mrf_reasoner.h"
#include "traxels.h"

using namespace Tracking;
using namespace std;
using namespace boost;

BOOST_AUTO_TEST_CASE( HypothesesGraph_build_hyp ) {
    HypothesesGraph graph;
    graph.add_node(13);
}

BOOST_AUTO_TEST_CASE( HypothesesGraph_build_hyp2 ) {

  std::cout << "Constructing HypothesesGraph" << std::endl;
  std::cout <<  std::endl;

  typedef HypothesesGraph::ArcIt ArcIt2;
  typedef HypothesesGraph::Arc Arc;
  typedef HypothesesGraph::NodeIt NodeIt;
  typedef HypothesesGraph::Node Node;
  using lemon::INVALID;

  HypothesesGraph g;
  HypothesesGraph* graph=&g;
  graph->add(node_traxel());

  std::cout << "Adding nodes" << std::endl;
  std::cout <<  std::endl;

  Node n1=g.add_node(0);
  Node n2=g.add_node(0);
  Node n3=g.add_node(0);
  g.add_node(0);
  Node m1=g.add_node(1);
  Node m2=g.add_node(1);
  Node m3=g.add_node(1);
  Node m4=g.add_node(1);
  Node o1=g.add_node(2);
  Node o2=g.add_node(2);
  Node o3=g.add_node(2);
  Node o4=g.add_node(2);

  std::cout << "Nodes:";
  for (NodeIt i(g); i!=INVALID; ++i)
    std::cout << " " << g.id(i);
  std::cout << std::endl;


  std::cout << "Addings arcs" << std::endl;
  std::cout <<  std::endl;

  g.addArc(n1,m2);
  g.addArc(n2,m1);
  g.addArc(n3,m3);
  g.addArc(n3,m4);
  g.addArc(m1,o1);
  g.addArc(m2,o2);
  g.addArc(m4,o3);
  g.addArc(m4,o4);

  std::cout << "Arcs:";
  for (ArcIt2 i(g); i!=INVALID; ++i)
  std::cout << " (" << g.id(g.source(i)) << "," << g.id(g.target(i)) << ")";
  std::cout << std::endl;
  std::cout <<  std::endl;

  std::cout << "Constructing Reasoner" << std::endl;
  std::cout <<  std::endl;


   Traxels empty;
   ConstantEnergy e(25);
    
   SingleTimestepTraxelMrf mrf(bind<double>(e, _1, empty, empty),		//detection
			        bind<double>(e, _1, empty, empty),		//non_detection
			        bind<double>(e, _1, empty, empty),		//appearance
			        bind<double>(e, _1, empty, empty),		//disappearance
				bind<double>(e, _1, _2, empty, empty),		//move
				bind<double>(e, _1, _2, _3, empty, empty)	//division
			       );

  std::cout << "Formulating Factors" << std::endl;
  std::cout <<  std::endl;

    cout << "-> workflow: formulating model" << endl; 
    mrf.formulate( *graph );
    cout << "-> workflow: infer" << endl; 
    mrf.infer();
    cout << "-> workflow: conclude" << endl; 
    mrf.conclude(*graph);
    prune_inactive(*graph);
}


BOOST_AUTO_TEST_CASE( HypothesesGraph_build_hyp_mrf ) {


    Traxel n1, n2, n3, n4, m1, m2, m3, m4, o1, o2, o3, o4, o5;
    feature_array comn1(feature_array::difference_type(3)); //float vector
    feature_array comn2(feature_array::difference_type(3));
    feature_array comn3(feature_array::difference_type(3));
    feature_array comn4(feature_array::difference_type(3));
    feature_array comm1(feature_array::difference_type(3)); 
    feature_array comm2(feature_array::difference_type(3));
    feature_array comm3(feature_array::difference_type(3));
    feature_array comm4(feature_array::difference_type(3));    
    feature_array como1(feature_array::difference_type(3)); 
    feature_array como2(feature_array::difference_type(3));
    feature_array como3(feature_array::difference_type(3));
    feature_array como4(feature_array::difference_type(3));
    feature_array como5(feature_array::difference_type(3));

    comn1[0] = 0;
    comn1[1] = 0;
    comn1[2] = 0;
    n1.features["com"] = comn1;
    n1.Id = 11;
    n1.Timestep = 0;

    comn2[0] = 1;
    comn2[1] = 1;
    comn2[2] = 1;
    n2.features["com"] = comn2;
    n2.Id = 12;
    n2.Timestep = 0;

    comn3[0] = 0;
    comn3[1] = 0;
    comn3[2] = 2;
    n3.features["com"] = comn3;
    n3.Id = 13;
    n3.Timestep = 0;

    comn4[0] = 2;
    comn4[1] = 0;
    comn4[2] = 0;
    n4.features["com"] = comn4;
    n4.Id = 14;
    n4.Timestep = 0;

    comm1[0] = 2;
    comm1[1] = 2;
    comm1[2] = 2;
    m1.features["com"] = comm1;
    m1.Id = 21;
    m1.Timestep = 1;

    comm2[0] = 1;
    comm2[1] = 0;
    comm2[2] = 1;
    m2.features["com"] = comm2;
    m2.Id = 22;
    m2.Timestep = 1;

    comm3[0] = 0;
    comm3[1] = 1;
    comm3[2] = 1;
    m3.features["com"] = comm3;
    m3.Id = 23;
    m3.Timestep = 1;

    comm4[0] = 3;
    comm4[1] = 3;
    comm4[2] = 3;
    m4.features["com"] = comm4;
    m4.Id = 24;
    m4.Timestep = 1;

    como1[0] = 4;
    como1[1] = 4;
    como1[2] = 4;
    o1.features["com"] = como1;
    o1.Id = 31;
    o1.Timestep = 2;

    como2[0] = 7;
    como2[1] = 19;
    como2[2] = 8;
    o2.features["com"] = como2;
    o2.Id = 32;
    o2.Timestep = 2;

    como3[0] = 0;
    como3[1] = 7;
    como3[2] = 1;
    o3.features["com"] = como3;
    o3.Id = 33;
    o3.Timestep = 2;

    como4[0] = 5;
    como4[1] = 0;
    como4[2] = 9;
    o4.features["com"] = como4;
    o4.Id = 34;
    o4.Timestep = 2;

    como5[0] = 0;
    como5[1] = 4;
    como5[2] = 0;
    o5.features["com"] = como5;
    o5.Id = 34;
    o5.Timestep = 2;

    TraxelStore ts;
    add(ts,n1);
    add(ts,n2);
    add(ts,n3);
    add(ts,n4);
    add(ts,m1);
    add(ts,m2);
    add(ts,m3);
    add(ts,m4);
    add(ts,o1);
    add(ts,o2);
    add(ts,o3);
    add(ts,o4);
    add(ts,o5);

    SingleTimestepTraxel_HypothesesBuilder builder(&ts);
    HypothesesGraph* graph = builder.build();

    Traxels empty;
    ConstantEnergy e1(10);	
    ConstantEnergy e2(90);	
    ConstantEnergy e3(70);	
    ConstantEnergy e4(50);	
    ConstantEnergy e5(20);		
    ConstantEnergy e6(5);	
    
    SingleTimestepTraxelMrf mrf(bind<double>(e1, _1, empty, empty),		//detection
			        bind<double>(e2, _1, empty, empty),		//non_detection
			        bind<double>(e3, _1, empty, empty),		//appearance
			        bind<double>(e4, _1, empty, empty),		//disappearance
				bind<double>(e5, _1, _2, empty, empty),		//move
				bind<double>(e6, _1, _2, _3, empty, empty)	//division
			       );
    cout << "-> workflow: formulating model" << endl; 
    mrf.formulate( *graph );
    cout << "-> workflow: infer" << endl; 
    mrf.infer();
    cout << "-> workflow: conclude" << endl; 
    mrf.conclude(*graph);
    prune_inactive(*graph);
}
