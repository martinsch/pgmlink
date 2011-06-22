#define BOOST_TEST_MODULE graphical_model_test

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <vigra/multi_array.hxx>

#include "energy.h"
#include "graphical_model.h"
#include "traxels.h"
#include "event.h"

using namespace Tracking;
using namespace std;
using namespace boost;
using namespace vigra;

/********
Missing tests:
- test not fully connected models
- test energies
**********/


void print_events(shared_ptr<std::vector< std::vector<Event> > > events) {
     for(size_t i = 0; i < events->size(); ++i){
	cout << "timestep " << i << ":\n";
	for(size_t j=0; j < (*events)[i].size(); j++) {
	    cout << "\t" << (*events)[i][j] << endl;
	}
    }
}

BOOST_AUTO_TEST_CASE( opengm_factor_table_ordering )
{
  OpengmMrf mrf;
  mrf.Space()->addDimension(2);
  mrf.Space()->addDimension(2);
  mrf.Space()->addDimension(2);
  size_t vi[] = {0,1,2};
  size_t vv[] = {1,2,3,4,5,6,7,8};
  OpengmMrf::ogmFactor f(*(mrf.Space()), vi, vi+3, vv, vv+8);
  BOOST_CHECK_EQUAL(f(0,0,0), 1);
  BOOST_CHECK_EQUAL(f(1,0,0), 2);
  BOOST_CHECK_EQUAL(f(0,1,0), 3);
  BOOST_CHECK_EQUAL(f(1,1,0), 4);
  BOOST_CHECK_EQUAL(f(0,0,1), 5);
  BOOST_CHECK_EQUAL(f(1,0,1), 6);
  BOOST_CHECK_EQUAL(f(0,1,1), 7);
  BOOST_CHECK_EQUAL(f(1,1,1), 8);
}

BOOST_AUTO_TEST_CASE( simple_move )
{
    //
    // tr11 --- tr21
    //

    // trXY: traxel Y at timestep X
    Traxel tr11, tr21;
    feature_array com11(feature_array::difference_type(3));
    feature_array com21(feature_array::difference_type(3));

    com11[0] = 0;
    com11[1] = 0;
    com11[2] = 0;
    tr11.features["com"] = com11;
    tr11.Id = 5;

    com21[0] = 1;
    com21[1] = 0;
    com21[2] = 0;
    tr21.features["com"] = com21;
    tr21.Id = 7;

    vector< vector<Traxel> > traxels;
    traxels.push_back(vector<Traxel>());
    traxels.push_back(vector<Traxel>());
    traxels[0].push_back(tr11);
    traxels[1].push_back(tr21);

    //----

    // objects for reuse
    shared_ptr<GraphicalModel> gm;
    ConstantEnergy div, mov, mismov, det, misdet;
    Event expected1, expected2;
    shared_ptr<std::vector< std::vector<Event> > > events;

    // cheap misdetections, cheap move -> expect disappearance
    mov.theEnergy = 0;
    mismov.theEnergy = 10;
    det.theEnergy = 10;
    misdet.theEnergy = 0;
    gm = GraphicalModel::from(traxels, div,mov,mismov,det,misdet);
    gm->infer();
    events = gm->events();
    BOOST_CHECK_EQUAL(events->size(), 1);
    BOOST_CHECK_EQUAL((*events)[0].size(), 2);
    expected1.type = Event::Disappearance;
    expected1.traxel_ids.push_back(tr11.Id);
    expected2.type = Event::Appearance;
    expected2.traxel_ids.push_back(tr21.Id);

    BOOST_CHECK((*events)[0][0] == expected1);
    BOOST_CHECK((*events)[0][1] == expected2);

    print_events(events);
    cout << "\n";

    // cheap detections, cheap move -> expect move
    mov.theEnergy = 1;
    mismov.theEnergy = 10;
    det.theEnergy = 1;
    misdet.theEnergy = 10;
    gm = GraphicalModel::from(traxels, div,mov,mismov,det,misdet);
    gm->infer();
    events = gm->events();
    BOOST_CHECK_EQUAL(events->size(), 1);
    BOOST_CHECK_EQUAL((*events)[0].size(), 1);
    expected1.type = Event::Move;
    expected1.traxel_ids.clear();
    expected1.traxel_ids.push_back(tr11.Id);
    expected1.traxel_ids.push_back(tr21.Id);
    BOOST_CHECK((*events)[0][0] == expected1);
    //BOOST_CHECK_EQUAL((*events)[0][0].energy, expected.energy);
    
    print_events(events);
    cout << "\n";

    // cheap detections, expensive move -> expect disappearance
    mov.theEnergy = 100;
    mismov.theEnergy = 10;
    det.theEnergy = 1;
    misdet.theEnergy = 10;
    gm = GraphicalModel::from(traxels, div,mov,mismov,det,misdet);
    gm->infer();
    events = gm->events();
    BOOST_CHECK_EQUAL(events->size(), 1);
    BOOST_CHECK_EQUAL((*events)[0].size(), 2);
    expected1.type = Event::Disappearance;
    expected1.traxel_ids.clear();
    expected1.traxel_ids.push_back(tr11.Id);
    expected2.traxel_ids.clear();
    expected2.type = Event::Appearance;
    expected2.traxel_ids.push_back(tr21.Id);

    BOOST_CHECK((*events)[0][0] == expected1);
    BOOST_CHECK((*events)[0][1] == expected2);
    
    print_events(events);
    cout << "\n";
}

BOOST_AUTO_TEST_CASE( four_pattern )
{
    //
    // tr11  tr21
    // tr12  tr22
    //

    // trXY: traxel Y at timestep X
    Traxel tr11, tr12, tr21, tr22;
    feature_array com11(feature_array::difference_type(3));
    feature_array com12(feature_array::difference_type(3));
    feature_array com21(feature_array::difference_type(3));
    feature_array com22(feature_array::difference_type(3));

    com11[0] = 0;
    com11[1] = 0;
    com11[2] = 0;
    tr11.features["com"] = com11;
    tr11.Id = 5;

    com12[0] = 0;
    com12[1] = 1;
    com12[2] = 0;
    tr12.features["com"] = com12;
    tr12.Id = 7;

    com21[0] = 1;
    com21[1] = 0;
    com21[2] = 0;
    tr21.features["com"] = com21;
    tr21.Id = 9;

    com22[0] = 1;
    com22[1] = 1;
    com22[2] = 0;
    tr22.features["com"] = com22;
    tr22.Id = 11;

    vector< vector<Traxel> > traxels;
    traxels.push_back(vector<Traxel>());
    traxels.push_back(vector<Traxel>());
    traxels[0].push_back(tr11);
    traxels[0].push_back(tr12);
    traxels[1].push_back(tr21);
    traxels[1].push_back(tr22);

    //----
    
    // objects for reuse
    shared_ptr<GraphicalModel> gm;
    ConstantEnergy div, mov, mismov, det, misdet;
    Event expected1, expected2, expected3, expected4;
    shared_ptr<std::vector< std::vector<Event> > > events;

    // cheap misdetections, cheap moves -> expect disappearances
    div.theEnergy = 100;
    mov.theEnergy = 0;
    mismov.theEnergy = 10;
    det.theEnergy = 100;
    misdet.theEnergy = 0;

    gm = GraphicalModel::from(traxels, div,mov,mismov,det,misdet);
    gm->infer();
    events = gm->events();
    BOOST_CHECK_EQUAL(events->size(), 1);
    BOOST_CHECK_EQUAL((*events)[0].size(), 4);
    expected1.type = Event::Disappearance;
    expected1.traxel_ids.clear();
    expected1.traxel_ids.push_back(tr11.Id);
    expected2.type = Event::Disappearance;
    expected2.traxel_ids.clear();
    expected2.traxel_ids.push_back(tr12.Id);
    expected3.type = Event::Appearance;
    expected3.traxel_ids.clear();
    expected3.traxel_ids.push_back(tr21.Id);
    expected4.type = Event::Appearance;
    expected4.traxel_ids.clear();
    expected4.traxel_ids.push_back(tr22.Id);

    BOOST_CHECK((*events)[0][0] == expected1);
    BOOST_CHECK((*events)[0][1] == expected2);
    BOOST_CHECK((*events)[0][2] == expected3);
    BOOST_CHECK((*events)[0][3] == expected4);
    print_events(events);
    cout << "\n";

    // cheap detections, cheap moves -> expect two moves
    div.theEnergy = 100;
    mov.theEnergy = 3;
    mismov.theEnergy = 10;
    det.theEnergy = 5;
    misdet.theEnergy = 50;

    gm = GraphicalModel::from(traxels, div,mov,mismov,det,misdet);
    gm->infer();
    events = gm->events();
    BOOST_CHECK_EQUAL(events->size(), 1);
    BOOST_CHECK_EQUAL((*events)[0].size(), 2);
    expected1.type = Event::Move;
    expected1.traxel_ids.clear();
    expected1.traxel_ids.push_back(tr11.Id);
    expected1.traxel_ids.push_back(tr22.Id);
    expected2.type = Event::Move;
    expected2.traxel_ids.clear();
    expected2.traxel_ids.push_back(tr12.Id);
    expected2.traxel_ids.push_back(tr21.Id);
    BOOST_CHECK((*events)[0][0] == expected1);
    BOOST_CHECK((*events)[0][1] == expected2);
    print_events(events);
    cout << "\n";

    // cheap division -> expect division and disappearance
    div.theEnergy = 0;
    mov.theEnergy = 200;
    mismov.theEnergy = 10;
    det.theEnergy = 5;
    misdet.theEnergy = 50;

    gm = GraphicalModel::from(traxels, div,mov,mismov,det,misdet);
    gm->infer();
    events = gm->events();
    BOOST_CHECK_EQUAL(events->size(), 1);
    BOOST_CHECK_EQUAL((*events)[0].size(), 2);
    expected1.type = Event::Division;
    expected1.traxel_ids.clear();
    expected1.traxel_ids.push_back(tr11.Id);
    expected1.traxel_ids.push_back(tr21.Id);
    expected1.traxel_ids.push_back(tr22.Id);
    expected2.type = Event::Disappearance;
    expected2.traxel_ids.clear();
    expected2.traxel_ids.push_back(tr12.Id);
    BOOST_CHECK((*events)[0][0] == expected1);
    BOOST_CHECK((*events)[0][1] == expected2);
    print_events(events);
    cout << "\n";
}

BOOST_AUTO_TEST_CASE( consecutive_moves ) {
    //
    // tr11 --- tr21 --- tr31
    //

    // trXY: traxel Y at timestep X
  Traxel tr11, tr21, tr31;
    feature_array com11(feature_array::difference_type(3));
    feature_array com21(feature_array::difference_type(3));
    feature_array com31(feature_array::difference_type(3));

    com11[0] = 0;
    com11[1] = 0;
    com11[2] = 0;
    tr11.features["com"] = com11;
    tr11.Id = 5;

    com21[0] = 1;
    com21[1] = 0;
    com21[2] = 0;
    tr21.features["com"] = com21;
    tr21.Id = 7;

    com31[0] = 2;
    com31[1] = 0;
    com31[2] = 0;
    tr31.features["com"] = com31;
    tr31.Id = 9;

    vector< vector<Traxel> > traxels;
    traxels.push_back(vector<Traxel>());
    traxels.push_back(vector<Traxel>());
    traxels.push_back(vector<Traxel>());
    traxels[0].push_back(tr11);
    traxels[1].push_back(tr21);
    traxels[2].push_back(tr31);

    //----
    // objects for reuse
    shared_ptr<GraphicalModel> gm;
    ConstantEnergy div, mov, mismov, det, misdet;
    Event expected1, expected2, expected3, expected4;
    shared_ptr<std::vector< std::vector<Event> > > events;

    // cheap detections, cheap moves -> expect two moves
    mov.theEnergy = 0;
    mismov.theEnergy = 10;
    det.theEnergy = 0;
    misdet.theEnergy = 10;
    gm = GraphicalModel::from(traxels, div,mov,mismov,det,misdet);
    gm->infer();
    events = gm->events();
    BOOST_CHECK_EQUAL(events->size(), 2);
    BOOST_CHECK_EQUAL((*events)[0].size(), 1);
    BOOST_CHECK_EQUAL((*events)[1].size(), 1);
    expected1.type = Event::Move;
    expected1.traxel_ids.clear();
    expected1.traxel_ids.push_back(tr11.Id);
    expected1.traxel_ids.push_back(tr21.Id);
    expected2.type = Event::Move;
    expected2.traxel_ids.clear();
    expected2.traxel_ids.push_back(tr21.Id);
    expected2.traxel_ids.push_back(tr31.Id);

    BOOST_CHECK((*events)[0][0] == expected1);
    BOOST_CHECK((*events)[1][0] == expected2);
    //BOOST_CHECK_EQUAL((*events)[0][0].energy, expected.energy);
    print_events(events);
    cout << "\n";

    // cheap detections, cheap mismoves -> expect disappearance and appearance
    mov.theEnergy = 10;
    mismov.theEnergy = 0;
    det.theEnergy = 0;
    misdet.theEnergy = 10;
    gm = GraphicalModel::from(traxels, div,mov,mismov,det,misdet);
    gm->infer();
    events = gm->events();
    BOOST_CHECK_EQUAL(events->size(), 2);
    BOOST_CHECK_EQUAL((*events)[0].size(), 2);
    BOOST_CHECK_EQUAL((*events)[1].size(), 2);
    expected1.type = Event::Disappearance;
    expected1.traxel_ids.clear();
    expected1.traxel_ids.push_back(tr11.Id);
    expected2.type = Event::Appearance;
    expected2.traxel_ids.clear();
    expected2.traxel_ids.push_back(tr21.Id);
    expected3.type = Event::Disappearance;
    expected3.traxel_ids.clear();
    expected3.traxel_ids.push_back(tr21.Id);
    expected4.type = Event::Appearance;
    expected4.traxel_ids.clear();
    expected4.traxel_ids.push_back(tr31.Id);

    BOOST_CHECK((*events)[0][0] == expected1);
    BOOST_CHECK((*events)[0][1] == expected2);
    BOOST_CHECK((*events)[1][0] == expected3);
    BOOST_CHECK((*events)[1][1] == expected4);

    print_events(events);
    cout << "\n";

    // cheap misdetections, cheap moves -> expect non-appearances
    mov.theEnergy = 0;
    mismov.theEnergy = 10;
    det.theEnergy = 10;
    misdet.theEnergy = 0;
    gm = GraphicalModel::from(traxels, div,mov,mismov,det,misdet);
    gm->infer();
    events = gm->events();
    BOOST_CHECK_EQUAL(events->size(), 2);
    BOOST_CHECK_EQUAL((*events)[0].size(), 2);
    BOOST_CHECK_EQUAL((*events)[1].size(), 2);
    expected1.type = Event::Disappearance;
    expected1.traxel_ids.clear();
    expected1.traxel_ids.push_back(tr11.Id);
    expected2.type = Event::Appearance;
    expected2.traxel_ids.clear();
    expected2.traxel_ids.push_back(tr21.Id);
    expected3.type = Event::Disappearance;
    expected3.traxel_ids.clear();
    expected3.traxel_ids.push_back(tr21.Id);
    expected4.type = Event::Appearance;
    expected4.traxel_ids.clear();
    expected4.traxel_ids.push_back(tr31.Id);

    BOOST_CHECK((*events)[0][0] == expected1);
    BOOST_CHECK((*events)[0][1] == expected2);
    BOOST_CHECK((*events)[1][0] == expected3);
    BOOST_CHECK((*events)[1][1] == expected4);
    print_events(events);
    cout << "\n";

  
}
// EOF

