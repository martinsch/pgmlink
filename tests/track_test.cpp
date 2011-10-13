#define BOOST_TEST_MODULE track_test

#include <vector>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/shared_ptr.hpp>

#include "track.h"
#include "traxels.h"

using namespace boost;
using namespace Tracking;
using namespace std;

Event mk_event(Event::EventType type, vector<unsigned int> ids) {
    Event e;
    e.type = type;
    e.traxel_ids = ids;
    return e;
}

Traxel mk_traxel(double x, double y, double z, unsigned int id, int timestep = 0) {
   Traxel t;
    feature_array com(3);
    com[0] = x;
    com[1] = y;
    com[2] = z;
    t.Id = id;
    t.Timestep = timestep;
    t.features["com"] = com;
    return t;
}

void print_events(std::vector< std::vector<Event> > events) {
     for(size_t i = 0; i < events.size(); ++i){
	cout << "timestep " << i << ":\n";
	for(size_t j=0; j < events[i].size(); j++) {
	    cout << "\t" << events[i][j] << endl;
	}
    }
}



BOOST_AUTO_TEST_CASE( same_energy__different_objects )
{
  shared_ptr<const ConstantEnergy> div(new ConstantEnergy(99));
  shared_ptr<const ConstantEnergy> mov(new ConstantEnergy(100));
  shared_ptr<const ConstantEnergy> app(new ConstantEnergy(1000));
  shared_ptr<const ConstantEnergy> dis(new ConstantEnergy(1001));
  AdaptiveEnergiesFormulation f(div, mov, dis, app);

  Track m(f);
    Traxels prev, curr;
    vector<Event> events;

    // empty world
    events = m(prev, curr);
    BOOST_CHECK(events.empty());

    // a disappearance
    prev[27] = mk_traxel(10, 11.3, 9.0, 27);
    events = m(prev, curr);
    BOOST_CHECK_EQUAL(events.size(), 1);
    BOOST_CHECK_EQUAL(events[0].type, Event::Disappearance);
    BOOST_CHECK_EQUAL(events[0].traxel_ids[0], 27);
    BOOST_CHECK_EQUAL(events[0].energy, 1001);
    prev.clear();
    curr.clear();

    // an appearance
    curr[3] = mk_traxel(1, 0.5, 1.5, 3);
    events = m(prev, curr);
    BOOST_CHECK_EQUAL(events.size(), 1);
    BOOST_CHECK_EQUAL(events[0].type, Event::Appearance);
    BOOST_CHECK_EQUAL(events[0].traxel_ids[0], 3);
    BOOST_CHECK_EQUAL(events[0].energy, 1000);
    prev.clear();
    curr.clear();

    // a move
    prev[14] = mk_traxel(1,1,1,14);
    curr[17] = mk_traxel(2,1,1,17);
    events = m(prev, curr);
    BOOST_CHECK_EQUAL(events.size(), 1);
    BOOST_CHECK_EQUAL(events[0].type, Event::Move);
    BOOST_CHECK_EQUAL(events[0].traxel_ids[0], 14);
    BOOST_CHECK_EQUAL(events[0].traxel_ids[1], 17);
    BOOST_CHECK_EQUAL(events[0].energy, 100);
    prev.clear();
    curr.clear();

    // a division
    prev[14] = mk_traxel(1,1,1,14);
    curr[17] = mk_traxel(2,1,1,17);
    curr[18] = mk_traxel(0,1,1,18);
    events = m(prev, curr);
    BOOST_CHECK_EQUAL(events.size(), 1);
    BOOST_CHECK_EQUAL(events[0].type, Event::Division);
    BOOST_CHECK_EQUAL(events[0].traxel_ids[0], 14);
    BOOST_CHECK_EQUAL(events[0].traxel_ids[1], 17);
    BOOST_CHECK_EQUAL(events[0].traxel_ids[2], 18);
    BOOST_CHECK_EQUAL(events[0].energy, 99);
    prev.clear();
    curr.clear();
}

BOOST_AUTO_TEST_CASE( same_objects__different_energies )
{
  shared_ptr<ConstantEnergy> div(new ConstantEnergy);
  shared_ptr<ConstantEnergy> mov(new ConstantEnergy);
  shared_ptr<ConstantEnergy> app(new ConstantEnergy);
  shared_ptr<ConstantEnergy> dis(new ConstantEnergy);

  Track m;
  Traxels prev, curr;
  vector<Event> events;
  
  prev[1] = mk_traxel(0,0,0,1);
  prev[2] = mk_traxel(0,1,0,2);
  curr[3] = mk_traxel(1,0,0,3);
  curr[4] = mk_traxel(1,1,0,4);
  
  // only (dis-)appearances
  div->theEnergy = 1000;
  mov->theEnergy = 1000;
  app->theEnergy = 0;
  dis->theEnergy = 0;
  m.Formulation(AdaptiveEnergiesFormulation(div, mov, dis, app));
  
  events = m(prev, curr);
  BOOST_CHECK_EQUAL(events.size(), 4);
  BOOST_CHECK_EQUAL(events[0].type, Event::Disappearance);
  BOOST_CHECK_EQUAL(events[1].type, Event::Disappearance);
  BOOST_CHECK_EQUAL(events[2].type, Event::Appearance);
  BOOST_CHECK_EQUAL(events[3].type, Event::Appearance);

  // two moves
  div->theEnergy = 0;
  mov->theEnergy = 100;
  app->theEnergy = 100000;
  dis->theEnergy = 100000;
  m.Formulation(AdaptiveEnergiesFormulation(div, mov, dis, app));
  
  events = m(prev, curr);
  BOOST_CHECK_EQUAL(events.size(), 2);
  BOOST_CHECK_EQUAL(events[0].type, Event::Move);
  BOOST_CHECK_EQUAL(events[1].type, Event::Move);

  // one disappearance, one division
  div->theEnergy = 0;
  mov->theEnergy = 100;
  app->theEnergy = 10;
  dis->theEnergy = 10;
  m.Formulation(AdaptiveEnergiesFormulation(div, mov, dis, app));
  
  events = m(prev, curr);
  BOOST_CHECK_EQUAL(events.size(), 2);
  BOOST_CHECK_EQUAL(events[0].type, Event::Division);
  BOOST_CHECK_EQUAL(events[1].type, Event::Disappearance);
}

// EOF

