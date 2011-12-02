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
using namespace opengm;

BOOST_AUTO_TEST_CASE( ogm_model_building )
{
   OpengmMrf::ogmGraphicalModel gm;
   gm.addVariable(2);
   gm.addVariable(2);
   gm.addVariable(2);

   size_t vv[] = {1,2,3,4,5,6,7,8};
   for(size_t variable1 = 0; variable1 < gm.numberOfVariables(); ++variable1)
   for(size_t variable2 = variable1 + 1; variable2 < gm.numberOfVariables(); ++variable2)
   for(size_t variable3 = variable2 + 1; variable3 < gm.numberOfVariables(); ++variable3) {
      const size_t shape[] = {
         gm.numberOfStates(variable1),
         gm.numberOfStates(variable2),
         gm.numberOfStates(variable3)
      };
      OpengmMrf::ExplicitFunctionType f(shape, shape + 3);
      for(size_t state1 = 0; state1 < gm.numberOfStates(variable1); ++state1)
      for(size_t state2 = 0; state2 < gm.numberOfStates(variable2); ++state2)
      for(size_t state3 = 0; state3 < gm.numberOfStates(variable3); ++state3) {
         f(0,0,0)=float(vv[0]);
	 f(1,0,0)=float(vv[1]);
	 f(0,1,0)=float(vv[2]);
	 f(1,1,0)=float(vv[3]);
	 f(0,0,1)=float(vv[4]);
	 f(1,0,1)=float(vv[5]);
	 f(0,1,1)=float(vv[6]);
	 f(1,1,1)=float(vv[7]);
      }
      OpengmMrf::FunctionIdentifier id = gm.addFunction(f);
      size_t variableIndexSequence[] = {variable1, variable2, variable3};
      gm.addFactor(id, variableIndexSequence, variableIndexSequence + 3);
      BOOST_CHECK_EQUAL(f(0,0,0), 1);
      BOOST_CHECK_EQUAL(f(1,0,0), 2);
      BOOST_CHECK_EQUAL(f(0,1,0), 3);
      BOOST_CHECK_EQUAL(f(1,1,0), 4);
      BOOST_CHECK_EQUAL(f(0,0,1), 5);
      BOOST_CHECK_EQUAL(f(1,0,1), 6);
      BOOST_CHECK_EQUAL(f(0,1,1), 7);
      BOOST_CHECK_EQUAL(f(1,1,1), 8);
   }
}
// EOF

