
#include <opengm/unittests/test.hxx>
#include <pgmlink/ext_opengm/feature.hxx>
#include <opengm/functions/explicit_function.hxx>

struct FeatureTest {
  void testWeightObserving() {
    typedef opengm::ExplicitFunction<double> Function;
    using opengm::learning::WeightObserver;
    using opengm::learning::ObservableWeight;


    ObservableWeight<Function>* w = new ObservableWeight<Function>( 2.3 );

    WeightObserver<Function>* o0 = new WeightObserver<Function>(new Function());
    WeightObserver<Function>* o1 = new WeightObserver<Function>(new Function(), w);
    o1->weight(1);

    // value observation
    OPENGM_TEST( !w->hasObserver(o0) );
    OPENGM_TEST( w->hasObserver(o1) );
    
    OPENGM_TEST_EQUAL( w->weight(), 2.3);
    OPENGM_TEST_EQUAL( o1->weight(), 1);
    w->notifyObservers();
    OPENGM_TEST_EQUAL( o1->weight(), 2.3);
    
    w->weight(4);
    OPENGM_TEST_EQUAL( o1->weight(), 4);    
    OPENGM_TEST_EQUAL( w->weight(), 4);    

    // (de)registration
    w->registerObserver(o0);
    OPENGM_TEST( w->hasObserver(o0) );    
    OPENGM_TEST( o0->observable() == NULL);
    w->removeObserver(o0);
    OPENGM_TEST( !w->hasObserver(o0) );
    
    o0->observable(w);
    OPENGM_TEST( w->hasObserver(o0) );    
    OPENGM_TEST_EQUAL( o0->observable(), w );

    o0->observable(NULL);
    OPENGM_TEST( w->hasObserver(o0) );    
    OPENGM_TEST( o0->observable() == NULL );

    // test random destruction
    delete o0;
    delete w;
    delete o1;
  }

  void run() {
    testWeightObserving();
  }
};

int main() {
   std::cout << "Feature test...  " << std::endl;
   {
      FeatureTest t;
      t.run();
   }
   std::cout << "done.." << std::endl;
   return 0;
}
