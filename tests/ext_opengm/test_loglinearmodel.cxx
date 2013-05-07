#include <vector>

#include <pgmlink/ext_opengm/loglinearmodel.hxx>
#include <opengm/operations/adder.hxx>
#include <opengm/unittests/test.hxx>
#include <opengm/functions/explicit_function.hxx>
#include <pgmlink/ext_opengm/decorator_weighted.hxx>
#include <opengm/utilities/metaprogramming.hxx>

template<class T, class I, class L>
struct LoglinearModelTest {
   void testAccessingWeights() {
     typedef opengm::LoglinearModel<
       T,
       opengm::meta::TypeList<opengm::ExplicitFunction<T>,
			      opengm::meta::TypeList<opengm::FunctionDecoratorWeighted<opengm::ExplicitFunction<T> >,
						     opengm::meta::ListEnd>
				        >
       > gm_type;
     gm_type gm;
     opengm::ExplicitFunction<T>* f = new opengm::ExplicitFunction<T>(5);
     opengm::FunctionDecoratorWeighted<opengm::ExplicitFunction<T> > wf(f);
     opengm::ExplicitFunction<T> ef(5);
     wf.weight(2.);
     OPENGM_TEST_EQUAL(gm.numberOfWeights(), 0);
     gm.incrementNumberOfWeights(2);
     gm.addFeature(wf, 0);
     OPENGM_TEST_EQUAL(gm.numberOfWeights(), 1);
     gm.incrementNumberOfWeights(3);
     wf.weight(7.); // should be overriden by the incrementNumberOfWeights assignment

     // add to same weight index
     wf.weight(4.);
     gm.addFeature(wf, 1);

     // should be ignored by the loglinear special interface
     // (both weighted and other types)
     gm.addFunction(ef);
     gm.addFunction(wf);

     std::vector<T> weights;
     gm.getWeights( weights );
     OPENGM_TEST_EQUAL(weights.size(), 2);
     OPENGM_TEST_EQUAL(weights[0], 2);
     OPENGM_TEST_EQUAL(weights[1], 3);

     weights = std::vector<T>(2);
     weights[0] = 20;
     weights[1] = 30;
     gm.setWeights( weights );
     weights.clear();
     gm.getWeights( weights );
     OPENGM_TEST_EQUAL(weights.size(), 2);
     OPENGM_TEST_EQUAL(weights[0], 20);
     OPENGM_TEST_EQUAL(weights[1], 30);
   }

   void testFeatureValues() {
     typedef opengm::LoglinearModel<
       T,
       opengm::meta::TypeList<opengm::ExplicitFunction<T>,
			      opengm::meta::TypeList<opengm::FunctionDecoratorWeighted<opengm::ExplicitFunction<T> >,
						     opengm::meta::ListEnd>
				        >
       > gm_type;
     gm_type gm(2);

     size_t singlesite_shape[] = {2};
     opengm::ExplicitFunction<T>* f1 = new opengm::ExplicitFunction<T>(singlesite_shape, singlesite_shape + 1, 2);
     opengm::FunctionDecoratorWeighted<opengm::ExplicitFunction<T> > wf1(f1);
     typename gm_type::FunctionIdentifier feat1 = gm.addFeature(wf1, 0);

     size_t pairwise_shape[] = {2, 2};
     opengm::ExplicitFunction<T>* f2 = new opengm::ExplicitFunction<T>(pairwise_shape, pairwise_shape + 2, 3);
     opengm::FunctionDecoratorWeighted<opengm::ExplicitFunction<T> > wf2(f2);
     typename gm_type::FunctionIdentifier feat2 = gm.addFeature(wf2, 1);

     std::vector<typename gm_type::IndexType> idx;
     idx.push_back(gm.addVariable(2));
     idx.push_back(gm.addVariable(2));
     
     gm.addFactor(feat1, idx.begin(), idx.begin() + 1 );
     gm.addFactor(feat1, idx.begin() + 1, idx.end() );
     gm.addFactor(feat2, idx.begin(), idx.end() );

     std::vector<typename gm_type::LabelType> labels(2, 0);
     std::vector<T> featVals(2, 0);
     gm.weightedFeatureSums( labels, featVals );
     OPENGM_TEST_EQUAL(featVals[0], 4);
     OPENGM_TEST_EQUAL(featVals[1], 3);
     OPENGM_TEST_EQUAL(featVals.size(), 2);
   }


   void run() {
      this->testAccessingWeights();
      this->testFeatureValues();
   }
};

int main() {
   {
      std::cout << "LoglinearModel test... " << std::flush;
      {
         LoglinearModelTest<float, unsigned int, unsigned short> t;
         t.run();
      }
      {
         LoglinearModelTest<double, unsigned long, unsigned short> t;
         t.run();
      }
      {
         LoglinearModelTest<double, unsigned short, unsigned short> t;
         t.run();
      }
      {
         LoglinearModelTest<double, size_t, size_t> t;
         t.run();
      }
      std::cout << "done." << std::endl;
   }

   return 0;
}
