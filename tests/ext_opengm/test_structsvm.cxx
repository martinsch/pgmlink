#include <vector>
#include <ostream>

#include <opengm/unittests/test.hxx>
#include <pgmlink/ext_opengm/loglinearmodel.hxx>
// #include <pgmlink/ext_opengm/structsvm.hxx>
#include <pgmlink/ext_opengm/loss_hamming.hxx>
#include "opengm/utilities/metaprogramming.hxx"

struct StructSvmTest {
   void testStructSvmDlib() {
     // typedef opengm::LoglinearModel<double,
     // 				    typename opengm::meta::TypeListGenerator<opengm::FunctionDecoratorWeighted<opengm::ExplicitFunction<double> >,
     // 									     opengm::HammingFunction<double> >::type
     //   > gm_type;
     // gm_type sample1(2), sample2(2);
     // sample1.addVariable(2);
     // sample1.addVariable(2);
     // sample2.addVariable(2);
     // sample2.addVariable(2);

     // size_t singlesite_shape[] = {2};
     // opengm::ExplicitFunction<double>* f1 = new opengm::ExplicitFunction<double>(singlesite_shape, singlesite_shape + 1, 2);
     // opengm::FunctionDecoratorWeighted<opengm::ExplicitFunction<double> > wf1(f1);

     // size_t pairwise_shape[] = {2, 2};
     // opengm::ExplicitFunction<double>* f2 = new opengm::ExplicitFunction<double>(pairwise_shape, pairwise_shape + 2, 3);
     // opengm::FunctionDecoratorWeighted<opengm::ExplicitFunction<double> > wf2(f2);

     // typename gm_type::FunctionIdentifier feat1, feat2;

     // wf1.innerFunction()->operator()(0) = 0.05;
     // wf1.innerFunction()->operator()(1) = 10;

     // wf2.innerFunction()->operator()(1, 1) = 3;
     // wf2.innerFunction()->operator()(1, 0) = 0.01;
     // wf2.innerFunction()->operator()(0, 1) = 0.01;
     // wf2.innerFunction()->operator()(0, 0) = 3;

     // feat1 = sample1.addFeature(wf1, 0);
     // feat2 = sample1.addFeature(wf2, 1);
     // {size_t idx[] = {0}; sample1.addFactor(feat1, idx, idx+1);}
     // {size_t idx[] = {1}; sample1.addFactor(feat1, idx, idx+1);}
     // {size_t idx[] = {0,1}; sample1.addFactor(feat2, idx, idx+2);}


     // wf1.innerFunction()->operator()(0) = 12;
     // wf1.innerFunction()->operator()(1) = 0.3;

     // // wf2.innerFunction()->operator()(1, 1) = 1;
     // // wf2.innerFunction()->operator()(1, 0) = 0.01;
     // // wf2.innerFunction()->operator()(0, 1) = 0.01;
     // // wf2.innerFunction()->operator()(0, 0) = 1.4;

     // feat1 = sample2.addFeature(wf1, 0);
     // feat2 = sample2.addFeature(wf2, 1);
     // {size_t idx[] = {0}; sample2.addFactor(feat1, idx, idx+1);}
     // {size_t idx[] = {1}; sample2.addFactor(feat1, idx, idx+1);}
     // {size_t idx[] = {0,1}; sample2.addFactor(feat2, idx, idx+2);}

     
     // std::vector<typename gm_type::LabelType> labels1;
     // labels1.push_back(1);
     // labels1.push_back(1);
     // std::vector<typename gm_type::LabelType> labels2;
     // labels2.push_back(0);
     // labels2.push_back(0);
     // std::vector<std::vector<typename gm_type::LabelType> > labels;
     // labels.push_back(labels1);
     // labels.push_back(labels2);

     // ///
     
     // std::vector<gm_type> samples;
     // samples.push_back(sample1); samples.push_back(sample2);
     // std::vector<typename gm_type::ValueType> weights;

     // opengm::StructSvmDlib< gm_type > structsvm(samples, labels);
     // structsvm.train(weights);
     // OPENGM_TEST_EQUAL(weights.size(), 2);
     
     // std::cout << weights[0] << "\n";
     // std::cout << weights[1] << "\n";
   }

   void run() {
      this->testStructSvmDlib();
   }
};

int main() {
  std::cout << "structsvm test... " << std::flush;
  StructSvmTest t;
  t.run();
  std::cout << "done." << std::endl;

  return 0;
}
