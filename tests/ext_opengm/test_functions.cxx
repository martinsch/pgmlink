#include "opengm/functions/constant.hxx"
#include "pgmlink/ext_opengm/decorator_weighted.hxx"
#include "pgmlink/ext_opengm/indicator_function.hxx"

#include <opengm/unittests/test.hxx>

template<class T>
struct FunctionsTest {
  template<class FUNCTION>
  void testProperties(const FUNCTION & f) {
    opengm::ShapeWalker<typename FUNCTION::FunctionShapeIteratorType > walker(f.functionShapeBegin(),f.dimension());
    for(size_t i=0;i<f.dimension();++i) {
      OPENGM_TEST(f.functionShapeBegin()[i] == f.shape(i));   
    }
      
    typedef typename  FUNCTION::ValueType ValueType;
    ValueType min=f(walker.coordinateTuple().begin());
    ValueType max=f(walker.coordinateTuple().begin());
    ValueType sum=static_cast<ValueType>(0);
    ValueType prod=static_cast<ValueType>(1);
    for(size_t i=0;i<f.size();++i) {
      ValueType tmp=f(walker.coordinateTuple().begin());
      min=tmp<min?tmp:min;
      max=tmp>max?tmp:max;
      sum+=tmp;
      prod*=tmp;
      ++walker;
    }
    OPENGM_TEST_EQUAL_TOLERANCE(min,f.min(),static_cast<ValueType>(0.0001));
    OPENGM_TEST_EQUAL_TOLERANCE(max,f.max(),static_cast<ValueType>(0.0001));
    OPENGM_TEST_EQUAL_TOLERANCE(f.minMax().max(),f.max(),static_cast<ValueType>(0.0001));
    OPENGM_TEST_EQUAL_TOLERANCE(f.minMax().min(),f.min(),static_cast<ValueType>(0.0001));
    OPENGM_TEST_EQUAL_TOLERANCE(sum,f.sum(),static_cast<ValueType>(0.0001));
    OPENGM_TEST_EQUAL_TOLERANCE(prod,f.product(),static_cast<ValueType>(0.0001));
  }

  void testIndicatorFunction() {
      std::cout << "  * IndicatorFunction" << std::endl;
      size_t shape[]={4,4,4};
      size_t indicate[]={1,2,3};
      opengm::IndicatorFunction<T> f(shape, shape+3, indicate, 17, 2);
      OPENGM_TEST_EQUAL(f.dimension(), 3);
      OPENGM_TEST_EQUAL(f.shape(0), 4);
      OPENGM_TEST_EQUAL(f.shape(1), 4);
      OPENGM_TEST_EQUAL(f.shape(2), 4);
      OPENGM_TEST_EQUAL(f.size(), 4*4*4);

      size_t arg2[]={1,1,3};
      OPENGM_TEST_EQUAL(f(indicate), 17);
      OPENGM_TEST_EQUAL(f(arg2), 2);
  }

  void testFunctionDecoratorWeighted() {
    std::cout << "  * FunctionDecoratorWeighted" << std::endl;
    double constant = 3;
    size_t shape[] = {2,3,1};
    size_t index[] = {1,1,0};

    opengm::ConstantFunction<double>* decorated = new opengm::ConstantFunction<double>(shape, shape+3, constant);
    opengm::FunctionDecoratorWeighted<opengm::ConstantFunction<double> > f =
      opengm::FunctionDecoratorWeighted<opengm::ConstantFunction<double> >( decorated);

    OPENGM_TEST_EQUAL(f.innerFunction(), decorated);
    OPENGM_TEST_EQUAL(f.weight(), 1.);
    OPENGM_TEST_EQUAL(f(index), constant);

    f.weight(2.);
    OPENGM_TEST_EQUAL(f.weight(), 2.);
    OPENGM_TEST_EQUAL(f(index), 2.*constant);     

    opengm::ConstantFunction<double>* another = new opengm::ConstantFunction<double>(shape, shape+3, constant);
    f.innerFunction( another );
    OPENGM_TEST(f.innerFunction() != decorated );     
    OPENGM_TEST_EQUAL(f.innerFunction(), another);     

    testProperties( f );
    decorated = NULL;

    //
    // test copy, assignment, and swap
    //
    double constant_other = 4;
    opengm::ConstantFunction<double>* decorated_other = new opengm::ConstantFunction<double>(shape, shape+3, constant_other);
    opengm::FunctionDecoratorWeighted<opengm::ConstantFunction<double> > f_other =
      opengm::FunctionDecoratorWeighted<opengm::ConstantFunction<double> >( decorated_other );

    // copy
    f.weight(3.);
    opengm::FunctionDecoratorWeighted<opengm::ConstantFunction<double> > f_copy =
      opengm::FunctionDecoratorWeighted<opengm::ConstantFunction<double> >( f );
    OPENGM_TEST_EQUAL( f_copy.weight(), f.weight() );
    OPENGM_ASSERT( f_copy.innerFunction() != f.innerFunction() );
    OPENGM_TEST_EQUAL( f_copy(index), f(index) );

    // swap
    opengm::swap(f, f_other);
    OPENGM_TEST_EQUAL( f_copy.weight(), f_other.weight() );
    OPENGM_ASSERT( f_copy.innerFunction() != f_other.innerFunction() );
    OPENGM_TEST_EQUAL( f_copy(index), f_other(index) );

    OPENGM_ASSERT( f_copy.weight() != f.weight() );
    OPENGM_ASSERT( f_copy.innerFunction() != f.innerFunction() );
    OPENGM_ASSERT( f_copy(index) != f(index) );

    opengm::swap(f, f_other);
    OPENGM_TEST_EQUAL( f_copy.weight(), f.weight() );
    OPENGM_ASSERT( f_copy.innerFunction() != f.innerFunction() );
    OPENGM_TEST_EQUAL( f_copy(index), f(index) );

    // assignment
    OPENGM_ASSERT( f.weight() != f_other.weight() );
    OPENGM_ASSERT( f.innerFunction() != f_other.innerFunction() );
    OPENGM_ASSERT( f(index) != f_other(index) );
    f = f_other;
    OPENGM_TEST_EQUAL( f.weight(), f_other.weight() );
    OPENGM_ASSERT( f.innerFunction() != f_other.innerFunction() );
    OPENGM_TEST_EQUAL( f(index), f_other(index) );
  }

   void run() {
      testIndicatorFunction();
      testFunctionDecoratorWeighted();
   }
};

int main() {
   std::cout << "Functions test...  " << std::endl;
   {
      FunctionsTest<int >t;
      t.run();
   }
   std::cout << "done.." << std::endl;
   return 0;
}
