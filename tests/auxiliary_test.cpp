#define BOOST_TEST_MODULE multi_hypotheses_test

// stl
#include <vector>
#include <iostream>

// boost
#include <boost/test/unit_test.hpp>

// vigra
#include <vigra/multi_array.hxx>

// pgmlink
#include "pgmlink/auxiliary.h"

using namespace pgmlink;


BOOST_AUTO_TEST_CASE( auxiliary_extract_intersects_3D ) {

  std::cout << "Testing extract_intersects for 3D arrays" << std::endl;

  vigra::Shape3 shape(3,3,1);
  vigra::MultiArray<3, unsigned> arr1(shape);
  vigra::MultiArray<3, unsigned> arr2(shape);
  
  arr1(0,0,0) = 1;
  arr1(1,0,0) = 1;
  arr1(2,0,0) = 1;
  arr1(0,1,0) = 1;
  arr1(1,1,0) = 2;
  arr1(2,1,0) = 2;
  arr1(0,2,0) = 0;
  arr1(1,2,0) = 2;
  arr1(2,2,0) = 2;

  arr2(0,0,0) = 1;
  arr2(1,0,0) = 1;
  arr2(2,0,0) = 2;
  arr2(0,1,0) = 2;
  arr2(1,1,0) = 3;
  arr2(2,1,0) = 3;
  arr2(0,2,0) = 4;
  arr2(1,2,0) = 4;
  arr2(2,2,0) = 5;

  std::map<unsigned, std::map<unsigned, unsigned> > intersects;

  extract_intersects<3, unsigned>(arr1, arr2, intersects);
  BOOST_CHECK_EQUAL(intersects[1][1], 2);
  BOOST_CHECK_EQUAL(intersects[1][2], 2);
  BOOST_CHECK_EQUAL(intersects[2][3], 2);
  BOOST_CHECK_EQUAL(intersects[2][4], 1);
  BOOST_CHECK_EQUAL(intersects[2][5], 1);
}


BOOST_AUTO_TEST_CASE( auxiliary_extract_intersects_2D ) {

  std::cout << "Testing extract_intersects for 2D arrays" << std::endl;

  vigra::Shape2 shape(3,3);
  vigra::MultiArray<2, unsigned> arr1(shape);
  vigra::MultiArray<2, unsigned> arr2(shape);
  
  arr1(0,0) = 1;
  arr1(1,0) = 1;
  arr1(2,0) = 1;
  arr1(0,1) = 1;
  arr1(1,1) = 2;
  arr1(2,1) = 2;
  arr1(0,2) = 0;
  arr1(1,2) = 2;
  arr1(2,2) = 2;

  arr2(0,0) = 1;
  arr2(1,0) = 1;
  arr2(2,0) = 2;
  arr2(0,1) = 2;
  arr2(1,1) = 3;
  arr2(2,1) = 3;
  arr2(0,2) = 4;
  arr2(1,2) = 4;
  arr2(2,2) = 5;

  std::map<unsigned, std::map<unsigned, unsigned> > intersects;

  extract_intersects<2, unsigned>(arr1, arr2, intersects);
  BOOST_CHECK_EQUAL(intersects[1][1], 2);
  BOOST_CHECK_EQUAL(intersects[1][2], 2);
  BOOST_CHECK_EQUAL(intersects[2][3], 2);
  BOOST_CHECK_EQUAL(intersects[2][4], 1);
  BOOST_CHECK_EQUAL(intersects[2][5], 1);
}
