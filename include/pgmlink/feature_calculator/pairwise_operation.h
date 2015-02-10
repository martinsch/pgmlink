#ifndef PAIRWISE_OPERATION_H
#define PAIRWISE_OPERATION_H

#include <boost/function.hpp>

#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator/base.h"
#include "pgmlink/feature_extraction.h"

namespace pgmlink {

namespace feature_extraction {

////
//// class PairwiseOperationCalculator
////
class PairwiseOperationCalculator : public FeatureCalculator {
 public:
  typedef boost::function<feature_array(const feature_array&, const feature_array&)> Operation;

  PairwiseOperationCalculator(Operation op, const std::string& name );
  virtual ~PairwiseOperationCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual const std::string& name() const;

 private:
  std::string name_;
  Operation operation_;
};

////
//// class TripletOperationCalculator
////
class TripletOperationCalculator : public FeatureCalculator {
 public:
  typedef boost::function<feature_array(const feature_array&, const feature_array&, const feature_array&)> Operation;

  TripletOperationCalculator(Operation op, const std::string& name );
  virtual ~TripletOperationCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const;
  virtual const std::string& name() const;

 private:
  std::string name_;
  Operation operation_;
};

} // namespace feature_extraction

} // namespace pgmlink

# endif // PAIRWISE_OPERATION_H
