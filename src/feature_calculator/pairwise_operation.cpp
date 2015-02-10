// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/feature_calculator/base.h"
#include "pgmlink/feature_calculator/pairwise_operation.h"

namespace pgmlink {

namespace feature_extraction {

PairwiseOperationCalculator::PairwiseOperationCalculator(PairwiseOperationCalculator::Operation op, const std::string& name):
	name_(name),
	operation_(op)
{}

PairwiseOperationCalculator::~PairwiseOperationCalculator() {}

feature_array PairwiseOperationCalculator::calculate(const feature_array& f1, const feature_array& f2) const {
  return operation_(f1, f2);
}

const std::string& PairwiseOperationCalculator::name() const {
  return name_;
}

TripletOperationCalculator::TripletOperationCalculator(TripletOperationCalculator::Operation op, const std::string& name):
  name_(name),
  operation_(op)
{}

TripletOperationCalculator::~TripletOperationCalculator() {}

feature_array TripletOperationCalculator::calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const {
  return operation_(f1, f2, f3);
}

const std::string& TripletOperationCalculator::name() const {
  return name_;
}

} // namespace feature_extraction

} // namespace pgmlink
