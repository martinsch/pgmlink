#pragma once
#ifndef OPENGM_LOSS_HAMMING_HXX
#define OPENGM_LOSS_HAMMING_HXX

#include <cassert>
#include <vector>
#include <sstream>

#include "opengm/opengm.hxx"
#include "opengm/functions/function_properties_base.hxx"

namespace opengm {

template<class T, class I = size_t, class L = size_t>
class HammingFunction
  : public FunctionBase<HammingFunction<T, I, L>, T, I, L> {
public:
  typedef T ValueType;
  typedef I IndexType;
  typedef L LabelType;

  HammingFunction() : trueLabel_(0), numberOfLabels_(0) {}
  HammingFunction( LabelType numberOfLabels, LabelType trueLabel ) : trueLabel_(trueLabel), numberOfLabels_(numberOfLabels ) {}

  template<class Iterator>
  ValueType operator()(Iterator) const; 
  size_t shape(const size_t) const;
  size_t dimension() const { return 1; }
  size_t size() const { return numberOfLabels_; }

private:
  L trueLabel_;
  L numberOfLabels_;
};

////////////////////
/* implementation */
////////////////////
template<class T, class I, class L>
template<class Iterator>
typename HammingFunction<T,I,L>::ValueType HammingFunction<T,I,L>::operator()(Iterator it) const {
  return (*it == trueLabel_ ? 0. : -1.) ;
}

template<class T, class I, class L>
size_t HammingFunction<T,I,L>::shape(const size_t idx) const { 
  if(idx > 1) {
    std::stringstream ss;
    ss << "HammingFunction::shape(): argument has to be 1 ( was " << idx << ")";
    throw RuntimeError(ss.str());
  }
  return numberOfLabels_; 
}


} // end namespace opengm

#endif
