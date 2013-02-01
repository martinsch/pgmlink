#pragma once
#ifndef OPENGM_INDICATORFUNCTION_HXX
#define OPENGM_INDICATORFUNCTION_HXX

#include <algorithm>
#include <numeric>
#include <vector>
#include "opengm/opengm.hxx"
#include "opengm/functions/function_properties_base.hxx"

namespace opengm {

/// Indicator Function
///
/// The indicator function assumes two values: a default value
/// (usually 0) for all arguments other than the indicated argument
/// and a another value for the indicated argument (usually 1).
///
/// \ingroup functions
template<class T, class I = size_t, class L = size_t>
class IndicatorFunction
  : public FunctionBase<IndicatorFunction<T, I, L>, T, I, L>
{
public:
  typedef T ValueType;
  typedef I IndexType;
  typedef L LabelType;
  
  template<class IT1, class IT2>
  IndicatorFunction(IT1 shapeBegin, IT1 shapeEnd,
		    IT2 indicate,
		    ValueType indicate_value = 1., ValueType default_value = 0.); 

  template<class ITERATOR> ValueType operator()(ITERATOR) const;  

  size_t dimension() const { return shape_.size(); }
  size_t shape(const IndexType i) const { return shape_[i]; }
  size_t size() const { if(shape_.size()){ return std::accumulate(shape_.begin(), shape_.end(), 1, Multiplier());} else return 0; }
  
  const std::vector<LabelType>& indicate() const { return indicate_; };
  void indicate_value( const ValueType& v ) { value_ = v; }
  ValueType indicate_value() const { return value_; }
  void default_value( const ValueType& v ) { default_ = v; }
  ValueType default_value() const { return default_; }
   
private:
  struct Multiplier {
    ValueType operator()( ValueType a, ValueType b ) { return a*b; }
  };

  std::vector<LabelType> shape_;
  std::vector<LabelType> indicate_;
  ValueType value_;
  ValueType default_;
};



/**/
/* implementation */
/**/
  template<class T, class I, class L>
  template<class IT1, class IT2>
  IndicatorFunction<T,I,L>::IndicatorFunction(IT1 shapeBegin, IT1 shapeEnd,
					      IT2 indicate,
					      ValueType indicate_value, ValueType default_value)
    : value_(indicate_value), default_(default_value) {
    for(; shapeBegin!=shapeEnd; ++shapeBegin, ++indicate) {
      shape_.push_back(*shapeBegin);
      indicate_.push_back(*indicate);
    }
  }

  template<class T, class I, class L>
  template<class Iterator>
  inline typename IndicatorFunction<T,I,L>::ValueType IndicatorFunction<T,I,L>::operator()(Iterator begin) const {
    return std::equal(indicate_.begin(), indicate_.end(), begin) ? value_ : default_;
  }
} // namespace opengm

#endif // #ifndef OPENGM_INDICATORFUNCTION_HXX
