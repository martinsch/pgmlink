#pragma once
#ifndef OPENGM_DECORATORWEIGHTED_HXX
#define OPENGM_DECORATORWEIGHTED_HXX

#include <algorithm>
#include "opengm/opengm.hxx"
#include "opengm/functions/function_properties_base.hxx"

namespace opengm {

/// Make any function a weighted function.
///
/// \ingroup functions
template<class FUNCTION>
class FunctionDecoratorWeighted
  : public FunctionBase<FunctionDecoratorWeighted<FUNCTION>,
		      typename FUNCTION::ValueType,
		      typename FUNCTION::IndexType,
		      typename FUNCTION::LabelType>
{
public:
  typedef FUNCTION InnerFunctionType;
  typedef typename FUNCTION::ValueType ValueType;
  typedef typename FUNCTION::IndexType IndexType;
  typedef typename FUNCTION::LabelType LabelType;
  
  FunctionDecoratorWeighted( InnerFunctionType*,
			     ValueType weight=1.);
  FunctionDecoratorWeighted( const FunctionDecoratorWeighted<FUNCTION>& );
  FunctionDecoratorWeighted& operator=( FunctionDecoratorWeighted<FUNCTION> );
  template <class FUNC>
  friend void swap(FunctionDecoratorWeighted<FUNC>&, FunctionDecoratorWeighted<FUNC>&);
  ~FunctionDecoratorWeighted() { delete inner_; };
  
  template<class Iterator> ValueType operator()(Iterator begin) const;
  size_t dimension() const { return inner_->dimension(); };
  size_t shape(const IndexType i) const { return inner_->shape( i ); };
  size_t size() const { return inner_->size(); };
  
  void innerFunction( InnerFunctionType* f ) { delete inner_; inner_ = f; };
  InnerFunctionType* innerFunction() { return inner_; };
  const InnerFunctionType* innerFunction() const { return inner_; };
  void weight( const ValueType& w ) { weight_ = w; };
  ValueType weight() const { return weight_; };
   
private:
   InnerFunctionType* inner_;
   ValueType weight_;
};

/// Constructor
/// \param f function to be decorated (takes ownership of pointer)
/// \param weight the decorated function is multiplied by weight 
template <class FUNCTION>
FunctionDecoratorWeighted<FUNCTION>::FunctionDecoratorWeighted
(
 InnerFunctionType* f, ValueType weight
)
:  inner_(f), weight_(weight)
{
  if( f == NULL ) {
    throw RuntimeError("FunctionDecoratorWeighted(): inner function is NULL");
  }
}

template <class FUNCTION>
FunctionDecoratorWeighted<FUNCTION>::FunctionDecoratorWeighted( const FunctionDecoratorWeighted<FUNCTION>& other )
  : weight_(other.weight_) {
  inner_ = new InnerFunctionType(*other.inner_);
}

template <class FUNCTION>
FunctionDecoratorWeighted<FUNCTION>& FunctionDecoratorWeighted<FUNCTION>::operator=( FunctionDecoratorWeighted<FUNCTION> other ) {
  swap(*this, other);
  return *this;
}

template <class FUNC>
void swap(FunctionDecoratorWeighted<FUNC>& lhs, FunctionDecoratorWeighted<FUNC>& rhs) {
    using std::swap; // enable ADL
    swap(lhs.inner_, rhs.inner_);
    swap(lhs.weight_, rhs.weight_);
}

template <class FUNCTION>
template<class Iterator>
inline typename FunctionDecoratorWeighted<FUNCTION>::ValueType
FunctionDecoratorWeighted<FUNCTION>::operator()
(
Iterator begin
) const {
  return weight_ * (*inner_)(begin);
}

} // namespace opengm

#endif // #ifndef OPENGM_DECORATORWEIGHTED_HXX
