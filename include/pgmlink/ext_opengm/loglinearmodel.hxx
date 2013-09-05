#pragma once
#ifndef OPENGM_LOGLINEARMODEL_HXX
#define OPENGM_LOGLINEARMODEL_HXX

#include <algorithm>
#include <cassert>
#include <map>
#include <stdexcept>
#include <utility>

#include "opengm/opengm.hxx"
#include "opengm/operations/adder.hxx"
#include "pgmlink/ext_opengm/decorator_weighted.hxx"
#include "opengm/graphicalmodel/graphicalmodel.hxx"
#include "opengm/utilities/metaprogramming.hxx"

namespace opengm {
/// \cond HIDDEN_SYMBOLS
namespace detail_graphical_model {
  template<size_t FUNCTION_INDEX, class FUNCTION_TYPE, class T, class FUNCTION_TYPE_LIST, class SPACE>
  struct WeightAccessor;
  template<size_t FUNCTION_TYPE_INDEX, class FUNCTION_TYPE, class LLM>
  struct InnerFunctionAccessor;
}
/// \endcond

/// \brief LoglinearModel
///
/// \ingroup graphical_models
template<
   class T,
   class FUNCTION_TYPE_LIST = meta::TypeList<FunctionDecoratorWeighted<ExplicitFunction<T> >, meta::ListEnd>, 
   class SPACE = opengm::DiscreteSpace<size_t, size_t>
>
class LoglinearModel 
  : public GraphicalModel<T, Adder, FUNCTION_TYPE_LIST, SPACE> 
{
public:
    typedef LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE> GraphicalModelType;
    typedef typename GraphicalModel<T, Adder, FUNCTION_TYPE_LIST, SPACE>::SpaceType SpaceType;
    typedef typename GraphicalModel<T, Adder, FUNCTION_TYPE_LIST, SPACE>::IndexType IndexType;
    typedef typename GraphicalModel<T, Adder, FUNCTION_TYPE_LIST, SPACE>::LabelType LabelType;
    typedef typename GraphicalModel<T, Adder, FUNCTION_TYPE_LIST, SPACE>::ValueType ValueType;
       
    typedef typename GraphicalModel<T, Adder, FUNCTION_TYPE_LIST, SPACE>::FunctionTypeList FunctionTypeList;
       
    typedef typename GraphicalModel<T, Adder, FUNCTION_TYPE_LIST, SPACE>::FunctionIdentifier FunctionIdentifier;
    typedef typename GraphicalModel<T, Adder, FUNCTION_TYPE_LIST, SPACE>::IndependentFactorType IndependentFactorType; 
    typedef typename GraphicalModel<T, Adder, FUNCTION_TYPE_LIST, SPACE>::FactorType FactorType;
    typedef typename GraphicalModel<T, Adder, FUNCTION_TYPE_LIST, SPACE>::OperatorType OperatorType;

    LoglinearModel( IndexType numberOfWeights = 0, ValueType val = 1 ) : weights_(numberOfWeights, val) {}
    template<class FUNCTION_TYPE_LIST_OTHER>
    LoglinearModel(const LoglinearModel<T, FUNCTION_TYPE_LIST_OTHER, SPACE>& m,
		   IndexType numberOfWeights = 0,
		   ValueType val = 1)
      : GraphicalModel<T, Adder, FUNCTION_TYPE_LIST, SPACE>::GraphicalModel(m), weights_(numberOfWeights, val) {}
    LoglinearModel(const SpaceType& s,
		   IndexType numberOfWeights = 0,
		   ValueType val = 1)
      : GraphicalModel<T, Adder, FUNCTION_TYPE_LIST, SPACE>::GraphicalModel(s), weights_(numberOfWeights, val) {}

    template<class FEATURE_FUNCTION>
      FunctionIdentifier
      addFeature(FunctionDecoratorWeighted<FEATURE_FUNCTION>, IndexType weightIndex );

    void setWeights( std::vector<ValueType> in );
    void getWeights( std::vector<ValueType>& out ) const;
    IndexType numberOfWeights() const { return static_cast<IndexType>(weights_.size()); }
    IndexType incrementNumberOfWeights( ValueType val=1 ) { weights_.push_back(val); return numberOfWeights() - 1; }
    IndexType increaseNumberOfWeights( IndexType by_n=1 ) { weights_.resize(weights_.size() + by_n, 1.); return numberOfWeights(); }

    void weightedFeatureSums( const std::vector<LabelType>& labelIndices, std::vector<ValueType>& out ) const;

private:
    std::vector<ValueType> weights_;
    std::map<FunctionIdentifier, IndexType> feature_ids_; // map function id to weight index


    template<size_t FUNCTION_INDEX, class FUNCTION_TYPE, class t, class function_type_list, class space>
    friend struct detail_graphical_model::WeightAccessor;

    template<size_t FUNCTION_TYPE_INDEX, class FUNCTION_TYPE, class LLM>
    friend struct detail_graphical_model::InnerFunctionAccessor;
};

/// \cond HIDDEN_SYMBOLS
namespace detail_graphical_model {
////
//// Accessing weights via TMP recursion
////

template<size_t FUNCTION_INDEX, class FUNCTION_TYPE, class T, class FUNCTION_TYPE_LIST, class SPACE>
struct WeightAccessor {
  typedef LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE> model_t;
  void get( const model_t* m, std::vector<T>& out ) const {
    w_next_.get(m, out);
  };
  void set( model_t* m, std::vector<T>& in ) {
    w_next_.set(m, in);
  };
  WeightAccessor<FUNCTION_INDEX-1, typename meta::TypeAtTypeList<FUNCTION_TYPE_LIST, FUNCTION_INDEX-1>::type, T, FUNCTION_TYPE_LIST, SPACE> w_next_;
};

template<size_t FUNCTION_INDEX, class FUNCTION_TYPE, class T, class FUNCTION_TYPE_LIST, class SPACE>
struct WeightAccessor<FUNCTION_INDEX, opengm::FunctionDecoratorWeighted<FUNCTION_TYPE>, T, FUNCTION_TYPE_LIST, SPACE> {
  typedef LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE> model_t;
  void get( const model_t* m, std::vector<T>& out ) const {
    typename model_t::FunctionIdentifier fid(0, FUNCTION_INDEX);
    for(size_t i=0; i<m->template functions<FUNCTION_INDEX>().size(); ++i){
      fid.functionIndex = i;
      if(m->feature_ids_.count(fid)) { // is weighted function an actual feature not just a function?
	out[ m->feature_ids_.find( fid )->second ] = m->template functions<FUNCTION_INDEX>()[i].weight();
      }
    }
    w_next_.get(m, out);
  };
  void set( model_t* m, std::vector<T>& in ) {
    typename model_t::FunctionIdentifier fid(0, FUNCTION_INDEX);
    for(size_t i=0; i<m->template functions<FUNCTION_INDEX>().size(); ++i){
      fid.functionIndex = i;
      if(m->feature_ids_.count(fid)) { // is weighted function an actual feature not just a function?
	m->template functions<FUNCTION_INDEX>()[i].weight( in[m->feature_ids_.find( fid )->second] );
      }
    }
    w_next_.set(m, in);
  };
  WeightAccessor<FUNCTION_INDEX-1, typename meta::TypeAtTypeList<FUNCTION_TYPE_LIST, FUNCTION_INDEX-1>::type, T, FUNCTION_TYPE_LIST, SPACE> w_next_;
};

template<class FUNCTION_TYPE, class T, class FUNCTION_TYPE_LIST, class SPACE>
struct WeightAccessor<0, FUNCTION_TYPE, T, FUNCTION_TYPE_LIST, SPACE> {
  typedef LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE> model_t;
  void get( const model_t* m, std::vector<T>& out ) const {
    assert(out.size() == m->numberOfWeights() );
  };
  void set( model_t* /*m*/, std::vector<T>& /*in*/ ) {
  };
};

template<class FUNCTION_TYPE, class T, class FUNCTION_TYPE_LIST, class SPACE>
struct WeightAccessor<0, opengm::FunctionDecoratorWeighted<FUNCTION_TYPE>, T, FUNCTION_TYPE_LIST, SPACE> {
  typedef LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE> model_t;
  void get( const model_t* m, std::vector<T>& out ) const {
    typename model_t::FunctionIdentifier fid(0, 0);
    for(size_t i=0; i<m->template functions<0>().size(); ++i){
      fid.functionIndex = i;
      if(m->feature_ids_.count(fid)) { // is weighted function an actual feature not just a function?
	out[ m->feature_ids_.find( fid )->second ] = m->template functions<0>()[i].weight();
      }
    }
  };
  void set( model_t* m, std::vector<T>& in ) {
    typename model_t::FunctionIdentifier fid(0, 0);
    for(size_t i=0; i<m->template functions<0>().size(); ++i){
      fid.functionIndex = i;
      if(m->feature_ids_.count(fid)) { // is weighted function an actual feature not just a function?
	m->template functions<0>()[i].weight( in[m->feature_ids_.find( fid )->second] );
      }
    }
  };
};



////
//// Accessing inner function values via TMP recursion
////
template <size_t FUNCTION_TYPE_INDEX, class FUNCTION_TYPE, class LLM>
struct InnerFunctionAccessor {
  template <class ITERATOR>
  static typename LLM::ValueType getFeatureValue
  (
   const LLM* m,
   ITERATOR iterator,
   const typename LLM::FunctionIdentifier::FunctionIndexType functionIndex,
   const typename LLM::FunctionIdentifier::FunctionTypeIndexType functionType
   ) {
    return InnerFunctionAccessor<FUNCTION_TYPE_INDEX -1,
    				   typename meta::TypeAtTypeList<typename LLM::FunctionTypeList, FUNCTION_TYPE_INDEX - 1>::type,
    				   LLM>::getFeatureValue(m, iterator, functionIndex, functionType);
  }
};
template <class FUNCTION_TYPE, class LLM>
struct InnerFunctionAccessor<0, FUNCTION_TYPE, LLM> {
  template <class ITERATOR>
  static typename LLM::ValueType getFeatureValue
  (
   const LLM* /*m*/,
   ITERATOR /*iterator*/,
   const typename LLM::FunctionIdentifier::FunctionIndexType /*functionIndex*/,
   const typename LLM::FunctionIdentifier::FunctionTypeIndexType /*functionType*/
   ) {
    throw std::runtime_error("InnerFunctionAccessor::getFeatureValue(): function not found");
  }
};
template <size_t FUNCTION_TYPE_INDEX, class FUNCTION_TYPE, class LLM>
struct InnerFunctionAccessor<FUNCTION_TYPE_INDEX, opengm::FunctionDecoratorWeighted<FUNCTION_TYPE>, LLM> {
  template <class ITERATOR>
  static typename LLM::ValueType getFeatureValue
  (
   const LLM* m,
   ITERATOR iterator,
   const typename LLM::FunctionIdentifier::FunctionIndexType functionIndex,
   const typename LLM::FunctionIdentifier::FunctionTypeIndexType functionType
   ) {
    if(functionType == FUNCTION_TYPE_INDEX) {
      return m->template functions<FUNCTION_TYPE_INDEX>()[functionIndex].innerFunction()->operator()(iterator);
    } else {
       return InnerFunctionAccessor<FUNCTION_TYPE_INDEX - 1,
				    typename meta::TypeAtTypeList<typename LLM::FunctionTypeList, FUNCTION_TYPE_INDEX - 1>::type,
				    LLM>::getFeatureValue(m, iterator, functionIndex, functionType);
    }
  }
};
template <class FUNCTION_TYPE, class LLM>
struct InnerFunctionAccessor<0, opengm::FunctionDecoratorWeighted<FUNCTION_TYPE>, LLM> {
  template <class ITERATOR>
  static typename LLM::ValueType getFeatureValue
  (
   const LLM* m,
   ITERATOR iterator,
   const typename LLM::FunctionIdentifier::FunctionIndexType functionIndex,
   const typename LLM::FunctionIdentifier::FunctionTypeIndexType functionType
   ) {
    if(functionType == 0) {
      return m->template functions<0>()[functionIndex].innerFunction()->operator()(iterator);
    } else {
      return 0;
    }
  }
};

} // namespace detail_graphical_model
/// \endcond


template<class T, class FUNCTION_TYPE_LIST, class SPACE>
template<class FEATURE_FUNCTION>
typename LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE>::FunctionIdentifier
LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE>::addFeature(FunctionDecoratorWeighted<FEATURE_FUNCTION> f, IndexType weightIndex) {
  if( !(weightIndex < this->numberOfWeights()) ) {
    throw std::invalid_argument("LoglinearModel::addFeatureWithRefReturn(): weightIndex out of range");
  }

  typedef typename LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE>::FunctionIdentifier FidType;
  f.weight( weights_[weightIndex] );
  FidType fid = this->addFunction(f);
  if(feature_ids_.empty()) {
    feature_ids_.insert( std::pair<FidType, size_t>(fid, weightIndex) );
  } else {
    feature_ids_.insert( --feature_ids_.end(), std::pair<FidType, size_t>(fid, weightIndex) );
  }
  return fid;
}

template<class T, class FUNCTION_TYPE_LIST, class SPACE>
void
LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE>::setWeights( std::vector<ValueType> in ) 
{
  if( in.size() != this->numberOfWeights() ) {
    throw std::invalid_argument("LoglinearModel::setWeights(): vector size not equal to number of weights");
  }
  weights_ = in;

  // update feature functions with new weight values
  typedef LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE> model_t;
  detail_graphical_model::WeightAccessor<model_t::NrOfFunctionTypes - 1,
					 typename meta::TypeAtTypeList<FUNCTION_TYPE_LIST, model_t::NrOfFunctionTypes -1>::type,
					 T, FUNCTION_TYPE_LIST, SPACE> w;
  w.set(this, in);
}

template<class T, class FUNCTION_TYPE_LIST, class SPACE>
void
LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE>::getWeights( std::vector<ValueType>& out ) const
{
  out.clear();
  out = weights_;

  // consistency check
  typedef LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE> model_t;
  detail_graphical_model::WeightAccessor<model_t::NrOfFunctionTypes - 1,
					 typename meta::TypeAtTypeList<FUNCTION_TYPE_LIST, model_t::NrOfFunctionTypes -1>::type,
					 T, FUNCTION_TYPE_LIST, SPACE> w;
  std::vector<ValueType> check( this->numberOfWeights() );
  w.get(this, check);
  assert(mismatch(out.begin(), out.end(), check.begin()).first == out.end());
}

template<class T, class FUNCTION_TYPE_LIST, class SPACE>
void 
LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE>::weightedFeatureSums( const std::vector<LabelType>& labelIndices, std::vector<ValueType>& out ) const
{
  assert(labelIndices.size() == this->numberOfVariables());
  out.resize(this->numberOfWeights());

  FunctionIdentifier fid;
  std::vector<size_t> factor_state;

  for(size_t j = 0; j < this->numberOfFactors(); ++j) {
    const FactorType& factor = (*this)[j];
    
    // check if factor is a feature
    fid.functionIndex = factor.functionIndex();
    fid.functionType = factor.functionType();
    if(feature_ids_.count(fid)) {
      // determine factor state from labelIndices
      factor_state.clear();
      factor_state.resize(factor.numberOfVariables());
      for(size_t i = 0; i < factor.numberOfVariables(); ++i) {
         assert(labelIndices[factor.variableIndex(i)] < factor.numberOfLabels(i));
         factor_state[i] = labelIndices[factor.variableIndex(i)];
      }
      // retrieve factor value at state
      typedef LoglinearModel<T, FUNCTION_TYPE_LIST, SPACE> model_t;
      detail_graphical_model::InnerFunctionAccessor<model_t::NrOfFunctionTypes - 1,
       						    typename meta::TypeAtTypeList<FunctionTypeList, model_t::NrOfFunctionTypes -1>::type,
       						    model_t> a;
      ValueType v = a.getFeatureValue(this, factor_state.begin(), factor.functionIndex(), factor.functionType());
      out[feature_ids_.find(fid)->second] += v;
    }
  }
}

} // namespace opengm
#endif // #ifndef OPENGM_LOGLINEARMODEL_HXX

