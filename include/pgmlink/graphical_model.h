#ifndef GRAPHICAL_MODEL_H
#define GRAPHICAL_MODEL_H

#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <opengm/functions/explicit_function.hxx>
#include <opengm/operations/adder.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/inference/lpcplex.hxx>
#include <vigra/multi_array.hxx>

#include "pgmlink/hypotheses.h"
#include "pgmlink/graph.h"
#include "pgmlink/util.h"

using boost::vecS;
using boost::bidirectionalS;
using boost::shared_ptr;

namespace pgmlink {
template <typename VALUE>
  class OpengmBinaryFactor;

class OpengmModel {
   public:
   typedef double Energy;
   typedef opengm::GraphicalModel<Energy, opengm::Adder> ogmGraphicalModel;
   typedef opengm::Factor<ogmGraphicalModel> ogmFactor;
   typedef opengm::Minimizer ogmAccumulator;
   typedef opengm::Inference<ogmGraphicalModel, ogmAccumulator> ogmInference;
   typedef opengm::meta::TypeAtTypeList<ogmGraphicalModel::FunctionTypeList, 0>::type ExplicitFunctionType;
   typedef ogmGraphicalModel::FunctionIdentifier FunctionIdentifier;
   
   OpengmModel();
   ~OpengmModel();

   ogmGraphicalModel* Model() {return model_; }
   const ogmGraphicalModel* Model() const {return model_; }
 
   private:
   OpengmModel(const OpengmModel&);
   OpengmModel& operator=(const OpengmModel&);
   
   ogmGraphicalModel* model_;
};



template <typename OGM_FUNCTION>
 class OpengmFactor {
 public:
  typedef OGM_FUNCTION FunctionType;

  OpengmFactor( const FunctionType&, const std::vector<size_t>& ogm_var_indices );
  template <typename ITER>
    OpengmFactor( const FunctionType&, ITER first_ogm_idx, ITER last_ogm_idx );

  typename FunctionType::ValueType get_value( std::vector<size_t> coords ) const;
  void add_to( OpengmModel& ) const;

  const FunctionType& function() const { return ogmfunction_; }
  const std::vector<size_t>& var_indices() const { return ogmfunction_; }

 protected:
  FunctionType ogmfunction_;
  std::vector<size_t> vi_;
  std::vector<size_t> order_;
 };


template <typename VALUE>
  class OpengmExplicitFactor : public OpengmFactor<opengm::ExplicitFunction<VALUE> > {
 public:
  OpengmExplicitFactor( const std::vector<size_t>& ogm_var_indices, VALUE init=0, size_t states_per_var=2 );
  template <typename ITER>
    OpengmExplicitFactor( ITER first_ogm_idx, ITER last_ogm_idx, VALUE init=0, size_t states_per_var=2 );  

  void set_value( std::vector<size_t> coords, VALUE v);

 private:
  void init_( VALUE init, size_t states_per_var );
};    



/******************/
/* Implementation */
/******************/

////
//// class OpengmFactor
////
 template <typename OGM_FUNCTION>
   OpengmFactor<OGM_FUNCTION>::OpengmFactor( const FunctionType& f, const std::vector<size_t>& ogm_var_indices )
   : ogmfunction_(f), vi_(ogm_var_indices) {
  // store sorted order of indices
  indexsorter::sort_indices( vi_.begin(), vi_.end(), order_ );
}

 template <typename OGM_FUNCTION>
   template< typename ITER >
   OpengmFactor<OGM_FUNCTION>::OpengmFactor( const FunctionType& f, ITER first_ogm_idx, ITER last_ogm_idx )
   : ogmfunction_(f), vi_(first_ogm_idx, last_ogm_idx) {
  // store sorted order of indices
  indexsorter::sort_indices( vi_.begin(), vi_.end(), order_ );
 }

 template <typename OGM_FUNCTION>
   typename OpengmFactor<OGM_FUNCTION>::FunctionType::ValueType OpengmFactor<OGM_FUNCTION>::get_value( std::vector<size_t> coords ) const {
   if( coords.size() != vi_.size() ) {
     throw std::invalid_argument("OpengmFactor::get_value(): coordinate dimension differs from factor dimension");
   }
   indexsorter::reorder( coords, order_ );
   return ogmfunction_(coords.begin());
 }

 template <typename OGM_FUNCTION>
   void OpengmFactor<OGM_FUNCTION>::add_to( OpengmModel& m ) const {
   std::vector<size_t> sorted_vi(vi_);
   std::sort(sorted_vi.begin(), sorted_vi.end());
   // opengm expects a monotonic increasing sequence
   if(!(m.Model()->isValidIndexSequence(sorted_vi.begin(), sorted_vi.end()))) {
      throw std::runtime_error("OpengmExplicitFactor::add_to(): invalid index sequence");
   }

   OpengmModel::FunctionIdentifier id=m.Model()->addFunction( ogmfunction_ );
   m.Model()->addFactor(id, sorted_vi.begin(), sorted_vi.end());
 }


 
////
//// class OpengmExplicitFactor
////
 template <typename VALUE>
   OpengmExplicitFactor<VALUE>::OpengmExplicitFactor( const std::vector<size_t>& ogm_var_indices, VALUE init, size_t states_per_var )
   : OpengmFactor<opengm::ExplicitFunction<VALUE> >(opengm::ExplicitFunction<VALUE>(), ogm_var_indices) {
  init_( init, states_per_var );
}
 template <typename VALUE>
   template< typename ITER >
   OpengmExplicitFactor<VALUE>::OpengmExplicitFactor( ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var )
   : OpengmFactor<opengm::ExplicitFunction<VALUE> >(opengm::ExplicitFunction<VALUE>(),first_ogm_idx, last_ogm_idx) {
  init_( init, states_per_var );
 }

 template <typename VALUE>  
   void OpengmExplicitFactor<VALUE>::set_value( std::vector<size_t> coords, VALUE v) {
   if( coords.size() != this->vi_.size() ) {
     throw std::invalid_argument("OpengmExplicitFactor::set_value(): coordinate dimension differs from factor dimension");
   }
   typename opengm::ExplicitFunction<VALUE>::iterator element( this->ogmfunction_ );
   size_t index;

   indexsorter::reorder( coords, this->order_ );
   this->ogmfunction_.coordinatesToIndex(coords.begin(), index);
   element[index] = v;
 }

 template <typename VALUE>
   void OpengmExplicitFactor<VALUE>::init_( VALUE init, size_t states_per_var ) {
   std::vector<size_t> shape( this->vi_.size(), states_per_var );
   this->ogmfunction_ = OpengmModel::ExplicitFunctionType( shape.begin(), shape.end(), init );
 }
} /* namespace pgmlink */

#endif /* GRAPHICAL_MODEL_H */
