/**
   @file
   @ingroup pgm
   @brief opengm extensions
*/

#ifndef GRAPHICAL_MODEL_H
#define GRAPHICAL_MODEL_H

#include <vector>

#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/graphicalmodel/loglinearmodel.hxx>
#include <opengm/inference/inference.hxx>
#include <opengm/functions/explicit_function.hxx>
#include <opengm/functions/decorator_weighted.hxx>
#include <opengm/functions/indicator_function.hxx>
#include <opengm/operations/adder.hxx>
#include <opengm/utilities/metaprogramming.hxx>

#include "pgmlink/hypotheses.h"
#include "pgmlink/graph.h"
#include "pgmlink/util.h"

namespace pgmlink {
namespace pgm {

//typedef opengm::GraphicalModel<double, opengm::Adder> OpengmModel;
typedef opengm::ExplicitFunction<double> ExplicitFunction;
typedef opengm::FunctionDecoratorWeighted< opengm::IndicatorFunction<double> > FeatureFunction;
typedef opengm::LoglinearModel<double,
  opengm::meta::TypeList<ExplicitFunction,
  opengm::meta::TypeList<FeatureFunction, opengm::meta::ListEnd > > > OpengmModel;

template <typename OGM_FUNCTION>
 class OpengmFactor {
 public:
  typedef OGM_FUNCTION FunctionType;

  OpengmFactor( const FunctionType&, const std::vector<size_t>& ogm_var_indices );
  template <typename ITER>
    OpengmFactor( const FunctionType&, ITER first_ogm_idx, ITER last_ogm_idx );

  typename FunctionType::ValueType get_value( std::vector<size_t> coords ) const;
  template <typename OGM_GRAPHICAL_MODEL>
    void add_to( OGM_GRAPHICAL_MODEL& ) const;

  FunctionType& function() { return ogmfunction_; }
  const FunctionType& function() const { return ogmfunction_; }
  const std::vector<size_t>& var_indices() const { return vi_; }

 protected:
  FunctionType ogmfunction_;
  std::vector<size_t> vi_;
  std::vector<size_t> order_;
 };


template <typename VALUE>
  class OpengmExplicitFactor : public OpengmFactor<opengm::ExplicitFunction<VALUE> > {
 public:
  typedef typename OpengmFactor<opengm::ExplicitFunction<VALUE> >::FunctionType FunctionType; 

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
   template <typename OGM_GRAPHICAL_MODEL>
   void OpengmFactor<OGM_FUNCTION>::add_to( OGM_GRAPHICAL_MODEL& m ) const {
   std::vector<size_t> sorted_vi(vi_);
   std::sort(sorted_vi.begin(), sorted_vi.end());
   // opengm expects a monotonic increasing sequence
   if(!(m.isValidIndexSequence(sorted_vi.begin(), sorted_vi.end()))) {
     throw std::runtime_error("OpengmExplicitFactor::add_to(): invalid index sequence");
   }
   
   typename OGM_GRAPHICAL_MODEL::FunctionIdentifier id=m.addFunction( ogmfunction_ );
   m.addFactor(id, sorted_vi.begin(), sorted_vi.end());
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
   this->ogmfunction_ = opengm::ExplicitFunction<VALUE>( shape.begin(), shape.end(), init );
 }

} /* namespace pgm */
} /* namespace pgmlink */

#endif /* GRAPHICAL_MODEL_H */
