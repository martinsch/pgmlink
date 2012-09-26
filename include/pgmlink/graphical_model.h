#ifndef GRAPHICAL_MODEL_H
#define GRAPHICAL_MODEL_H

#include <vector>
#include <utility>
#include <string>
#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <opengm/operations/adder.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/inference/lpcplex.hxx>
#include <vigra/multi_array.hxx>

#include "pgmlink/hypotheses.h"
#include "pgmlink/graph.h"

using boost::vecS;
using boost::bidirectionalS;
using boost::shared_ptr;

namespace Tracking {

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



template <typename VALUE_T>
  class OpengmBinaryFactor {
 public:
  OpengmBinaryFactor( const std::vector<size_t>& ogm_var_indices, VALUE_T init=0 );

  OpengmBinaryFactor& set_value( const std::vector<size_t> coords, VALUE_T v);
  VALUE_T get_value( const std::vector<size_t> coords ) const;

  const OpengmModel::ExplicitFunctionType& OgmFunction() const;

 private:
  OpengmModel::ExplicitFunctionType ogmfunction_;
  std::vector<size_t> vi_;
};    



/******************/
/* Implementation */
/******************/

////
//// class OpengmBinaryFactor
////
 template <typename VALUE_T>
   OpengmBinaryFactor<VALUE_T>::OpengmBinaryFactor( const std::vector<size_t>& ogm_var_indices, VALUE_T init ) : vi_(ogm_var_indices) {
  std::vector<size_t> shape( vi_.size(), 2);
  ogmfunction_ = OpengmModel::ExplicitFunctionType( shape.begin(), shape.end(), init );
}

 template <typename VALUE_T>  
   OpengmBinaryFactor<VALUE_T>& OpengmBinaryFactor<VALUE_T>::set_value( const std::vector<size_t> coords, VALUE_T v) {
   if( coords.size() != vi_.size() ) {
     throw std::invalid_argument("OpengmBinaryFactor::set_value(): coordinate dimension differs from factor dimension");
   }
   OpengmModel::ExplicitFunctionType::iterator element( ogmfunction_ );
   size_t index;
   ogmfunction_.coordinatesToIndex(coords.begin(), index);
   element[index] = v;
   return *this;
 }

 template <typename VALUE_T>
   VALUE_T OpengmBinaryFactor<VALUE_T>::get_value( const std::vector<size_t> coords ) const {
   if( coords.size() != vi_.size() ) {
     throw std::invalid_argument("OpengmBinaryFactor::get_value(): coordinate dimension differs from factor dimension");
   }
   OpengmModel::ExplicitFunctionType::iterator element( ogmfunction_ );
   size_t index;
   ogmfunction_.coordinatesToIndex(coords.begin(), index);
   return element[index];
 }

 template <typename VALUE_T>
   const OpengmModel::ExplicitFunctionType& OpengmBinaryFactor<VALUE_T>::OgmFunction() const {
   return ogmfunction_;
 }
 

} /* namespace Tracking */

#endif /* GRAPHICAL_MODEL_H */
