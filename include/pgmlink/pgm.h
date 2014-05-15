/**
   @file
   @ingroup pgm
   @brief opengm extensions
*/

#ifndef PGMLINK_PGM_H
#define PGMLINK_PGM_H

#include <vector>

#include <pgmlink/ext_opengm/loss_hamming.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <pgmlink/ext_opengm/loglinearmodel.hxx>
#include <opengm/inference/inference.hxx>
#include <opengm/functions/explicit_function.hxx>
#include <opengm/functions/function_properties_base.hxx>
#include <pgmlink/ext_opengm/decorator_weighted.hxx>
#include <pgmlink/ext_opengm/indicator_function.hxx>
#include <opengm/operations/adder.hxx>
#include <opengm/utilities/metaprogramming.hxx>

#include "pgmlink/hypotheses.h"
#include "pgmlink/graph.h"
#include "pgmlink/util.h"
#include "event.h"

namespace pgmlink {
  namespace pgm {
    using boost::shared_ptr;
    


  template <typename VALUE>
     class OpengmEventExplicitFunction:
             public opengm::FunctionBase<OpengmEventExplicitFunction<VALUE>,
                                         VALUE, /*INDEX_TYPE*/ size_t, /*LABEL_TYPE*/ size_t> {

     public:
         typedef VALUE ValueType;
         typedef size_t LabelType;
         typedef size_t IndexType;

         typedef typename opengm::FunctionBase<OpengmEventExplicitFunction<VALUE>, VALUE, /*INDEX_TYPE*/ size_t, /*LABEL_TYPE*/ size_t> FunctionBaseType;

         typedef std::pair<std::string, std::string> FeatureName;
         typedef std::map<FeatureName,pgmlink::feature_array> FeatureOperatorMap;
         typedef typename std::vector<VALUE> WeightVector;
         typedef std::pair<size_t,FeatureName > WeightFeaturePair;
         typedef std::map<size_t,std::vector<WeightFeaturePair > > EventConfigurationMap;
         typedef std::map<pgmlink::Event::EventType,EventConfigurationMap > EventMap;

         template <class SHAPE_ITERATOR>
         OpengmEventExplicitFunction(SHAPE_ITERATOR shapeBegin, SHAPE_ITERATOR shapeEnd, const VALUE & value,
                                     const FeatureOperatorMap& features, WeightVector* weights,
                                     const EventMap& event_map) :
             FunctionBaseType(),
             feature_map_(features), weight_vector_(weights), event_map_(event_map), value_(value), shape_(shapeBegin,shapeEnd)
             {   }

         template <typename ITERATOR>
         VALUE operator()(ITERATOR it) const {
             std::vector<size_t> configuration;
             for(size_t i = 0; i < this->dimension(); ++i) {
                 configuration.push_back(it[i]);
             }
             return get_energy_of_configuration(configuration);
         }


         size_t shape(const size_t var_idx) const {
            return shape_[var_idx];
         }

         size_t size() const {
             size_t res = 1;
             for(std::vector<std::size_t>::const_iterator it = shape_.begin(); it != shape_.end(); ++it ) {
                 res *= (*it);
             }
             return res;
         }

         size_t dimension() const {
             return shape_.size();
         }

         size_t numberOfParameters() const {
             std::set<size_t> weight_indexes;
             for (EventMap::const_iterator ev_type_it = event_map_.begin(); ev_type_it != event_map_.end(); ++ev_type_it) {
                 for (EventConfigurationMap::const_iterator ev_config_it = ev_type_it->second.begin();
                      ev_config_it != ev_type_it->second.end(); ++ev_config_it) {
                     for (std::vector<WeightFeaturePair>::const_iterator wtf_it = ev_config_it->second.begin();
                          wtf_it != ev_config_it->second.end(); ++wtf_it) {
                         weight_indexes.insert(wtf_it->first);
                     }
                 }
             }
             return weight_indexes.size();
         }

         IndexType parameterIndex(const size_t /*paramNumber*/) const {
             // returns the global weight index for the local paramNumber
             throw std::runtime_error("not implemented");
         }

     protected:
         virtual VALUE get_energy_of_configuration(const std::vector<size_t>&) const {
             throw std::runtime_error("not implemented");
         }

           // returns inner product of <w,f(x)> for given event with its specific configuration
         VALUE get_event_energy(Event::EventType event_name, size_t event_configuration) const {
            VALUE energy = 0;
            EventMap::const_iterator event_it = event_map_.find(event_name);
            if (event_it == event_map_.end()) {
                throw std::runtime_error("event not found");
            }
            const EventConfigurationMap& ev_config_map = event_it->second;

            EventConfigurationMap::const_iterator ev_config_it = ev_config_map.find(event_configuration);
            if (ev_config_it == ev_config_map.end()) {
                throw std::runtime_error("event configuration not found");
            }

            const std::vector<WeightFeaturePair>& weight_feature_index_pairs = ev_config_it->second;

            for(std::vector<WeightFeaturePair>::const_iterator it = weight_feature_index_pairs.begin();
                it != weight_feature_index_pairs.end(); ++it) {
                const pgmlink::feature_array& feats = feature_map_.at(it->second); // access by at() for const
                for(pgmlink::feature_array::const_iterator feat_it = feats.begin(); feat_it != feats.end(); ++feat_it) {
                    energy += ((*weight_vector_)[it->first]) * (*feat_it);
                }
            }

            return energy;
         }


       FeatureOperatorMap feature_map_;
       WeightVector* weight_vector_;
       EventMap event_map_;
       VALUE value_;
       std::vector<std::size_t> shape_;
     };




    //typedef opengm::GraphicalModel<double, opengm::Adder> OpengmModel;
    typedef opengm::ExplicitFunction<double> ExplicitFunction;
    typedef opengm::FunctionDecoratorWeighted< opengm::IndicatorFunction<double> > FeatureFunction;
    typedef opengm::HammingFunction<double> LossFunction;
    typedef opengm::LoglinearModel<double,
      opengm::meta::TypeList<ExplicitFunction,
      opengm::meta::TypeList<FeatureFunction,
      opengm::meta::TypeList<LossFunction, opengm::meta::ListEnd > > > > OpengmModel;

    class OpengmModelDeprecated {
    public:
      typedef double Energy;
      typedef OPENGM_TYPELIST_2(ExplicitFunction, pgmlink::pgm::OpengmEventExplicitFunction<double>) ogmFunctionsTypelist;
      typedef opengm::GraphicalModel<Energy, opengm::Adder, ogmFunctionsTypelist> ogmGraphicalModel;
      typedef opengm::Factor<ogmGraphicalModel> ogmFactor;
      typedef opengm::Minimizer ogmAccumulator;
      typedef opengm::Inference<ogmGraphicalModel, ogmAccumulator> ogmInference;
      typedef opengm::meta::TypeAtTypeList<ogmGraphicalModel::FunctionTypeList, 0>::type ExplicitFunctionType;
      typedef ogmGraphicalModel::FunctionIdentifier FunctionIdentifier;
      
      OpengmModelDeprecated() {
        model_ = new ogmGraphicalModel();
      }
      ~OpengmModelDeprecated() {
        delete model_;
      }
      
      ogmGraphicalModel* Model() {return model_; }
      const ogmGraphicalModel* Model() const {return model_; }
      
      std::vector<size_t> get_number_of_labels(const std::vector<size_t>& var_idxes) {
          std::vector<size_t> result;
          for(std::vector<size_t>::const_iterator it = var_idxes.begin(); it != var_idxes.end(); ++it) {
              result.push_back(model_->numberOfLabels(*it));
          }
          return result;
      }

    private:
      OpengmModelDeprecated(const OpengmModelDeprecated&);
      OpengmModelDeprecated& operator=(const OpengmModelDeprecated&);
            
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
      template <typename OGM_GRAPHICAL_MODEL>
	void add_to( OGM_GRAPHICAL_MODEL& ) const;

      FunctionType& function() { return ogmfunction_; }
      const FunctionType& function() const { return ogmfunction_; }
      const std::vector<size_t>& var_indices() const { return vi_; }
      const std::vector<size_t>& var_order() const { return order_; }

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
      template <typename ITER>
      OpengmExplicitFactor( ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, std::vector<size_t> states_vars );
  
      void set_value( std::vector<size_t> coords, VALUE v);
  
    private:
      void init_( VALUE init, size_t states_per_var );
      void init_( VALUE init, std::vector<size_t> states_vars );
    };    


    template <typename VALUE>
      class OpengmWeightedFeature
      : public OpengmFactor< opengm::FunctionDecoratorWeighted< opengm::IndicatorFunction<VALUE> > >
    {
    public:
      typedef typename opengm::FunctionDecoratorWeighted< opengm::IndicatorFunction<VALUE> > FunctionType;

      template<class IT1, class IT2>
	OpengmWeightedFeature(const std::vector<size_t>& ogm_var_indices,
			      IT1 shapeBegin, IT1 shapeEnd,
			      IT2 indicate,
			      VALUE indicate_value = 1., 
			      VALUE weight = 1.,
			      VALUE default_value = 0);

      template <typename OGM_LOGLINEARMODEL>
	void add_as_feature_to( OGM_LOGLINEARMODEL&, typename OGM_LOGLINEARMODEL::IndexType weight_index ) const;

      void weight( VALUE w );
      VALUE weight() const;

      void indicate_value( VALUE v );
      VALUE indicate_value() const;
    
      void default_value( VALUE v );
      VALUE default_value() const;
    };


    /**
       @brief Accessing entries of a Factor/Function that was already added to a graphical model.

       Manages a pointer to an element of an array-like opengm function (usually opengm::ExplicitFunction).
       Validity of the pointer is ensured by owning a smart pointer to the full model.

       Use this class to modify factor elements of an already instantiated opengm graphical model.
    */
    class FactorEntry {
    public:
    FactorEntry() : entry_(NULL) {}
    FactorEntry( shared_ptr<OpengmModel> m, /**< has to be valid */
		 OpengmModel::ValueType* entry /**< has to point into the opengm model to ensure the same lifetime */
		 ) :
      m_(m), entry_(entry) {}
      
      void set( OpengmModel::ValueType );
      OpengmModel::ValueType get() const;

      shared_ptr<OpengmModel> model() const { return m_; }

    private:
      shared_ptr<OpengmModel> m_;
      OpengmModel::ValueType* entry_;
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
   indexsorter::reorder( sorted_vi, order_);
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
    template< typename ITER >
    OpengmExplicitFactor<VALUE>::OpengmExplicitFactor( ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, std::vector<size_t> states_vars )
    : OpengmFactor<opengm::ExplicitFunction<VALUE> >(opengm::ExplicitFunction<VALUE>(),first_ogm_idx, last_ogm_idx) {
   init_( init, states_vars );
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

 template <typename VALUE>
    void OpengmExplicitFactor<VALUE>::init_( VALUE init, std::vector<size_t> states_vars ) {
    this->ogmfunction_ = opengm::ExplicitFunction<VALUE>( states_vars.begin(), states_vars.end(), init );
  }



////
//// class OpengmWeightedFeature
////
template<class VALUE>
  template<class IT1, class IT2>
  OpengmWeightedFeature<VALUE>::OpengmWeightedFeature(
			  const std::vector<size_t>& ogm_var_indices,
		          IT1 shapeBegin, IT1 shapeEnd,
			  IT2 indicate,
			  VALUE indicate_value, 
			  VALUE weight,
			  VALUE default_value)
  : OpengmFactor<FunctionType >(FunctionType(new opengm::IndicatorFunction<VALUE>( shapeBegin,
										   shapeEnd,
										   indicate,
										   indicate_value,
										   default_value),
					     weight)
				, ogm_var_indices) {
   // replace function with a sorted variant
   std::vector<size_t> sorted_shape( shapeBegin, shapeEnd );
   indexsorter::reorder( sorted_shape, this->order_ );
   std::vector<size_t> sorted_indicate( this->ogmfunction_.innerFunction()->indicate() );
   indexsorter::reorder( sorted_indicate, this->order_ );

   this->ogmfunction_ = FunctionType(new opengm::IndicatorFunction<VALUE>( sorted_shape.begin(),
									   sorted_shape.end(),
									   sorted_indicate.begin(),
									   indicate_value,
									   default_value),
				     weight);
 }

template<class VALUE>
  template <typename OGM_LOGLINEARMODEL>
  void OpengmWeightedFeature<VALUE>::add_as_feature_to( OGM_LOGLINEARMODEL& m, typename OGM_LOGLINEARMODEL::IndexType weight_index ) const {
  std::vector<size_t> sorted_vi(this->vi_);
  std::sort(sorted_vi.begin(), sorted_vi.end());
  // opengm expects a monotonic increasing sequence
  if(!(m.isValidIndexSequence(sorted_vi.begin(), sorted_vi.end()))) {
    throw std::runtime_error("OpengmExplicitFactor::add_to(): invalid index sequence");
  }
  
  typename OGM_LOGLINEARMODEL::FunctionIdentifier id=m.addFeature( this->ogmfunction_, weight_index );
  m.addFactor(id, sorted_vi.begin(), sorted_vi.end());
 }

template<class VALUE>
  void OpengmWeightedFeature<VALUE>::weight( VALUE w ) {
  this->ogmfunction_.weight(w);
 }

template<class VALUE>
  VALUE OpengmWeightedFeature<VALUE>::weight() const {
  return this->ogmfunction_.weight();
 }

template<class VALUE>
  void OpengmWeightedFeature<VALUE>::indicate_value( VALUE v ) {
  this->ogmfunction_.innerFunction()->indicate_value(v);
 }
  
template<class VALUE>
  VALUE OpengmWeightedFeature<VALUE>::indicate_value() const {
  return this->ogmfunction_.innerFunction()->indicate_value();
 }
  
template<class VALUE>    
  void OpengmWeightedFeature<VALUE>::default_value( VALUE v ) {
  this->ogmfunction_.innerFunction()->default_value(v);
 }

template<class VALUE>
  VALUE OpengmWeightedFeature<VALUE>::default_value() const {
  return this->ogmfunction_.innerFunction()->default_value();
 }


} /* namespace pgm */
} /* namespace pgmlink */

#endif /* PGMLINK_PGM_H */
