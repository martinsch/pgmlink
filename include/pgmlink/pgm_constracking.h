#ifndef PGM_CONSTRACKING_H
#define PGM_CONSTRACKING_H

#include "pgm.h"
#include "event.h"

namespace pgm {

template <typename VALUE>
   class OneEventExplicitFunction: public OpengmEventExplicitFunction<VALUE> {

       template <class SHAPE_ITERATOR>
       OneEventExplicitFunction(SHAPE_ITERATOR shapeBegin, SHAPE_ITERATOR shapeEnd, const VALUE & value,
                                   const FeatureVector& features, const WeightVector& weights,
                                   const EventMap& event_map):
           OpengmEventExplicitFunction(shapeBegin, shapeEnd, value, features, weights, event_map) {
           if(event_map.size() != 1) {
               throw std::exception("only one event is allowed in the one event function");
           }
       }

       virtual VALUE get_energy_of_configuration(const std::vector<size_t>& configuration) {
           return get_event_energy(event_map_.begin()->first);
       }
   };


   //TODO: Specific event functions for outgoing, etc.


} // namespace pgm
#endif // PGM_CONSTRACKING_H
