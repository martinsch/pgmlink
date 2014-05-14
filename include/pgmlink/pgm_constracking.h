#ifndef PGM_CONSTRACKING_H
#define PGM_CONSTRACKING_H

#include "pgm.h"
#include "event.h"

namespace pgmlink{
namespace pgm {

// ------------------------------------------------------
template <typename VALUE>
class UnaryEventExplicitFunction: public OpengmEventExplicitFunction<VALUE> {

    template <class SHAPE_ITERATOR>
    OneEventExplicitFunction(SHAPE_ITERATOR shapeBegin, SHAPE_ITERATOR shapeEnd, const VALUE & value,
                             const FeatureMap& features, const WeightVector& weights,
                             const EventMap& event_map):
        OpengmEventExplicitFunction(shapeBegin, shapeEnd, value, features, weights, event_map)
    {
        if(event_map.size() != 1) {
            throw std::exception("only one event is allowed in the one event function");
        }
    }

    virtual VALUE get_energy_of_configuration(const std::vector<size_t>& configuration)
    {
        assert( configuration.size() == 1 && "only unaries are allowed" );

        return get_event_energy(event_map_.begin()->first, configuration[0]);
    }
};


// ------------------------------------------------------
template <typename VALUE>
class PairwiseDetectionExplicitFunction: public OpengmEventExplicitFunction<VALUE>
{
public:
    typedef std::vector< std::pair<Event::EventType, size_t> > EventConfigurationVector;

    template <class SHAPE_ITERATOR>
    PairwiseDetectionExplicitFunction(SHAPE_ITERATOR shapeBegin,
                                      SHAPE_ITERATOR shapeEnd,
                                      const VALUE & value,
                                      const FeatureMap& features,
                                      const WeightVector& weights,
                                      const EventMap& event_map):
        OpengmEventExplicitFunction(shapeBegin, shapeEnd, value, features, weights, event_map)
    {}

    virtual VALUE get_energy_of_configuration(const std::vector<size_t>& configuration)
    {
        assert( configuration.size() == 2 && "only pairwise interactions are allowed");

        VALUE energy = 0;
        EventConfigurationVector event_configurations = determine_event_configuration(configuration);

        for(EventConfigurationVector::iterator event_config_it = event_configurations.begin();
            event_config_it != event_configurations.end();
            ++event_config_it)
        {
            energy += get_event_energy(event_config_it->first, event_config_it->second);
        }

        return energy;
    }

    EventConfigurationVector determine_event_configurations(const std::vector<size_t>& configuration)
    {
        EventConfigurationVector event_configurations;

        if(configuration[0] == configuration[1])
        {
            event_configurations.push_back(std::make_pair(Event::Merger, configuration[0]));
        }
        else if(configuration[0] < configuration[1] && configuration[0] == 0)
        {
            event_configurations.push_back(std::make_pair(Event::Appearance, configuration[1]));
        }
        else if(configuration[0] > configuration[1] && configuration[1] == 0)
        {
            event_configurations.push_back(std::make_pair(Event::Disappearance, configuration[0]));
        }

        return event_configurations;
    }
};
//TODO: Specific event functions for outgoing, etc.



} // namespace pgm
} // namespace pgmlink

#endif // PGM_CONSTRACKING_H
