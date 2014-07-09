/**
   @file
   @ingroup tracking
   @brief tracking events like Move and Division
*/

#ifndef EVENT_H
#define EVENT_H
#include <vector>
#include <ostream>
#include <stdexcept>
#include <stdint.h>

#include "pgmlink/pgmlink_export.h"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

// template class PGMLINK_EXPORT std::vector<unsigned int>;

namespace pgmlink {

class Event 
{
 public:
    PGMLINK_EXPORT Event() 
    : n_features_(0) 
    { 
      type = Void;
    }
    
    enum EventType {Move, Division, Appearance, Disappearance, Merger, ResolvedTo, MultiFrameMove, Void};
    EventType type;
    
    /**
       \brief The traxels participating in the event.
       
       A (Dis-)appearance is described by one traxel id, a move by two
       traxel ids and a division by three (motherId, daughterId1,
       daughterId2). A Void is just an empty vector. ResolvedTo contains the
       id of the merger node followed by pairs of <new id, respective COM>.
    */
    std::vector<std::size_t> traxel_ids;
    
    /**
       \brief The energy assigned to this Event happening.
       
       It's calculated as the scalar product of weights and features.
       If there are no features the energy will be 0 by definition.
    */
    PGMLINK_EXPORT double energy() const;
    
    PGMLINK_EXPORT unsigned int number_of_features() const { return n_features_; }
    PGMLINK_EXPORT Event& number_of_features( unsigned int );

    PGMLINK_EXPORT const std::vector<double>& features() const { return features_; }
    PGMLINK_EXPORT Event& features( const std::vector<double>& );
    
    PGMLINK_EXPORT const std::vector<double>& weights() const { return weights_; }
    PGMLINK_EXPORT Event& weights( const std::vector<double>& );
    
    // Beware: doesn't compare by energy!
    PGMLINK_EXPORT bool operator==(const Event& other) const;
    PGMLINK_EXPORT bool operator!=(const Event& other) const;
    PGMLINK_EXPORT bool operator<(const Event& other) const;
    PGMLINK_EXPORT bool operator>(const Event& other) const;
    friend std::ostream& operator<< (std::ostream &out, const Event &e);
    
  private:
    uint32_t n_features_;
    std::vector<double> weights_;
    std::vector<double> features_;

  private:
    // serialization interface
    friend class boost::serialization::access;

    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/)
    {
        ar & type;
        ar & traxel_ids;
        ar & n_features_;
        ar & weights_;
        ar & features_;
    }
  };

struct EventsStatistics 
{
    PGMLINK_EXPORT EventsStatistics() 
    : n_total(0), n_mov(0), n_div(0), n_app(0), n_dis(0) 
    {}

    PGMLINK_EXPORT EventsStatistics& operator+=( const EventsStatistics& rhs ) 
    {
      n_total += rhs.n_total;
      n_mov += rhs.n_mov;
      n_div += rhs.n_div;
      n_app += rhs.n_app;
      n_dis += rhs.n_dis;
      return *this;
    }
    size_t n_total, n_mov, n_div, n_app, n_dis;
};

inline EventsStatistics operator+( EventsStatistics lhs, const EventsStatistics& rhs ) {
  lhs += rhs;
  return lhs;
}
inline std::ostream& operator<<(std::ostream& os, const EventsStatistics& obj) { 
  os << "#total: " << obj.n_total
     << ", #move: " << obj.n_mov
     << ", #division: " << obj.n_div
     << ", #app.: " << obj.n_app
     << ", #disapp.: " << obj.n_dis;
  return os;
} 

template<class event_it>
  EventsStatistics collect_events_statistics(event_it begin, event_it end, bool ignore_unknown_events=false) {
  EventsStatistics stats;
  for(event_it event=begin; event!=end; ++event) {
    switch(event->type) {
    case Event::Move:
  stats.n_mov += 1;
  stats.n_total += 1;
  break;
    case Event::Division:
  stats.n_div += 1;
  stats.n_total += 1;
  break;
    case Event::Appearance:
  stats.n_app += 1;
  stats.n_total += 1;
  break;
    case Event::Disappearance:
  stats.n_dis += 1;
  stats.n_total += 1;
  break;
    default:
  if(!ignore_unknown_events) {
    throw std::runtime_error("collect_events_statistics(): unknown event type encountered");
  }
  break;
    }
  }
  return stats;
}

} /* namespace pgmlink */

#endif /* EVENT_H */
