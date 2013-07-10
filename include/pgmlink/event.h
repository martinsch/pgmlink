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

#include "pgmlink/pgmlink_export.h"



template class PGMLINK_EXPORT std::vector<unsigned int>;

namespace pgmlink {
  class PGMLINK_EXPORT Event {
  public:
  Event() : n_features_(0) { type = Void; };
    
    enum EventType {Move, Division, Appearance, Disappearance, Void};
    EventType type;
    
    /**
       \brief The traxels participating in the event.
       
       A (Dis-)appearance is described by one traxel id, a move by two
       traxel ids and a division by three (motherId, daughterId1,
       daughterId2). A Void is just an empty vector.
     */
    std::vector<unsigned int> traxel_ids;

    /**
       \brief The energy assigned to this Event happening.
     
       It's calculated as the scalar product of weights and features.
       If there are no features the energy will be 0 by definition.
     */
    double energy() const;

    unsigned int number_of_features() const { return n_features_; }
    Event& number_of_features( unsigned int );

    const std::vector<double>& features() const { return features_; }
    Event& features( const std::vector<double>& );

    const std::vector<double>& weights() const { return weights_; }
    Event& weights( const std::vector<double>& );

    // Beware: doesn't compare by energy!
    bool operator==(const Event& other) const;
    bool operator!=(const Event& other) const;
    friend std::ostream& operator<< (std::ostream &out, const Event &e);

  private:
    unsigned int n_features_;
    std::vector<double> weights_;
    std::vector<double> features_;
  };

  struct PGMLINK_EXPORT EventsStatistics {
  EventsStatistics() : n_total(0), n_mov(0), n_div(0), n_app(0), n_dis(0) {}
    EventsStatistics& operator+=( const EventsStatistics& rhs ) {
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
