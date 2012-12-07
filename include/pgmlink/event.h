/**
   @file
   @ingroup tracking
   @brief tracking events like Move and Division
*/

#ifndef EVENT_H
#define EVENT_H
#include <vector>
#include <ostream>

namespace pgmlink {
  class Event {
  public:
  Event() : n_features_(0) { type = Void; };
    
	 enum EventType {Move, Division, Appearance, Disappearance, Merger, Void};
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
} /* namespace pgmlink */

#endif /* EVENT_H */
