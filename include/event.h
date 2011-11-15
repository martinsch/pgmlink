#ifndef EVENT_H
#define EVENT_H
#include <vector>
#include <ostream>

namespace Tracking {
   struct Event {
	Event() : energy(0) { type = Void; };

	enum EventType {Move, Division, Appearance, Disappearance, Void};
	EventType type;

	/**
	* The traxels participating in the event.
	*
	* A (Dis-)appearance is described by one traxel id, a move by two traxel ids and
	* a division by three (motherId, daughterId1, daughterId2). A Void is just an empty vector.
	*/
	std::vector<unsigned int> traxel_ids;

	/**
	* The energy assigned to this Event happening.
	*/
	double energy;

	// Beware: doesn't compare by energy!
	bool operator==(const Event& other) const;
	bool operator!=(const Event& other) const;
	friend std::ostream& operator<< (std::ostream &out, const Event &e);
    };
} /* namespace Tracking */

#endif /* EVENT_H */
