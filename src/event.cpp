#include "event.h"
#include <string>

using namespace std;

namespace Tracking {

///
/// class Event
///
bool Event::operator==(const Event& other) const {
    bool same = false;
    if(
	other.type == type &&
	other.traxel_ids.size() == traxel_ids.size() &&
	equal(traxel_ids.begin(), traxel_ids.end(), other.traxel_ids.begin())
    ) {
	same = true;
    }
    return same;
}

bool Event::operator!=(const Event& other) const {
    return !(*this == other);
}

ostream& operator<< (ostream &out, const Event &e)
{
    string type;
    switch(e.type) {
	case Event::Move:
	    type = "Move";
	    break;
	case Event::Division:
	    type = "Division";
	    break;
	case Event::Appearance:
	    type = "Appearance";
	    break;
	case Event::Disappearance:
	    type = "Disappearance";
	    break;
	case Event::Void:
	    type = "Voide";
	    break;
	default:
	    type = "unknown";
	    break;
    }
    out << "(" << type << ", traxel_ids:";
    for(size_t i = 0; i < e.traxel_ids.size(); ++i) {
	out << " " << e.traxel_ids[i]; 
    }
    out << ", energy: " << e.energy << ")";
    return out;
}

} /* namespace Tracking */
