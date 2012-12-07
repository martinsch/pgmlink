#include "pgmlink/event.h"
#include <cassert>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace pgmlink {

  ///
  /// class Event
  ///
  double Event::energy() const {
    assert(weights_.size() == features_.size());
    double energy = 0;
    for(size_t i=0; i != weights_.size(); ++i) {
      energy += weights_[i] * features_[i];
    }
    return energy;
  }

  Event& Event::number_of_features( unsigned int n ) {
    n_features_ = n;
    weights_.resize(n, 0);
    features_.resize(n, 0);
    return *this;
  }

  Event& Event::features( const std::vector<double>& v ) {
    if(v.size() != n_features_) {
      stringstream msg;
      msg << "Event::features(): expected argument length [" << v.size() << " ] equal to number of features [ " << n_features_ << " ]";
      throw invalid_argument(msg.str());
    }
    features_ = v;
    return *this;
  }

  Event& Event::weights(const std::vector<double>& v ) {
    if(v.size() != n_features_) {
      stringstream msg;
      msg << "Event::features(): expected argument length [" << v.size() << " ] equal to number of features [ " << n_features_ << " ]";
      throw invalid_argument(msg.str());
    }
    weights_ = v;
    return *this;
  }

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
	    type = "Void";
	    break;
	case Event::Merger:
		type = "Merger";
		break;
	default:
	    type = "unknown";
	    break;
    }
    out << "(" << type << ", traxel_ids:";
    for(size_t i = 0; i < e.traxel_ids.size(); ++i) {
	out << " " << e.traxel_ids[i]; 
    }
    out << ", energy: " << e.energy() << ")";
    return out;
}

} /* namespace pgmlink */
