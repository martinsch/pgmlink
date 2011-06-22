#include <traxels.h>
#include <synthetic_data.h>
#include <track.h>
#include <event.h>
#include <iostream>

using namespace Tracking;
using namespace std;

void print_traxels(std::vector< std::vector<Traxel> > traxels) {
     for(size_t i = 0; i < traxels.size(); ++i){
	cout << "timestep " << i << ":\n";
	for(size_t j=0; j < traxels[i].size(); j++) {
	    cout << "\t" << traxels[i][j] << endl;
	}
    }
}

void print_events(std::vector< std::vector<Event> > events) {
     for(size_t i = 0; i < events.size(); ++i){
	cout << "timestep " << i << ":\n";
	for(size_t j=0; j < events[i].size(); j++) {
	    cout << "\t" << events[i][j] << endl;
	}
    }
}

int main() {
    TraxelStore ts;
    cout << "Filling traxelstore...";
    fill(ts, 10, 100);
    cout << "done!" << endl;

    vector<vector<Traxel> > nested;
    nested = nested_vec_from(ts);
    //print_traxels(nested);
    MultiTrack track;
    vector<vector<Event> > output = track(ts);
    print_events(output);

    return 0;
}
