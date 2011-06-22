#include "synthetic_data.h"

namespace Tracking {
    void fill(TraxelStore& ts, int n_traxels, int n_timesteps) {
	unsigned int id = 0;
	for(int t = 0; t<n_timesteps; ++t){
	    for(int n = 0; n<n_traxels; ++n, ++id) {
		add(ts, Traxel(id, t));
	    }
	}
    }
}
