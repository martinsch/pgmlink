#include <map>
#include <cassert>
#include <iterator>
#include <ANN/ANN.h>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <pgmlink/nearest_neighbors.h>
#include "pgmlink/traxels.h"

namespace pgmlink{
using namespace std;
using namespace boost;


NearestNeighborSearch::~NearestNeighborSearch() {
	if(points_ != NULL) {
	    annDeallocPts(points_);
	    points_ = NULL;
	}
}



map<unsigned int, double> NearestNeighborSearch::knn_in_range( const Traxel& query, double radius, unsigned int knn , const bool reverse) {
    if( radius < 0 ) {
	throw "knn_in_range: radius has to be non-negative.";
    }

    map<unsigned int, double> return_value;

    // empty search space?
    if(points_ == NULL && kd_tree_.get() == NULL) {
      return return_value;
    }

    // allocate
    ANNpoint query_point( NULL );
    query_point = this->point_from_traxel(query, reverse);
    if( query_point == NULL ) {
	throw "query point allocation failure";
    }

    // search
    try {
	scoped_array<ANNidx> nn_indices( new ANNidx[knn] );
	scoped_array<ANNdist> nn_distances( new ANNdist[knn] );

	const int points_in_range = kd_tree_->annkFRSearch( query_point, radius*radius, knn,
                                         nn_indices.get(), nn_distances.get());

	if( points_in_range < 0 ) {
	    throw "knn_in_range: ANN search return negative number of nearest neighbors";
	}

	// construct return value
	
	// there may be less points in range, than nearest neighbors demanded
	const int actual = (static_cast<int>(knn) < points_in_range) ? knn : points_in_range;
	for( ANNidx i = 0; i < actual; ++i) {
	    return_value[ point_idx2traxel_id_[nn_indices[i]] ] = nn_distances[i];
	}
    } catch(...) {
	if( query_point != NULL) {
	    annDeallocPt( query_point );
	}
	throw;
    }

    // clean up
    if( query_point != NULL) {
	annDeallocPt( query_point );
    }

    return return_value;
}


std::map<unsigned int, double> NearestNeighborSearch::knn( const Traxel& query, unsigned int knn, const bool reverse ) {
   map<unsigned int, double> return_value;

    // empty search space?
    if(points_ == NULL && kd_tree_.get() == NULL) {
      return return_value;
    }

    // allocate
    ANNpoint query_point( NULL );
    query_point = this->point_from_traxel(query, reverse);
    if( query_point == NULL ) {
	throw "query point allocation failure";
    }

    // search
    try {
	scoped_array<ANNidx> nn_indices( new ANNidx[knn] );
	scoped_array<ANNdist> nn_distances( new ANNdist[knn] );

	kd_tree_->annkSearch( query_point, knn,
                              nn_indices.get(), nn_distances.get());

	// construct return value
	
	// there may be less points in range, than nearest neighbors demanded
	for( ANNidx i = 0; i < knn; ++i) {
	    return_value[ point_idx2traxel_id_[nn_indices[i]] ] = nn_distances[i];
	}
    } catch(...) {
	if( query_point != NULL) {
	    annDeallocPt( query_point );
	}
	throw;
    }

    // clean up
    if( query_point != NULL) {
	annDeallocPt( query_point );
    }

    return return_value;
}



unsigned int NearestNeighborSearch::count_in_range( const Traxel& query, double radius , const bool reverse) {
    if( radius < 0 ) {
	throw "count_in_range: radius has to be non-negative.";
    }

    // empty search space?
    if(points_ == NULL && kd_tree_.get() == NULL) {
      return 0;
    }

    // allocate
    ANNpoint query_point( NULL );
    query_point = this->point_from_traxel(query, reverse);
    if( query_point == NULL ) {
	throw "query point allocation failure";
    }

    // search
    int points_in_range;
    try {
	// search with 0 nearest neighbors -> returns just a range count
	points_in_range = kd_tree_->annkFRSearch( query_point, radius*radius, 0 );

	if( points_in_range < 0 ) {
	    throw "knn_in_range: ANN search return negative number of nearest neighbors";
	}
    } catch(...) {
	if( query_point != NULL) {
	    annDeallocPt( query_point );
	}
	throw;
    }

    // clean up
    if( query_point != NULL) {
	annDeallocPt( query_point );
    }

    return points_in_range;
}




ANNpoint NearestNeighborSearch::point_from_traxel( const Traxel& traxel , const bool reverse) {
    assert(dim_ == 3);
    ANNpoint point = annAllocPt( dim_ );

    if (reverse) {
        point[0] = traxel.X();
            point[1] = traxel.Y();
            point[2] = traxel.Z();
            LOG(logDEBUG4) << "NearestNeighborSearch::point_from_traxel (reverse): " << traxel <<
                    " point = " << point[0] << "," << point[1] << "," << point[2];
    } else {
    	point[0] = traxel.X_corr();
		point[1] = traxel.Y_corr();
		point[2] = traxel.Z_corr();
		LOG(logDEBUG4) << "NearestNeighborSearch::point_from_traxel: " << traxel <<
		                " point = " << point[0] << "," << point[1] << "," << point[2];
    }return point;
}


}
