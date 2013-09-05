/**
   @file
   @ingroup tracking
   @brief nearest neighbors of traxels
*/


#ifndef NEAREST_NEIGHBORS_H
#define NEAREST_NEIGHBORS_H
#include <map>
#include <ANN/ANN.h>
#include <boost/shared_ptr.hpp>

#include "pgmlink/pgmlink_export.h"

namespace pgmlink {
    class Traxel;

    class PGMLINK_EXPORT NearestNeighborSearch {
	public:
	    template <typename InputIt>
	    NearestNeighborSearch( InputIt traxel_begin,
				   InputIt traxel_end);
	    ~NearestNeighborSearch();
	
	     /**
	      * Returns (traxel id, distance*distance) map.
	      */
	    std::map<unsigned int, double> knn_in_range( const Traxel& query, double radius, unsigned int knn );
	    unsigned int count_in_range( const Traxel& query, double radius );

	private:
	    /**
	     * Ctor helper: define points and association between traxels and points
	     */
	    template <typename InputIt>
	    void define_point_set( InputIt traxel_begin, InputIt traxel_end );
	    ANNpoint point_from_traxel( const Traxel& traxel );

	    std::map<unsigned int, unsigned int> point_idx2traxel_id_;
	
	    const int dim_;
	
	    ANNpointArray points_;
	    boost::shared_ptr<ANNkd_tree> kd_tree_;
    };

} /* namespace pgmlink */



/****
 Implementation
 ****/
#include <cassert>
#include <iterator>
#include <boost/scoped_array.hpp>
#include "pgmlink/traxels.h"



namespace pgmlink {
using namespace std;
using namespace boost;

template <typename InputIt>
NearestNeighborSearch::NearestNeighborSearch(InputIt traxel_begin, InputIt traxel_end) : 
  dim_(3), points_(NULL) {
  size_t size(distance(traxel_begin, traxel_end));

  if(size > 0) {
    this->define_point_set( traxel_begin, traxel_end );
    try {
	kd_tree_ = boost::shared_ptr<ANNkd_tree>( new ANNkd_tree( points_, size, dim_ ) );
    } catch(...) {
	annDeallocPts(points_);
	points_ = NULL;
	throw;
    }
  }
}

template <typename InputIt>
void NearestNeighborSearch::define_point_set( InputIt traxel_begin, InputIt traxel_end ) {
    // allocate memory for kd-tree nodes
    size_t traxel_number = distance(traxel_begin, traxel_end);
    points_ = annAllocPts( traxel_number, dim_ );
    if( points_ == NULL ){
      throw "Allocation of points for kd-tree failed";
    }

    // fill the nodes with coordinates
    try {
	point_idx2traxel_id_.clear();
	size_t i = 0;
	for( InputIt traxel = traxel_begin; traxel != traxel_end; ++traxel, ++i) {
	  ANNpoint point = points_[i];
	  assert(dim_ == 3);
	  point[0] = traxel->X();
	  point[1] = traxel->Y();
	  point[2] = traxel->Z();

	  // save point <-> traxel association
	  point_idx2traxel_id_[i] = traxel->Id;
	}
    } catch(...) {
	annDeallocPts(points_);
	points_ = NULL;
	point_idx2traxel_id_.clear();
	throw;
    }
}

} /* namespace pgmlink */

#endif
