#ifndef PGMLINK_UTIL_H
#define PGMLINK_UTIL_H

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

namespace Tracking {

  namespace indexsorter {
    template <typename pair_t>
      bool compare( pair_t lhs, pair_t rhs ){
      return (lhs.second < rhs.second);
    }

    template <typename iterator>
      void sort_indices(iterator first, iterator last, std::vector<size_t>& indices) {
      typedef std::pair< size_t, typename std::iterator_traits<iterator>::value_type> pair_t;
      std::vector< pair_t > seq;
      size_t idx = 0;
      for(iterator it = first; it!=last; ++it) {
	seq.push_back( pair_t(idx, *it) );
	idx += 1;
      }
      std::sort( seq.begin(), seq.end(), (bool(*)(pair_t,pair_t))compare );

      for( typename std::vector<pair_t>::iterator it = seq.begin(); it != seq.end(); ++it ){
      	indices.push_back( it->first );
      }
    }
  } /* namespace indexsorter */

} /* namespace Tracking */
#endif /* PGMLINK_UTIL_H */
