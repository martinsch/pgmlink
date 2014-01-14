// stl
#include <vector>
#include <map>
#include <utility>

// vigra
#include <vigra/multi_array.hxx> // vigra::MultiArrayView
#include <vigra/multi_iterator_coupled.hxx> // vigra::CoupledScanOrderIterator

// pgmlink
#include <pgmlink/traxels.h> // pgmlink::feature_array
#include <pgmlink/log.h> // LOG

namespace pgmlink {
template <int N, typename T>
void extract_intersects( const vigra::MultiArrayView<N, T> arr1,
                         const vigra::MultiArrayView<N, T> arr2,
                         std::map<T, std::map<T, unsigned> >& intersects);


}

// --------- IMPLEMENTATION  --------- //

namespace pgmlink {
template <int N, typename T>
void extract_intersects( const vigra::MultiArrayView<N, T> arr1,
                         const vigra::MultiArrayView<N, T> arr2,
                         std::map<T, std::map<T, unsigned> >& intersects) {
  LOG(logDEBUG) << "extract_intersects() -- entered";
  typedef typename vigra::CoupledIteratorType<N, T, T>::type CoupledIterator;
  CoupledIterator start = vigra::createCoupledIterator(arr1, arr2);
  CoupledIterator end = start.getEndIterator();
  for (CoupledIterator it = start; it != end; ++it) {
    const T& v1 = it.get<1>();
    const T& v2 = it.get<2>();
    if (v1 != 0 && v2 != 0) {
      intersects[v1][v2] += 1;
    }
  }
}

}
