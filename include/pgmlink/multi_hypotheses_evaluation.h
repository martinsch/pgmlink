#ifndef MULTI_HYPOTHESES_EVALUATION_H
#define MULTI_HYPOTHESES_EVALUATION_H

// stl
#include <stdexcept>
#include <map>
#include <utility>

// vigra
#include <vigra/multi_array.hxx>
#include <vigra/multi_iterator_coupled.hxx>

namespace pgmlink {
typedef unsigned LabelType;
typedef unsigned IntersectCountType;
template <typename T=LabelType, typename U=LabelType, typename V=IntersectCountType>
struct IntersectCountMap {
  typedef std::map<T, std::map<U, V> > type;
};

template <int N, typename T, typename U>
typename IntersectCountMap<T, U>::type get_intersect_count(vigra::MultiArrayView<N, T>, vigra::MultiArrayView<N, U>);

template <typename T, typename U>
std::pair<IntersectCountType, IntersectCountType> calculate_intersect_union(const typename IntersectCountMap<T,U,IntersectCountType>::type&,
                                                                            T,
                                                                            U,
                                                                            IntersectCountType,
                                                                            IntersectCountType);


/* IMPLEMENTATION */
template <int N, typename T, typename U>
typename IntersectCountMap<T, U>::type get_intersect_count(vigra::MultiArrayView<N, T> image1, vigra::MultiArrayView<N, U> image2) {
  if (image1.shape() != image2.shape()) {
    throw std::runtime_error("shape mismatch!");
  }
  typename IntersectCountMap<T, U>::type counts;
  typedef typename vigra::CoupledIteratorType<N, T, U>::type Iterator;
  Iterator start = vigra::createCoupledIterator(image1, image2);
  Iterator end = start.getEndIterator();
  for (Iterator it = start; it != end; ++start) {
    if (it.get<1>() != T()) {
      if(it.get<2>() != U()) {
        counts[it.get<1>()][it.get<2>()] += 1;
      }
    }
  }
  return counts;
}


template <typename T, typename U>
std::pair<IntersectCountType, IntersectCountType> calculate_intersect_union(const typename IntersectCountMap<T,U>::type& intersect_counts,
                                                                            T region1_label,
                                                                            U region2_label,
                                                                            IntersectCountType region1_count,
                                                                            IntersectCountType region2_count) {
  typename IntersectCountMap<T,U>::type::const_iterator it = intersect_counts.find(region1_label);
  if (it == intersect_counts.end()) {
    return std::make_pair(0, 1);
  }
  typename std::map<U, IntersectCountType>::const_iterator it2 = it->second.find(region2_label);
  if (it2 == it->second.end()) {
    return std::make_pair(0, 1);
  }
  return std::make_pair(it2->second, region1_count + region2_count - it2->second);
}

} // namespace pgmlink


#endif /* MULTI_HYPOTHESES_EVALUATION_H */

