#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

// stl
// temp
#include <iostream>
#include <string>
#include <sstream>

// boost
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/python.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

// vigra
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

// pgmlink
#include "pgmlink/multi_hypotheses_graph.h"
#include "pgmlink/multi_hypotheses_tracking.h"
#include "pgmlink/classifier_auxiliary.h"
#include "pgmlink/feature.h"
#include "pgmlink/log.h"
#include "pgmlink/multi_hypotheses_evaluation.h"

// temp
#include "pgmlink/pgm_multi_hypotheses.h"
#include "pgmlink/reasoner_multi_hypotheses.h"


using namespace pgmlink;
using namespace boost::python;

MultiHypothesesGraph::ContainedRegionsMap& addContainedRegionsMap(MultiHypothesesGraph* g) {
  g->add(node_regions_in_component());
  return g->get(node_regions_in_component());
}

MultiHypothesesGraph::ContainedRegionsMap& getContainedRegionsMap(MultiHypothesesGraph* g) {
  return g->get(node_regions_in_component());
}

inline object pass_through(object const& o) { return o; }

template <typename ITERABLE_VALUE_MAP>
struct IterableValueMap_ValueIterator {
  typedef ITERABLE_VALUE_MAP map_type;
  IterableValueMap_ValueIterator( const map_type& m ) {
    it_ = m.beginValue();
    end_ = m.endValue();
  }
  typename map_type::Value next() {
    if (it_ == end_) {
      PyErr_SetString(PyExc_StopIteration, "No more data.");
      boost::python::throw_error_already_set();
    }
    const typename map_type::Value& result = *it_;
    ++it_;
    return result;
  }

  static IterableValueMap_ValueIterator<map_type>
  values( const map_type& m ) {
    return IterableValueMap_ValueIterator<map_type>( m );
  }

  static void
  wrap( const char* python_name) {
    class_<IterableValueMap_ValueIterator<map_type> >(python_name, init<const map_type&>(args("iterable_value_map")) )
        .def("next", &IterableValueMap_ValueIterator::next)
        .def("__iter__", &pass_through)
        ;
  }

  typename map_type::ValueIt it_;
  typename map_type::ValueIt end_;
};


struct PyTrackingOptions {
  MultiHypothesesTracking::Options options;
  PyTrackingOptions() {
    options.weights["forbidden"] = 0.;
    options.weights["timeout"] = 1e+75;
    options.weights["gap"] = 0.01;
    options.weights["opportunity"] = 1000;

    options.forward_backward = true;
    options.with_maximum_arcs = false;
  }
  
  void set(const std::string& name, double value) {
    options.weights[name] = value;
  }

  void set_rf(const std::string& name, const vigra::RandomForest<double>& rf) {
    options.classifiers[name] = rf;
  }

  void set_path(const std::string& name, const std::string& path) {
    options.paths[name] = path;
  }

  void add_feature(const std::string& classifier_type, const std::string& feature_name, const std::string& traxel_feature_name) {
    options.feature_lists[classifier_type].push_back(std::make_pair(feature_name, traxel_feature_name));
  }

  void with_divisions(bool check) {
    options.with_divisions = check;
  }

  void with_constraints(bool check) {
    options.with_constraints = check;
  }

  void with_detection_vars(bool check) {
    options.with_detection_vars = check;
  }

  void with_classifiers(bool check) {
    options.with_classifiers = check;
  }

  void with_constant_classifiers(bool check) {
    options.with_constant_classifiers = check;
  }

  void with_maximal_conflict_cliques(bool check) {
    options.with_maximal_conflict_cliques = check;
  }

  void with_classifier_count_precomputed(bool check) {
    options.classifier_count_precomputed = check;
  }

  void forward_backward(bool check) {
    options.forward_backward = check;
  }

  void with_constant_classifier_fallback(bool check) {
    options.constant_classifier_fallback = check;
  }

  void with_hierarchical_count_factor(bool check) {
    options.hierarchical_count_factor = check;
  }

  void limit_outgoing_arcs(unsigned limit) {
    options.with_maximum_arcs = true;
    options.weights["arc_limit"] = limit;
  }

  void with_counting_incoming_factor(bool check) {
      options.counting_incoming_factor = check;
    }

  std::string sanity_check() {
    return std::string("To be implemented");
  }
  
};


class PyMultiHypothesesTracking {
 public:
  PyMultiHypothesesTracking(const PyTrackingOptions& options)
      : tracker_(MultiHypothesesTracking(options.options))
  {}
  std::vector<std::vector<Event> > operator()(MultiHypothesesTraxelStore& ts) {
    std::vector<std::vector<Event> > ret = *tracker_(ts);
    return ret;
  }
 private:
  MultiHypothesesTracking tracker_;
};

template <typename T>
class PySharedPtr {
  PySharedPtr() : ptr_(new T) {}
  PySharedPtr(boost::shared_ptr<T> other) : ptr_(other) {}
 private:
  boost::shared_ptr<T> ptr_;
};



struct traxelstore_pickle_suite : pickle_suite {
    static std::string getstate( const MultiHypothesesTraxelStore& g ) {
      std::stringstream ss;
      boost::archive::text_oarchive oa(ss);
      oa & g;
      return ss.str();
    }

    static void setstate( MultiHypothesesTraxelStore& g, const std::string& state ) {
      std::stringstream ss(state);
      boost::archive::text_iarchive ia(ss);
      ia & g;
    }
};


class FeatureExtractorCollection {
 public:
  FeatureExtractorCollection()
      : name_calculator_mapping_(AvailableCalculators::get()) {

  }
  

  void add_to_vector(std::string calculator_name, std::string feature_name) {
    std::map<std::string, boost::shared_ptr<FeatureCalculator> >::const_iterator it = name_calculator_mapping_.find(calculator_name);
    if (it == name_calculator_mapping_.end()) {
      throw std::runtime_error("Calculator \"" + calculator_name + "\" not available!");
    } else {
      extractors_.push_back(FeatureExtractor(it->second, feature_name));
      calculators_.push_back(it->second);
    }
  }


  std::vector<FeatureExtractor> get_extractors() {
    return extractors_;
  }

  
  std::vector<boost::shared_ptr<FeatureCalculator> > get_calculators() {
    return calculators_;
  }

  
 private:
  const std::map<std::string, boost::shared_ptr<FeatureCalculator> >& name_calculator_mapping_;
  std::vector<FeatureExtractor> extractors_;
  std::vector<boost::shared_ptr<FeatureCalculator> > calculators_;
};

// wrap feature calculation w/traxels
vigra::NumpyArray<2, feature_type> calculate_1t(const std::vector<FeatureExtractor>& extractors,
                                                const Traxel& t1) {
  feature_array res;
  for (std::vector<FeatureExtractor>::const_iterator ex = extractors.begin();
       ex != extractors.end();
       ++ex) {
    feature_array temp = ex->extract(t1);
    res.insert(res.end(), temp.begin(), temp.end());
    LOG(logDEBUG3) << "Calculated " << ex->name();
  }
  vigra::NumpyArray<2, feature_type> ret(vigra::Shape2(1, res.size()));
  std::copy(res.begin(), res.end(), ret.begin());
  return ret;
}
vigra::NumpyArray<2, feature_type> calculate_2t(const std::vector<FeatureExtractor>& extractors,
                                                const Traxel& t1,
                                                const Traxel& t2) {
  feature_array res;
  for (std::vector<FeatureExtractor>::const_iterator ex = extractors.begin();
       ex != extractors.end();
       ++ex) {
    feature_array temp = ex->extract(t1, t2);
    res.insert(res.end(), temp.begin(), temp.end());
    LOG(logDEBUG3) << "Calculated " << ex->name();
  }
  vigra::NumpyArray<2, feature_type> ret(vigra::Shape2(1, res.size()));
  std::copy(res.begin(), res.end(), ret.begin());
  return ret;
}


vigra::NumpyArray<2, feature_type> calculate_3t(const std::vector<FeatureExtractor>& extractors,
                                                const Traxel& t1,
                                                const Traxel& t2,
                                                const Traxel& t3) {
  feature_array res;
  for (std::vector<FeatureExtractor>::const_iterator ex = extractors.begin();
       ex != extractors.end();
       ++ex) {
    feature_array temp = ex->extract(t1, t2, t3);
    res.insert(res.end(), temp.begin(), temp.end());
    LOG(logDEBUG3) << "Calculated " << ex->name();
  }
  vigra::NumpyArray<2, feature_type> ret(vigra::Shape2(1, res.size()));
  std::copy(res.begin(), res.end(), ret.begin());
  return ret;
}


// wrap feature calculation w/o traxels
/* vigra::NumpyArray<2, feature_type> calculate_2f(const std::vector<boost::shared_ptr<FeatureCalculator> >& calculators,
                                                const feature_array& f1,
                                                const feature_array& f2) {
  feature_array ret;
  for (std::vector<boost::shared_ptr<FeatureCalculator> >::const_iterator it = calculators.begin();
       it != calculators.end();
       ++it) {
    feature_array res = (*it)->calculate(f1, f2);
    ret.insert(ret.end(), res.begin(), res.end());
  }
  vigra::NumpyArray<2, feature_type> features(vigra::Shape2(1, ret.size()));
  std::copy(ret.begin(), ret.end(), features.begin());
  return features;
}


vigra::NumpyArray<2, feature_type> calculate_2n(const std::vector<boost::shared_ptr<FeatureCalculator> >& calculators,
                                                const vigra::NumpyArray<2, feature_type> n1,
                                                const vigra::NumpyArray<2, feature_type> n2) {
  feature_array f1(n1.shape()[1]);
  feature_array f2(n2.shape()[1]);
  assert(f1.size() == f2.size());
  std::copy(n1.begin(), n1.end(), f1.begin());
  std::copy(n2.begin(), n2.end(), f2.begin());
  return calculate_2f(calculators, f1, f2);
}


vigra::NumpyArray<2, feature_type> calculate_3f(const std::vector<boost::shared_ptr<FeatureCalculator> >& calculators,
                                                const feature_array& f1,
                                                const feature_array& f2,
                                                const feature_array& f3) {
  feature_array ret;
  for (std::vector<boost::shared_ptr<FeatureCalculator> >::const_iterator it = calculators.begin();
       it != calculators.end();
       ++it) {
    feature_array res = (*it)->calculate(f1, f2, f3);
    ret.insert(ret.end(), res.begin(), res.end());
  }
  vigra::NumpyArray<2, feature_type> features(vigra::Shape2(1, ret.size()));
  std::copy(ret.begin(), ret.end(), features.begin());
  return features;
}


vigra::NumpyArray<2, feature_type> calculate_3n(const std::vector<boost::shared_ptr<FeatureCalculator> >& calculators,
                                                const vigra::NumpyArray<2, feature_type> n1,
                                                const vigra::NumpyArray<2, feature_type> n2,
                                                const vigra::NumpyArray<2, feature_type> n3) {
  feature_array f1(n1.shape()[1]);
  feature_array f2(n2.shape()[1]);
  feature_array f3(n3.shape()[1]);
  assert(f1.size() == f2.size());
  assert(f2.size() == f3.size());
  std::copy(n1.begin(), n1.end(), f1.begin());
  std::copy(n2.begin(), n2.end(), f2.begin());
  std::copy(n3.begin(), n3.end(), f3.begin());
  return calculate_3f(calculators, f1, f2, f3);
  } */



// evaluation
template <int N, typename T, typename U>
typename IntersectCountMap<T, U>::type py_get_intersect_count(vigra::NumpyArray<N, T> image1, vigra::NumpyArray<N, U> image2) {
  return get_intersect_count<N, T, U>(image1, image2);
}


template <typename T, typename U>
std::pair<IntersectCountType, IntersectCountType> py_calculate_intersect_union(const typename IntersectCountMap<T,U>::type& intersect_counts,
                                                                               T region1_label,
                                                                               U region2_label,
                                                                               IntersectCountType region1_count,
                                                                               IntersectCountType region2_count) {
  return calculate_intersect_union(intersect_counts, region1_label, region2_label, region1_count, region2_count);
}


void export_multi_hypotheses() {
  IterableValueMap_ValueIterator<MultiHypothesesGraph::ContainedRegionsMap>::wrap("ContainedRegionsMap_ValueIt");
  // IterableValueMap_ValueIterator<node_traxel>::wrap("NodeTraxelMap_ValueIt");

  class_<MultiHypothesesGraph::ContainedRegionsMap, boost::noncopyable>("ContainedRegionsMap",
                                                                        init<const MultiHypothesesGraph&>(args("multi_hypotheses_graph"))
                                                                        )
      .def("__getitem__",
           &MultiHypothesesGraph::ContainedRegionsMap::operator[],
           return_internal_reference<>())
      .def("get_value",
           &MultiHypothesesGraph::ContainedRegionsMap::get_value,
           return_internal_reference<>())
      .def("values", &IterableValueMap_ValueIterator<MultiHypothesesGraph::ContainedRegionsMap>::values)
      ;

  // handle function overloading
  MultiHypothesesGraph::Node (HypothesesGraph::*addnode1)(const int)
      = &MultiHypothesesGraph::add_node;
  MultiHypothesesGraph::Node (HypothesesGraph::*addnode2)(const std::vector<int>)
      = &MultiHypothesesGraph::add_node;

  class_<MultiHypothesesGraph, boost::noncopyable>("MultiHypothesesGraph")
      .def("addNode", addnode1)
      .def("addNode", addnode2)
      ;


  ////
  //// MultiHypothesesTraxelStore
  ////
  class_<MultiHypothesesTraxelStore>("MultiHypothesesTraxelStore")
      .def("add", &MultiHypothesesTraxelStore::add,
           return_internal_reference<>())
      .def("startComponent", &MultiHypothesesTraxelStore::start_component,
           return_internal_reference<>())
      .def("make_string", &MultiHypothesesTraxelStore::print)
      .def("__str__", &MultiHypothesesTraxelStore::print)
      .def("addConflictMap", &MultiHypothesesTraxelStore::add_conflict_map)
      .def("addConflictMap", &MultiHypothesesTraxelStore::add_signed_conflict_map)
      .def_pickle(traxelstore_pickle_suite())
      ;


  ////
  //// TrackingOptions
  ////
  class_<PyTrackingOptions, boost::noncopyable>("TrackingOptions", init<>())
      .def("set", &PyTrackingOptions::set, return_internal_reference<>())
      .def("set", &PyTrackingOptions::set_rf, return_internal_reference<>())
      .def("set", &PyTrackingOptions::set_path, return_internal_reference<>())
      .def("add", &PyTrackingOptions::add_feature, return_internal_reference<>())
      .def("withDivisions", &PyTrackingOptions::with_divisions, return_internal_reference<>())
      .def("withConstraints", &PyTrackingOptions::with_constraints, return_internal_reference<>())
      .def("withDetectionVars", &PyTrackingOptions::with_detection_vars, return_internal_reference<>())
      .def("withClassifiers", &PyTrackingOptions::with_classifiers, return_internal_reference<>())
      .def("withConstantClassifiers", &PyTrackingOptions::with_constant_classifiers, return_internal_reference<>())
      .def("withMaximalConflictCliques", &PyTrackingOptions::with_maximal_conflict_cliques, return_internal_reference<>())
      .def("withConstantClassifierFallback", &PyTrackingOptions::with_constant_classifier_fallback, return_internal_reference<>())
      .def("withClassifierCountPrecomputed", &PyTrackingOptions::with_classifier_count_precomputed, return_internal_reference<>())
      .def("withHierarchicalCountFactor", &PyTrackingOptions::with_hierarchical_count_factor, return_internal_reference<>())
      .def("withCountingIncomingFactor", &PyTrackingOptions::with_counting_incoming_factor, return_internal_reference<>())
      .def("forwardBackward", &PyTrackingOptions::forward_backward, return_internal_reference<>())
      .def("sanityCheck", &PyTrackingOptions::sanity_check)
      .def("limitOutgoingArcs", &PyTrackingOptions::limit_outgoing_arcs)
      ;


  ////
  //// MultiHypothesesTracking
  ////
  class_<PyMultiHypothesesTracking, boost::noncopyable>("MultiHypothesesTracker",
                                                        init<const PyTrackingOptions&>()[with_custodian_and_ward<1, 2>()]
                                                        )
      .def("track", &PyMultiHypothesesTracking::operator())
      .def("__call__", &PyMultiHypothesesTracking::operator())
      ;


  ////
  //// EventsPointer
  ////
  class_<boost::shared_ptr<std::vector<std::vector<Event> > > >("EventsPointer",
                                                                init<boost::shared_ptr<std::vector<std::vector<Event> > > >()
                                                                );


  ////
  //// FeatureExtractorCollection
  ////
  class_<std::vector<FeatureExtractor> >("FeatureExtractorVector");
  class_<std::vector<boost::shared_ptr<FeatureCalculator> > >("FeatureCalculatorVector");
  class_<FeatureExtractorCollection, boost::noncopyable>("FeatureExtractorCollection")
      .def("addToVector", &FeatureExtractorCollection::add_to_vector)
      .def("getExtractors", &FeatureExtractorCollection::get_extractors)
      .def("getCalculators", &FeatureExtractorCollection::get_calculators)
      ;

  def ("extractFeatures", &calculate_1t);
  def ("extractFeatures", &calculate_2t);
  def ("extractFeatures", &calculate_3t);
  /* def("calculateFeature", &calculate_2n);
  def("calculateFeature", &calculate_3n);
  def("calculateFeature", &calculate_2f);
  def("calculateFeature", &calculate_3f); */


  ////
  //// Evaluation
  ////
  class_<std::pair<IntersectCountType, IntersectCountType> >("CountPair")
      .def_readwrite("first", &std::pair<IntersectCountType, IntersectCountType>::first)
      .def_readwrite("second", &std::pair<IntersectCountType, IntersectCountType>::second)
      ;

  class_<IntersectCountMap<> >("IntersectCountMap");

  def("getIntersectCount", &py_get_intersect_count<2, LabelType, LabelType>);
  def("getIntersectCount", &py_get_intersect_count<3, LabelType, LabelType>);
  def("getIntersectCount", &py_get_intersect_count<4, LabelType, LabelType>);

  def("calculateIntersectUnion", &py_calculate_intersect_union<LabelType, LabelType>);

}
  



