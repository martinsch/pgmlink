#ifndef CLASSIFIER_AUXILIARY_H
#define CLASSIFIER_AUXILIARY_H

// stl
#include <vector>
#include <string>
#include <algorithm>

// boost
#include <boost/shared_ptr.hpp>

// vigra

// pgmlink
#include "pgmlink/feature.h"
#include "pgmlink/traxels.h"

namespace pgmlink {
////
//// class FeatureCalculator
////
class FeatureCalculator {
 public:
  static const std::string name_;
  static const unsigned length;

  virtual ~FeatureCalculator();
  virtual feature_array calculate(const feature_array& f1) const;
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const;
  virtual const std::string& name() const;

  bool operator==(const FeatureCalculator& other);
};


class IdentityCalculator : public FeatureCalculator {
 public:
  static const std::string name_;
  static const unsigned length;

  virtual ~IdentityCalculator();
  virtual feature_array calculate(const feature_array& f1) const;
  virtual const std::string& name() const;

};


class SquaredDifferenceCalculator : public FeatureCalculator {
 public:
  static const std::string name_;
  static const unsigned length;

  virtual ~SquaredDifferenceCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual const std::string& name() const;
};


class AbsoluteDifferenceCalculator : public FeatureCalculator {
 public:
  static const std::string name_;
  static const unsigned length;

  virtual ~AbsoluteDifferenceCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual const std::string& name() const;
};


class SquareRootSquaredDifferenceCalculator : public FeatureCalculator {
 public:
  static const std::string name_;
  static const unsigned length;

  virtual ~SquareRootSquaredDifferenceCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual const std::string& name() const;
};


class RatioCalculator : public FeatureCalculator {
 public:
  static const std::string name_;
  static const unsigned length;

  virtual ~RatioCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const;
  virtual const std::string& name() const;
};


////
//// class ParentRatios
////
template <typename Modifier>
class ParentRatios : public FeatureCalculator {
 public:
  static const std::string name_;
  static const unsigned length;

  virtual ~ParentRatios();
  feature_array calculate(const feature_array& parent, const feature_array& child1, const feature_array& child2) const;
  virtual const std::string& name() const;
};


////
//// class ParentSquaredDifferences
////
template <typename Modifier>
class ParentSquaredDifferences : public FeatureCalculator {
 public:
  static const std::string name_;
  static const unsigned length;

  virtual ~ParentSquaredDifferences();
  feature_array calculate(const feature_array& parent, const feature_array& child1, const feature_array& child2) const;
  virtual const std::string& name() const;
};


////
//// class MaxParentRatio
////
class Max {
 public:
  feature_type operator()(feature_array::const_iterator begin, feature_array::const_iterator end) {
    return *std::max_element(begin, end);
  }
};
typedef ParentRatios<Max> MaxParentRatio;

////
//// class MinParentRatio
////
class Min {
 public:
  feature_type operator()(feature_array::const_iterator begin, feature_array::const_iterator end) {
    return *std::min_element(begin, end);
  }
};
typedef ParentRatios<Min> MinParentRatio;


////
//// class MeanParentRatio
////
class Mean {
 public:
  feature_type operator()(feature_array::const_iterator begin, feature_array::const_iterator end) {
    assert (begin != end);
    return std::accumulate(begin, end, 0)/(end-begin);
  }
};
typedef ParentRatios<Mean> MeanParentRatio;


////
//// class MaxParentSquaredDifference
////
typedef ParentSquaredDifferences<Max> MaxParentSquaredDifference;


////
//// class MinParentSquaredDifference
////
typedef ParentSquaredDifferences<Min> MinParentSquaredDifference;


////
//// class MeanParentSquaredDifference
////
typedef ParentSquaredDifferences<Mean> MeanParentSquaredDifference;



////
//// class RatioParentSquaredDifference
class Ratio {
 public:
  feature_type operator()(feature_array::const_iterator begin, feature_array::const_iterator end) {
    assert(end - begin == 2);
    feature_array::const_iterator second = end - 1;
    if (*begin == *second) {
      return 1.;
    } else {
      return *begin < *second ? *begin/(*second) : (*second)/(*begin);
    }
  }
};
typedef ParentSquaredDifferences<Ratio> RatioParentSquaredDifference;





/* class IntensityRatioCalculator : public FeatureCalculator {
 public:
  static const std::string name_;
  static const unsigned length;

  virtual ~IntensityRatioCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2) const;
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const;
  virtual const std::string& name() const;
}; */


/* class ChildrenMeanParentIntensityRatioCalculator : public FeatureCalculator {
  public:
  static const std::string name_;
  static const unsigned length;

  virtual ~ChildrenMeanParentIntensityRatioCalculator();
  virtual feature_array calculate(const feature_array& f1, const feature_array& f2, const feature_array& f3) const;
  virtual const std::string& name() const;
}; */
  

////
//// class FeatureExtractor
////
class FeatureExtractor {
 public:
  FeatureExtractor(boost::shared_ptr<FeatureCalculator> calculator, const std::string& feature_name);
  feature_array extract(const Traxel& t1) const;
  feature_array extract(const Traxel& t1, const Traxel& t2) const;
  feature_array extract(const Traxel& t1, const Traxel& t2, const Traxel& t3) const;
  boost::shared_ptr<FeatureCalculator> calculator() const;
  std::string name() const;

 private:
  boost::shared_ptr<FeatureCalculator> calculator_;
  std::string feature_name_;
};


class AvailableCalculators {
 public:
  static const std::map<std::string, boost::shared_ptr<FeatureCalculator> >& get();
 private:
  static std::map<std::string, boost::shared_ptr<FeatureCalculator> > features_;
};


/* IMPLEMENTATIONS */


////
//// class ParentRatios
////
template <typename Modifier>
const std::string ParentRatios<Modifier>::name_ = "ParentRatio";

template <typename Modifier>
const unsigned ParentRatios<Modifier>::length = 1;

template <typename Modifier>
ParentRatios<Modifier>::~ParentRatios() {

}

template <typename Modifier>
feature_array ParentRatios<Modifier>::calculate(const feature_array& parent, const feature_array& child1, const feature_array& child2) const {
  feature_array ratios(2, 0.);
  if (child1[0] == parent[0]) {
    ratios[0] = 1.;
  } else {
    ratios[0] = child1[0] < parent[0] ? child1[0]/parent[0] : parent[0]/child1[0];
  }
  if (child2[0] == parent[0]) {
    ratios[1] = 1.;
  } else {
    ratios[1] = child2[0] < parent[0] ? child2[0]/parent[0] : parent[0]/child2[0];
  }
  Modifier mod;
  return feature_array(1, mod(ratios.begin(), ratios.end()));
}


template <typename Modifier>
const std::string& ParentRatios<Modifier>::name() const {
  return ParentRatios<Modifier>::name_;
}


////
//// class ParentSquaredDifferences
////
template <typename Modifier>
const std::string ParentSquaredDifferences<Modifier>::name_ = "ParentSquaredDifferencesRatio";

template <typename Modifier>
const unsigned ParentSquaredDifferences<Modifier>::length = 1;

template <typename Modifier>
ParentSquaredDifferences<Modifier>::~ParentSquaredDifferences() {

}

template <typename Modifier>
feature_array ParentSquaredDifferences<Modifier>::calculate(const feature_array& parent, const feature_array& child1, const feature_array& child2) const {
  feature_array dist(2, 0.);
  SquaredDistance distance;
  dist[0] = distance(parent, child1);
  dist[1] = distance(parent, child2);
  Modifier mod;
  return feature_array(1, mod(dist.begin(), dist.end()));
}


template <typename Modifier>
const std::string& ParentSquaredDifferences<Modifier>::name() const {
  return ParentSquaredDifferences<Modifier>::name_;
}


} /* namesapce pgmlink */

#endif /* CLASSIFIER_H */
