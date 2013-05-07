#pragma once
#ifndef OPENGM_FEATURE_HXX
#define OPENGM_FEATURE_HXX

#include <algorithm>
#include <set>
#include "pgmlink/ext_opengm/decorator_weighted.hxx"
#include "opengm/opengm.hxx"

namespace opengm {
namespace learning {

template<class FUNCTION> class ObservableWeight;

template<class FUNCTION>
class WeightObserver : public FunctionDecoratorWeighted<FUNCTION> {
public:
  typedef FUNCTION FunctionType;
  WeightObserver(FunctionType*, ObservableWeight<FUNCTION>* w = NULL);
  WeightObserver(const WeightObserver<FUNCTION>&);
  WeightObserver& operator=( WeightObserver<FUNCTION> );

  template <class FUNC>
  friend void swap(WeightObserver<FUNC>&, WeightObserver<FUNC>&);
  ~WeightObserver();

  void observable(ObservableWeight<FUNCTION>*);
  ObservableWeight<FUNCTION>* observable();

private:
  ObservableWeight<FUNCTION>* observable_;
};

template<class FUNCTION>
class ObservableWeight {
public:
  typedef FUNCTION FunctionType;
  typedef typename FUNCTION::ValueType ValueType;

  ObservableWeight( ValueType weight=1. ) : weight_(weight) {};
  ~ObservableWeight();

  void registerObserver(WeightObserver<FUNCTION>*);
  void removeObserver(WeightObserver<FUNCTION>*);
  bool hasObserver(WeightObserver<FUNCTION>*);

  void weight( const ValueType& );
  ValueType weight() const;
  void notifyObservers() const;

private:
  ObservableWeight( const ObservableWeight<FUNCTION>& ) {};
  ObservableWeight& operator=(ObservableWeight<FUNCTION>) {};

  std::set<WeightObserver<FUNCTION>*> observers_;
  ValueType weight_;
};


/**/
/* Implementation */
/**/

////
//// class WeightObserver
////
template<class FUNCTION>
WeightObserver<FUNCTION>::WeightObserver(FunctionType* f,
					 ObservableWeight<FUNCTION>* observable)
  : FunctionDecoratorWeighted<FunctionType>(f), observable_(observable) {
  // register itself
  if(observable_) {
    observable_->registerObserver( this );
  }
}

template <class FUNCTION>
WeightObserver<FUNCTION>::WeightObserver( const WeightObserver<FUNCTION>& other )
  : FunctionDecoratorWeighted<FUNCTION>( other ), observable_(other.observable_) {
  if(observable_) {
    observable_->registerObserver( this );
  }
}

template <class FUNCTION>
WeightObserver<FUNCTION>& WeightObserver<FUNCTION>::operator=( WeightObserver<FUNCTION> other ) {
  swap(*this, other);
  return *this;
}

template <class FUNC>
void swap(WeightObserver<FUNC>& lhs, WeightObserver<FUNC>& rhs) {
  if(lhs.observable_) {
    lhs.observable_->removeObserver(&lhs);
  }
  if(rhs.observable_) {
    rhs.observable_->removeObserver(&rhs);
  }

  using std::swap; // enable ADL
  swap(static_cast<FunctionDecoratorWeighted<FUNC> >(lhs),static_cast<FunctionDecoratorWeighted<FUNC> >(rhs));
  swap(lhs.observable_, rhs.observable_);

  if(lhs.observable_) {
    lhs.observable_->registerObserver(&lhs);
  }
  if(rhs.observable_) {
    rhs.observable_->registerObserver(&rhs);
  }
}

template<class FUNCTION>
WeightObserver<FUNCTION>::~WeightObserver() {
  if(observable_) {
    observable_->removeObserver(this);
  }
  observable_=NULL;
}

template<class FUNCTION>
void WeightObserver<FUNCTION>::observable(ObservableWeight<FUNCTION>* o) {
  if(observable_ && o) {
    throw RuntimeError("WeightObserver::observable(): instance has already an observable");
  }
  observable_=o;
  if(observable_) {
    observable_->registerObserver(this);
  }
}

template<class FUNCTION>
inline ObservableWeight<FUNCTION>* WeightObserver<FUNCTION>::observable() {
  return observable_;
}

////
//// class ObservableWeight
////
template<class FUNCTION>
void ObservableWeight<FUNCTION>::registerObserver(WeightObserver<FUNCTION>* o) {
  if(o==NULL) {
    throw RuntimeError("ObservableWeight::registerObserver(): observer is NULL");
  }
  observers_.insert(o);
}

template<class FUNCTION>
ObservableWeight<FUNCTION>::~ObservableWeight() {
  for(typename std::set<WeightObserver<FUNCTION>*>::iterator it=observers_.begin(); it!=observers_.end(); ++it) {
    (*it)->observable(NULL);
  }
}

template<class FUNCTION>
void ObservableWeight<FUNCTION>::removeObserver(WeightObserver<FUNCTION>* o) {
  if(o) {
    typename std::set<WeightObserver<FUNCTION>*>::iterator e = observers_.find(o);
    if(e != observers_.end()) {
      (*e)->observable(NULL);
      observers_.erase(e);
    }
  }
}

template<class FUNCTION>
bool ObservableWeight<FUNCTION>::hasObserver(WeightObserver<FUNCTION>* o) {
  if(o) {
    return observers_.count(o);
  } else {
    throw RuntimeError("ObservableWeight::hasObserver(): observer is NULL");
  }
}

template<class FUNCTION>
inline void ObservableWeight<FUNCTION>::weight( const ValueType& w ) {
  weight_ = w;
  this->notifyObservers();
}

template<class FUNCTION>
inline typename ObservableWeight<FUNCTION>::ValueType ObservableWeight<FUNCTION>::weight() const {
  return weight_;
}

template<class FUNCTION>
void ObservableWeight<FUNCTION>::notifyObservers() const {
  for(typename std::set<WeightObserver<FUNCTION>*>::iterator it=observers_.begin(); it!=observers_.end(); ++it) {
    (*it)->weight(weight_);
  }
}


}} // end namespaces

#endif
