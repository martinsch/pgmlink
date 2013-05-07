#include "pgmlink/pgm.h"

using namespace std;

namespace pgmlink {
  namespace pgm {
    ////
    //// class FactorEntry
    ////
    void FactorEntry::set( OpengmModel::ValueType v ) {
      if( entry_ == NULL ) {
	throw runtime_error("FactorEntry: is invalid");
      }
      *entry_ = v;
    }

    OpengmModel::ValueType FactorEntry::get() const {
      if( entry_ == NULL ) {
	throw runtime_error("FactorEntryPtr: is invalid");
      }
      return *entry_;
    } 
  } /* namespace pgm */
} /* namespace pgmlink */
