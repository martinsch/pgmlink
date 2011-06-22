#ifndef TRACKLETS_H
#define TRACKLETS_H

#include <map>
#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <ostream>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/composite_key.hpp>

namespace Tracking {
//
// feature data structures
//
typedef float feature_type;
typedef std::vector<feature_type> feature_array;
typedef std::map<std::string,feature_array> FeatureMap;



//
// retrieve spatial coordinates from features
//
class Locator {
 public:
  Locator(std::string fn) : feature_name_(fn) {};
  virtual Locator* clone() = 0;
  virtual ~Locator() {};

  virtual bool is_applicable(const FeatureMap& m) const { return m.count(feature_name_)==1; };
  virtual double X(const FeatureMap&) const = 0;
  virtual double Y(const FeatureMap&) const = 0;
  virtual double Z(const FeatureMap&) const = 0;
 protected:
  std::string feature_name_;
  double coordinate_from(const FeatureMap&, size_t idx) const;
};

class ComLocator : public Locator {
 public:
 ComLocator() : Locator("com") {};
  virtual ComLocator* clone() { return new ComLocator(*this); };
  double X(const FeatureMap& m) const { return coordinate_from(m, 0); };
  double Y(const FeatureMap& m) const { return coordinate_from(m, 1); };
  double Z(const FeatureMap& m) const { return coordinate_from(m, 2); };
};

class IntmaxposLocator : public Locator {
 public:
 IntmaxposLocator() : Locator("intmaxpos") {};
  virtual IntmaxposLocator* clone() { return new IntmaxposLocator(*this); };
  double X(const FeatureMap& m) const { return coordinate_from(m, 1); };
  double Y(const FeatureMap& m) const { return coordinate_from(m, 2); };
  double Z(const FeatureMap& m) const { return coordinate_from(m, 3); };  
};



//
// Traxel datatype
//
 class Traxel {
 public:
   // construction / assignment
   //takes ownership of locator pointer
   Traxel(unsigned int id = 0, int timestep = 0, FeatureMap fmap = FeatureMap(), Locator* l = new ComLocator()) :
	   Id(id), Timestep(timestep), features(fmap), locator_(l) {};
   Traxel(const Traxel& other);
   Traxel& operator=(const Traxel& other);
   ~Traxel() { delete locator_; };
   const Traxel& set_locator(Locator*);
   
   // fields
   unsigned int Id; // id of connected component (aka "label")
   int Timestep; // traxel occured after
   FeatureMap features;
   
   // position according to locator
   double X() const;
   double Y() const;
   double Z() const;

   // relation to other traxels
   double distance_to(const Traxel& other) const;
   double angle(const Traxel& leg1, const Traxel& leg2) const;
   friend std::ostream& operator<< (std::ostream &out, const Traxel &t);
 private:
   Locator* locator_;
 };
 // compare by (time,id) (Traxels can be used as keys (for instance in a std::map) )
 bool operator<(const Traxel& t1, const Traxel& t2);
 bool operator>(const Traxel& t1, const Traxel& t2);
 

 //
 // Traxel collections
 //
 
 // Traxel map
 typedef std::map<unsigned int, Traxel> Traxels;
 template<typename InputIterator>
   Traxels traxel_map_from_traxel_sequence(InputIterator begin, InputIterator end);
 
 // Traxel key-value store
 using namespace boost;
 using namespace boost::multi_index;

 //tags
 struct by_timestep {};
 struct by_timeid {};

 typedef multi_index_container<
   Traxel,
     indexed_by<
       ordered_non_unique<tag<by_timestep>, member<Traxel,int,&Traxel::Timestep> >,
       hashed_unique<
	   tag<by_timeid>,
	   composite_key<
	       Traxel,
	       member<Traxel, int, &Traxel::Timestep>,
	       member<Traxel, unsigned int, &Traxel::Id>
	   > 
	> 
     >
   > 
   TraxelStore;
 typedef TraxelStore::index<by_timestep>::type
   TraxelStoreByTimestep;
 typedef TraxelStore::index<by_timeid>::type
   TraxelStoreByTimeid;

 // TraxelStore functions 
 std::set<TraxelStoreByTimestep::key_type>
   timesteps(const TraxelStore&);

 TraxelStore& add(TraxelStore&, const Traxel&);

 template<typename InputIt>
   TraxelStore& add(TraxelStore&, InputIt begin, InputIt end);

 std::vector<std::vector<Traxel> > nested_vec_from(const TraxelStore&);



/**/
/* implementation */
/**/
template<typename InputIterator>
Traxels traxel_map_from_traxel_sequence(InputIterator begin, InputIterator end) {
    Traxels ret;
    for(;begin != end; ++begin) {
	ret[begin->Id] = *begin;
    }
    return ret;
}

 template<typename InputIt>
   TraxelStore& add(TraxelStore& ts, InputIt begin, InputIt end) {
   ts.get<by_timestep>().insert(begin, end);
   return ts;
 }

} /* namespace Tracking */


#endif /* TRACKLETS_H */

