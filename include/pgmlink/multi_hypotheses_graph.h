#ifndef MULTI_HYPOTHESES_GRAPH
#define MULTI_HYPOTHESES_GRAPH

// stl
#include <string>
#include <vector>
#include <set>
#include <sstream>
#include <cassert>

// boost
#include <boost/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/export.hpp>

// lemon
#include <lemon/list_graph.h>
#include <lemon/maps.h>

// pgmlink
#include "pgmlink/graph.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/traxels.h"



namespace pgmlink {


typedef std::map<int, std::vector<std::vector<unsigned> > > ConflictMap;



class ClassifierStrategy;


  ////
  //// IterableEditableValueMap
  ////
  template <typename Graph, typename Key, typename Value>
  class IterableEditableValueMap : public lemon::IterableValueMap<Graph, Key, Value> {
  private:
  public:
    Value& get_value(const Key& key);
    explicit IterableEditableValueMap(const Graph& graph,
                                      const Value& value = Value());
  };

  template <typename Graph, typename Key, typename Value>
  IterableEditableValueMap<Graph, Key, Value>::IterableEditableValueMap(const Graph& graph,
                                                                        const Value& value) :
    lemon::IterableValueMap<Graph,Key, Value>(graph, value) {
    
  }

  template <typename Graph, typename Key, typename Value>
  Value& IterableEditableValueMap<Graph, Key, Value>::get_value(const Key& key) {
    return lemon::IterableValueMap<Graph, Key, Value>::Parent::operator[](key).value;
  }



////
//// node_move_features
////
struct node_move_features{};

template <typename Graph>
struct property_map<node_move_features, Graph> {
  typedef IterableEditableValueMap<Graph, typename Graph::Node, std::map<Traxel, std::map<Traxel, feature_array> > > type;
  static const std::string name;
};

template <typename Graph>
const std::string property_map<node_move_features, Graph>::name = "node_move_features";


////
//// node_division_features
////
struct node_division_features{};

template <typename Graph>
struct property_map<node_division_features, Graph> {
  typedef IterableEditableValueMap<Graph, typename Graph::Node, std::map<Traxel, std::map<std::pair<Traxel, Traxel>, feature_array> > > type;
  static const std::string name;
};

template <typename Graph>
const std::string property_map<node_division_features, Graph>::name = "node_division_features";


////
//// node_count_features
////
struct node_count_features{};

template <typename Graph>
struct property_map<node_count_features, Graph> {
  typedef IterableEditableValueMap<Graph, typename Graph::Node, feature_type> type;
  static const std::string name;
};

template <typename Graph>
const std::string property_map<node_count_features, Graph>::name = "node_count_features";


////
//// node_conflict_sets
////
struct node_conflict_sets{};

template <typename Graph>
struct property_map<node_conflict_sets, Graph> {
  typedef IterableEditableValueMap<Graph, typename Graph::Node, std::vector<std::vector<unsigned> > > type;
  static const std::string name;
};

template <typename Graph>
const std::string property_map<node_conflict_sets, Graph>::name = "node_conflict_sets";


////
//// node_connected_component
////
struct node_connected_component{};

template <typename Graph>
struct property_map<node_connected_component, Graph> {
  typedef lemon::IterableValueMap<Graph, typename Graph::Node, std::pair<int, unsigned> > type;
  static const std::string name;
};

template <typename Graph>
const std::string property_map<node_connected_component, Graph>::name = "node_connected_component";




class MultiHypothesesGraph : public HypothesesGraph {
 public:
  typedef property_map<node_traxel, base_graph>::type TraxelMap;
  typedef property_map<node_division_features, base_graph>::type DivisionFeatureMap; // delete this?
  typedef property_map<node_move_features, base_graph>::type MoveFeatureMap; // delete this?
  typedef property_map<node_count_features, base_graph>::type CountFeatureMap; // delete this?
  typedef property_map<node_conflict_sets, base_graph>::type ConflictSetMap; // delete this?
  typedef property_map<node_connected_component, base_graph>::type ConnectedComponentMap;


  MultiHypothesesGraph();

  void remove_traxel_features();

  // if classification is done w/o saving, this can be removed
  void add_classifier_features(ClassifierStrategy* move,
                               ClassifierStrategy* division,
                               ClassifierStrategy* count,
                               ClassifierStrategy* detection);

  void add_conflicts(boost::shared_ptr<ConflictMap > conflicts);

  const ConflictMap& get_conflicts() const;


  
 private:
  // boost serialize
  friend class boost::serialization::access;
  template< typename Archive >
  void save( Archive&, const unsigned int /*version*/ ) const;
  template< typename Archive >
  void load( Archive&, const unsigned int /*version*/ );
  BOOST_SERIALIZATION_SPLIT_MEMBER()

  unsigned maximum_timestep_;
  // store conflicts by timestep
  boost::shared_ptr<ConflictMap > conflicts_;
};




////
//// SingleTimestepTraxel_MultiHypothesesbuilder
////
class PGMLINK_EXPORT SingleTimestepTraxel_MultiHypothesesBuilder : public SingleTimestepTraxel_HypothesesBuilder {
protected:
  virtual HypothesesGraph* construct() const;
  virtual HypothesesGraph* add_nodes(HypothesesGraph* graph) const;

};


// lemon graph format (lgf) serialization
PGMLINK_EXPORT void write_lgf( const MultiHypothesesGraph&, std::ostream& os=std::cout,
                               bool with_n_traxel=false );
PGMLINK_EXPORT void read_lgf( MultiHypothesesGraph&, std::istream& is=std::cin,
                              bool with_n_traxel=false);


template< typename Archive >
void MultiHypothesesGraph::save( Archive& ar, const unsigned int /*version*/ ) const {
  LOG(logDEBUG) << "MultiHypothesesGraph::save entered";
  ar & maximum_timestep_;
  ar & timesteps_;

  bool with_n_traxel = false;
  try {
    this->get(node_traxel());
    with_n_traxel = true;
  } catch( std::runtime_error& e ) {}

  ar & with_n_traxel;

  std::string lgf;
  {
    std::stringstream ss;
    write_lgf(*this, ss, with_n_traxel);
    lgf = ss.str();
  }
  ar & lgf;
}

template< typename Archive >
void MultiHypothesesGraph::load( Archive& ar, const unsigned int /*version*/ ) {
  LOG(logDEBUG) << "MultiHypothesesGraph::load entered";
  ar & maximum_timestep_;
  ar & timesteps_;

  bool with_n_traxel;
  ar & with_n_traxel;

  std::string lgf;
  ar & lgf;
  {
    std::stringstream ss(lgf);
    read_lgf(*this, ss, with_n_traxel);
  }
}





 
/* IMPLEMENTATIONS */  
    
}


#endif /* MULTI_HYPOTHESES_GRAPH */




