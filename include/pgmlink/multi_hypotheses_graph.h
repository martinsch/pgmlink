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
#include <boost/serialization/utility.hpp>
#include <boost/serialization/shared_ptr.hpp>

// lemon
#include <lemon/list_graph.h>
#include <lemon/maps.h>

// pgmlink
#include "pgmlink/graph.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/traxels.h"



namespace pgmlink {

typedef std::vector<unsigned> ConflictSet;

typedef std::vector<ConflictSet > ConflictSetVector;

typedef std::map<int, ConflictSetVector > ConflictMap;

class MultiHypothesesGraph;

typedef boost::shared_ptr<MultiHypothesesGraph> MultiHypothesesGraphPtr;



class ClassifierStrategy;


  

////
//// node_move_features
////
struct node_move_features{};

template <typename Graph>
struct property_map<node_move_features, Graph> {
  typedef IterableEditableValueMap<Graph, typename Graph::Node, std::map<unsigned, feature_array> > type;
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
  typedef IterableEditableValueMap<Graph, typename Graph::Node, std::map<std::pair<unsigned, unsigned>, feature_array> > type;
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
  typedef property_map<node_division_features, base_graph>::type DivisionFeatureMap;
  typedef property_map<node_move_features, base_graph>::type MoveFeatureMap;
  typedef property_map<node_count_features, base_graph>::type CountFeatureMap;
  typedef property_map<node_connected_component, base_graph>::type ConnectedComponentMap;


  MultiHypothesesGraph();

  void remove_traxel_features();

  // if classification is done w/o saving, this can be removed
  void add_classifier_features(ClassifierStrategy* move,
                               ClassifierStrategy* division,
                               ClassifierStrategy* count,
                               ClassifierStrategy* detection);

  void add_cardinalities();

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
  boost::shared_ptr<std::map<int, std::vector<std::vector<unsigned> > > > conflicts_node_;
};




////
//// SingleTimestepTraxel_MultiHypothesesBuilder
////
class PGMLINK_EXPORT SingleTimestepTraxel_MultiHypothesesBuilder : public SingleTimestepTraxel_HypothesesBuilder {
public:
  SingleTimestepTraxel_MultiHypothesesBuilder(const TraxelStore* ts, const Options& o);
  virtual MultiHypothesesGraphPtr build_multi_hypotheses_graph() const;
  virtual HypothesesGraph* build() const;
protected:
  virtual HypothesesGraph* construct() const;
  virtual HypothesesGraph* add_nodes(HypothesesGraph* graph) const;
  virtual HypothesesGraph* add_edges(HypothesesGraph* graph) const;

};


////
//// MultiHypothesesTraxelStore
////

struct MultiHypothesesTraxelStore {
  void add(const Traxel& trax, unsigned /*obsolete: component_id*/);
  const Traxel& get(int timestep, unsigned component_id, unsigned traxel_id) const;
  void add_conflict_map(int timestep, const std::map<int, std::vector<std::vector<int> > >& new_conflicts);
  std::string print();
  TraxelStore ts;
  boost::shared_ptr<ConflictMap > conflicts;

 private:
  friend class boost::serialization::access;
  template< typename Archive >
  void serialize( Archive& ar, const unsigned int /*version*/ );
};


////
//// Implementations
////

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
  ar & conflicts_;
  ar & conflicts_node_;

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
  ar & conflicts_;
  ar & conflicts_node_;

  bool with_n_traxel;
  ar & with_n_traxel;

  std::string lgf;
  ar & lgf;
  {
    std::stringstream ss(lgf);
    read_lgf(*this, ss, with_n_traxel);
  }
}


template< typename Archive >
void MultiHypothesesTraxelStore::serialize( Archive& ar, const unsigned int /*version*/ ) {
  ar & ts;
  ar & conflicts;
}





 
/* IMPLEMENTATIONS */  
    
}


#endif /* MULTI_HYPOTHESES_GRAPH */




